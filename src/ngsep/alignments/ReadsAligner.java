/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.alignments;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Logger;

import ngsep.alignments.io.ReadAlignmentFileWriter;
import ngsep.genome.ReferenceGenomeFMIndex;
import ngsep.main.CommandsDescriptor;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;

/**
 * Program to align reads to a reference genome
 * @author German Andrade
 * @author Jorge Duitama
 */
public class ReadsAligner {

	private Logger log = Logger.getLogger(ReadsAligner.class.getName());
	public static final double DEF_MIN_PROPORTION_KMERS = 0.7;
	static final int SEARCH_KMER_LENGTH = 15;
	private double minProportionKmers = DEF_MIN_PROPORTION_KMERS;
	private boolean onlyPositiveStrand = false;
	
	private ReferenceGenomeFMIndex fMIndex;

	public static final int MAX_SPACE_BETWEEN_KMERS = 200;

	public static void main(String[] args) throws Exception 
	{
		ReadsAligner instance = new ReadsAligner();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String fMIndexFile = args[i++];
		String readsFile = args[i++];
		String outFile = args[i++];
		instance.fMIndex = ReferenceGenomeFMIndex.loadFromBinaries(fMIndexFile);
		QualifiedSequenceList sequences = instance.fMIndex.getSequencesMetadata();
		
		try (PrintStream out = new PrintStream(outFile);
			ReadAlignmentFileWriter writer = new ReadAlignmentFileWriter(sequences, out)){
			instance.alignReads(readsFile, writer);
		}
	}
	
	/**
	 * @return the minProportionKmers
	 */
	public double getMinProportionKmers() {
		return minProportionKmers;
	}
	
	/**
	 * @param minProportionKmers the minProportionKmers to set
	 */
	public void setMinProportionKmers(double minProportionKmers) {
		this.minProportionKmers = minProportionKmers;
	}
	
	/**
	 * @param minProportionKmers the minProportionKmers to set
	 */
	public void setMinProportionKmers(Double minProportionKmers) {
		this.setMinProportionKmers(minProportionKmers.doubleValue());
	}

	/**
	 * Aligns readsFile with the fMIndexFile
	 * @param fMIndexFile Binary file with the serialization of an FMIndex
	 * @param readsFile Fastq file with the reads to align
	 * @param out
	 * @throws IOException
	 */
	public void alignReads( String readsFile, ReadAlignmentFileWriter writer) throws IOException {
		 
		
		int totalReads = 0;
		int readsAligned = 0;
		int uniqueAlignments=0;
		long time = System.currentTimeMillis();
		try (FastqFileReader reader = new FastqFileReader(readsFile)) {
			Iterator<RawRead> it = reader.iterator();
			while(it.hasNext()) {
				RawRead read = it.next();
				List<ReadAlignment> alns = alignRead(read);
				System.out.println("Alignments for: "+read.getName()+" "+alns.size());
				for(ReadAlignment aln:alns) writer.write(aln);
				int numAlns = alns.size();
				totalReads++;
				if(numAlns>0) readsAligned++;
				if(numAlns==1) uniqueAlignments++;
			}
		}
		log.info("Total reads: "+totalReads);
		log.info("Reads aligned: "+readsAligned);
		log.info("Unique alignments: "+uniqueAlignments);
		log.info("Overall alignment rate: "+(100.0*readsAligned/(double)totalReads)+"%");
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		log.info("Time: "+seconds+" seconds");

	}



	public List<ReadAlignment> alignRead(RawRead read) {
		List<ReadAlignment> alignments = search(read);
		//TODO: Choose better quality and sort by quality
		int qual = alignments.size()==1?40:0;
		int i=0;
		for (ReadAlignment aln: alignments) {
			aln.setReadName(read.getName());
			aln.setAlignmentQuality((short) qual);
			if(i>0) aln.setFlags(aln.getFlags()+ReadAlignment.FLAG_SECONDARY);
			i++;
		}
		return alignments;
	}

	public List<ReadAlignment> search (RawRead read) {
		return kmerBasedInexactSearchAlgorithm(read);
	}

	public List<ReadAlignment> exactSearch (RawRead read) {
		return fMIndex.search(read.getSequenceString(),!onlyPositiveStrand);
	}

	/**
	 * First approach to allow inexact search
	 * It iterates the Sequences of the genome and if there is at least MIN_ACCURACY percentage of the kmers
	 * it allow the alignment with the first an the last position of the kmers ocurrence
	 * @return 
	 */
	private List<ReadAlignment> kmerBasedInexactSearchAlgorithm (RawRead read) {
		List<ReadAlignment> alns = new ArrayList<>();
		DNAMaskedSequence readSeq = (DNAMaskedSequence)read.getCharacters();
		String qual = read.getQualityScores();
		alns.addAll(kmerBasedInexactSearchAlgorithm(readSeq, qual));
		
		if(!onlyPositiveStrand) {
			readSeq = readSeq.getReverseComplement();
			qual = new StringBuilder(qual).reverse().toString();
			List<ReadAlignment> alnsR = kmerBasedInexactSearchAlgorithm(readSeq, qual);
			for (ReadAlignment aln:alnsR) {
				aln.setFlags(aln.getFlags()+ReadAlignment.FLAG_READ_REVERSE_STRAND);
			}
			alns.addAll(alnsR);
		}
		return alns;
	}
	/**
	 * First approach to allow inexact search
	 * It iterates the Sequences of the genome and if there is at least MIN_ACCURACY percentage of the kmers
	 * it allow the alignment with the first an the last position of the kmers ocurrence
	 * Ony tries to align the given quey in the positive strand
	 * @return List<ReadAlignment>
	 */
	private List<ReadAlignment> kmerBasedInexactSearchAlgorithm (CharSequence query, String qualityScores) 
	{
		List<CharSequence> kmers = selectKmers(query);
		

		List<ReadAlignment> finalAlignments =  new ArrayList<>();
		if(kmers==null) return finalAlignments;
		//System.out.println("Query: "+query.toString()+" kmers: "+kmers.size());
		int kmersCount=kmers.size();
		Map<String, List<KmerAlignment>> seqHits = getSequenceHits(kmers);
		//System.out.println("Sequences: "+seqHits.size());
		
		//Processing part
		KmerAlignmentComparator cmp = KmerAlignmentComparator.getInstance();
		Set<String> keys= seqHits.keySet();

		for (String sequenceName:keys)
		{
			List<KmerAlignment> alns = seqHits.get(sequenceName);
			//System.out.println("KMer alns sequence: "+sequenceName+" : "+alns.size());
			
			Collections.sort(alns,cmp);
			
			Stack<KmerAlignment> stack = new Stack<KmerAlignment>();
			for (int i = 0; i < alns.size(); i++) 
			{
				KmerAlignment actual=alns.get(i);				
				if(!stack.isEmpty() && !isKmerAlignmentConsistent(stack.peek(), actual))
				{
					insert(finalAlignments, kmersCount, sequenceName, stack, query, qualityScores);
					stack.clear();
				}
				stack.push(actual);
			}
			insert(finalAlignments, kmersCount, sequenceName, stack, query, qualityScores);
			
		}
		//System.out.println("Found "+finalAlignments.size()+" alignments for query: "+query.toString());
		return finalAlignments;
	}
	
	private List<CharSequence> selectKmers(CharSequence characters) {
		List<CharSequence> kmers = new ArrayList<>();
		int n = characters.length();
		int lastPos = 0;
		for (int i = 0; i+SEARCH_KMER_LENGTH <= n; i+=SEARCH_KMER_LENGTH) {
			CharSequence kmer = characters.subSequence(i, i+SEARCH_KMER_LENGTH);
			if (DNASequence.isDNA(kmer.toString())) kmers.add(kmer);
			lastPos = i;
		}
		if(n-SEARCH_KMER_LENGTH > lastPos) {
			CharSequence kmer = characters.subSequence(n-SEARCH_KMER_LENGTH, n);
			if (DNASequence.isDNA(kmer.toString())) kmers.add(kmer);
		}
		
		return kmers;
	}

	private boolean isKmerAlignmentConsistent(KmerAlignment topAln, KmerAlignment nextAln) 
	{
		boolean negativeStrand = topAln.isNegativeStrand();
		if(negativeStrand != nextAln.isNegativeStrand()) return false;
		if(negativeStrand) {
			if(nextAln.getKmerNumber()>=topAln.getKmerNumber()) return false;
		} else {
			if(nextAln.getKmerNumber()<=topAln.getKmerNumber()) return false;
		}
		if(nextAln.getFirst()-topAln.getFirst()>MAX_SPACE_BETWEEN_KMERS) return false;
		return true;
	}

	
	private void insert(List<ReadAlignment> finalAlignments, int kmersCount, String sequenceName,Stack<KmerAlignment> stack, CharSequence query, String qualityScores) 
	{
		double prop = (double) stack.size()/kmersCount;
		if(prop<minProportionKmers) return;
		KmerAlignment[] arr=new KmerAlignment[stack.size()];
		stack.toArray(arr);
		int first = arr[0].getReadAlignment().getFirst();
		int last = arr[arr.length-1].getReadAlignment().getLast();
		//Instead of just add the sequence we are going to use smith waterman
		CharSequence refSeq = fMIndex.getSequence(sequenceName, first-1, last);
		if(refSeq == null) return;
		//String result = smithWatermanLocalAlingment(refSeq, query.toString());
		//System.out.println(result);
		ReadAlignment aln = new ReadAlignment(sequenceName, first, first+query.length()-1, query.length(), arr[0].getReadAlignment().getFlags());
		aln.setReadCharacters(query);
		aln.setQualityScores(qualityScores);
		aln.setCigarString(query.length()+"M");
		finalAlignments.add(aln);
		
	}

	/**
	 * It basically get the no overlapping kmers of  SEARCH_KMER_LENGTH length in the @param read 
	 * and find each alignment of each kmer using the @param fMIndex, saves the alignments in hashmap
	 * @param kmers
	 * @return Map with key SequenceName and value a List of alignments that has the kmer value.
	 */
	private Map<String, List<KmerAlignment>> getSequenceHits(List<CharSequence> kmers) {

		Map<String,List<KmerAlignment>> seqHits =  new HashMap<String,List<KmerAlignment>>();		
		for (int i=0;i<kmers.size();i++) 
		{
			String kmer = kmers.get(i).toString();
			//System.out.println("Number: "+i+" Kmer: "+kmer);
			//Where is located the kmer in exact way
			exactKmerSearch(i, kmer, seqHits);
		}
		return seqHits;
	}
	private void exactKmerSearch(int kmerNumber, String kmer, Map<String, List<KmerAlignment>> seqHits) {
		
		List<ReadAlignment> regions=fMIndex.search(kmer);
		//System.out.println("Number: "+kmerNumber+" Kmer: "+kmer+" hits: "+regions.size());
		for(ReadAlignment aln:regions)
		{
			KmerAlignment kmerAlignment = new KmerAlignment(kmerNumber, aln);
			List<KmerAlignment> seqAlns = seqHits.get(aln.getSequenceName());
			if(seqAlns==null) {
				seqAlns = new ArrayList<>();
				seqHits.put(aln.getSequenceName(), seqAlns);
			}
			seqAlns.add(kmerAlignment);
		}
	}
	
	
	
	
	private String smithWatermanLocalAlingment(CharSequence reference, CharSequence sequence) 
	{
//		System.out.println("ref");
//		System.out.println(reference);
//		
//		System.out.println("seq");
//		System.out.println(sequence);
		//Stack containing reference
		Stack<String> stackReference = new Stack<>();

		//Stack containing secuence
		Stack<String> stackSequence = new Stack<>();

		//Matrix with 0 if charAt(i)==charAt(j) otherwise 1
		int[][] diagonals = new int[reference.length()][sequence.length()];


		for (int i = 0; i < diagonals.length; i++) 
		{
			char actualReference=reference.charAt(i);
			stackReference.push(actualReference+"");

			for (int j = 0; j < diagonals[i].length; j++) 
			{
				char actualSequence=sequence.charAt(j);
				if(i==0)
					stackSequence.push(actualSequence+"");

				diagonals[i][j]= actualReference==actualSequence ? 0:1;
			}
		}

		//Matrix for dynamic programming saves the lowest weight from 0,0 to i,j
		int[][] lowestMatrix = new int[diagonals.length+1][diagonals[0].length+1]; 

		for (int i = 0; i < lowestMatrix.length; i++) 
		{	
			for (int j = 0; j < lowestMatrix[0].length; j++) 
			{
				//base case, the cost to reach 0,0 is 0
				if(i==0 && j==0)
				{
					lowestMatrix[i][j]=0;
				}
				//Base case: delete a base
				//the cost to reach 0,j is 1+ cost(0,j-1)
				//This an horizontal move  -->
				else if(i==0)
				{
					lowestMatrix[i][j]= 1 + lowestMatrix[i][j-1];
				}
				//Base case: insert a base
				//the cost to reach i,0 is 1 + cost(i-1,0)
				//This an vertical move |
				//						v	
				else if(j==0)
				{
					lowestMatrix[i][j]= 1 + lowestMatrix[i-1][j];
				}
				//This case can be reached from upper (i-1,j)
				//From left (i,j-1)
				//or from upper left diagonal (i-1,j-1)
				else
				{
					int[] options= {
							1 						+ lowestMatrix[i-1][j], 	// arrive from up is 1 + cost(i-1,j)
							1 						+ lowestMatrix[i][j-1], 	// arrive from left is 1 + cost(i-1,j)
							diagonals[i-1][j-1] 	+ lowestMatrix[i-1][j-1]	// arrive from diagonal is 1 if there is different characters
																				// and 0 if are equals, +cost(i-1,j-1)
					};
					int min = Integer.MAX_VALUE;
					for (int k = 0; k < options.length; k++) 
					{
						if(options[k]<min)
						{
							min = options[k];
						}
					}

					lowestMatrix[i][j]=min;
				}
			}
		}
		//At this point we have the minimum cost to reach the corner (A.length-1,A[0].length-1)
		//Now we have to go back and remember the decisions
		
		//Stack containing the alignment of the reference
		Stack<String> referenceAlingment = new Stack<>();

		//Stack containing the alignment of the sequence
		Stack<String> sequenceAlignment = new Stack<>();

		//Stores the actual position of the go back algorithm.
		//It starts in in the right bottom corner, where is the minimum cost to reach (A.length-1,A[0].length-1)
		int[] actualPositionBackTrack={lowestMatrix.length-1,lowestMatrix[0].length-1};

		//The algorithm keeps going back until reach origin
		while(!(actualPositionBackTrack[0]==0 && actualPositionBackTrack[1]==0))
		{
			//We want the path with lowest cost
			try
			{
				//We check if the lower comes from above
				if( lowestMatrix [actualPositionBackTrack[0]-1] [actualPositionBackTrack[1]] < lowestMatrix [actualPositionBackTrack[0]][actualPositionBackTrack[1]-1] && lowestMatrix [actualPositionBackTrack[0]-1] [actualPositionBackTrack[1]] < lowestMatrix [actualPositionBackTrack[0]-1] [actualPositionBackTrack[1]-1])
				{
					int[] a = { actualPositionBackTrack[0]-1,actualPositionBackTrack[1]};

					//Add - to sequenceAlignment
					
					sequenceAlignment.push("-");

					//Push next character to stackReference
					referenceAlingment.push(stackReference.pop());


					// update r, the actual position
					actualPositionBackTrack=a;
				}
				//We check if the lower comes from left
				else if( lowestMatrix [actualPositionBackTrack[0]] [actualPositionBackTrack[1]-1] < lowestMatrix [actualPositionBackTrack[0]-1][actualPositionBackTrack[1]] && lowestMatrix [actualPositionBackTrack[0]] [actualPositionBackTrack[1]-1] < lowestMatrix [actualPositionBackTrack[0]-1][actualPositionBackTrack[1]-1])
				{
					int[] a = { actualPositionBackTrack[0],actualPositionBackTrack[1]-1};

					//Add - to referenceAlingment
					//referenceAlingment.push("-");

					//Push next character to sequenceAlignment
					sequenceAlignment.push(stackSequence.pop());

					// update r, the actual position
					actualPositionBackTrack=a;
				}
				//Check diagonal (left upper)
				else
				{
					int[] a = { actualPositionBackTrack[0]-1,actualPositionBackTrack[1]-1};

					//Push next character to referenceAlingment
					referenceAlingment.push(stackReference.pop());

					//Push next character to sequenceAlignment
					sequenceAlignment.push(stackSequence.pop());

					// update r, the actual position
					actualPositionBackTrack=a;
				}
			}
			catch (Exception e) 
			{
//				e.printStackTrace();
				//A negative position is reached

				//If there is a path from above
				if(actualPositionBackTrack[0]-1>=0)
				{
					int[] a = { actualPositionBackTrack[0]-1,actualPositionBackTrack[1]};

					//Add - to sequenceAlignment
					
					sequenceAlignment.push("-");

					//Push next character to referenceAlingment
					referenceAlingment.push(stackReference.pop());

					// update r, the actual position
					actualPositionBackTrack=a;
				}
				//There must be a path from left
				else
				{
					int[] a = { actualPositionBackTrack[0],actualPositionBackTrack[1]-1};

					//Add - to referenceAlingment
					//referenceAlingment.push("-");

					//Push next character to sequenceAlignment
					sequenceAlignment.push(stackSequence.pop());

					// update r, the actual position
					actualPositionBackTrack=a;
				}
			}

		}
		//Words are in stacks just is needed to turn them around
		String p1="";
		while(!referenceAlingment.isEmpty() && !sequenceAlignment.isEmpty())
		{
			p1+=referenceAlingment.pop();
		}	
		return p1;
	}
}
