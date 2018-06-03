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
import java.util.Set;
import java.util.Stack;

import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.main.CommandsDescriptor;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersCounter;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;

/**
 * Program to align reads to a refernece genome
 * @author German Andrade
 * @author Jorge Duitama
 */
public class ReadsAligner {

	public static final double DEF_MIN_PROPORTION_KMERS = 0.7;
	static final int SEARCH_KMER_LENGTH = 15;
	private double minProportionKmers = DEF_MIN_PROPORTION_KMERS;

	public static final int MAX_SPACE_BETWEEN_KMERS = 200;

	public static void main(String[] args) throws Exception 
	{
		ReadsAligner instance = new ReadsAligner();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String fMIndexFile = args[i++];
		String readsFile = args[i++];
		String outFile = args[i++];
		try (PrintStream out = new PrintStream(outFile)){
			instance.alignReads(fMIndexFile, readsFile, out);
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
	public void alignReads( String fMIndexFile, String readsFile, PrintStream out) throws IOException {
		FMIndex fMIndex = FMIndex.loadFromBinaries(fMIndexFile);
		int totalReads = 0;
		int readsAligned = 0;
		int uniqueAlignments=0;
		long time = System.currentTimeMillis();
		try (FastqFileReader reader = new FastqFileReader(readsFile)) {
			Iterator<RawRead> it = reader.iterator();
			while(it.hasNext()) {
				RawRead read = it.next();
				
				int numAlns = alignRead(fMIndex, read, out);
				totalReads++;
				if(numAlns>0) readsAligned++;
				if(numAlns==1) uniqueAlignments++;
			}
		}
		System.out.println("Total reads: "+totalReads);
		System.out.println("Reads aligned: "+readsAligned);
		System.out.println("Unique alignments: "+uniqueAlignments);
		System.out.println("Overall alignment rate: "+(100.0*readsAligned/(double)totalReads)+"%");
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		System.out.println("Time: "+seconds+" seconds");

	}



	private int alignRead(FMIndex fMIndex, RawRead read, PrintStream out) {
		List<ReadAlignment> alignments = search(fMIndex, read);
		int i=0;
		for (ReadAlignment aln: alignments) {
			String readSeq = read.getSequenceString();
			String qual = read.getQualityScores();
			if(aln.isNegativeStrand()) {
				readSeq = DNAMaskedSequence.getReverseComplement(readSeq).toString();
				qual = new StringBuilder(qual).reverse().toString();
			}
			if(i>0) aln.setFlags(aln.getFlags()+ReadAlignment.FLAG_SECONDARY);
			out.println(
					//1.query name
					read.getName()+"\t"+

					//2.Flag
					aln.getFlags()+"\t"+

					//3.reference sequence name
					aln.getSequenceName()+"\t"+

					//4.POS
					aln.getFirst()+"\t"+

					//5.MAPQ
					"255\t"+

					//6.CIGAR
					read.getLength()+"M\t"+

					//7. RNEXT
					"*\t"+

					//8. PNEXT
					"0\t"+

					//9. TLEN
					"0\t"+

					//10. SEQ
					readSeq+"\t"+

					//11. QUAL
					qual

					);
			i++;
		}
		return alignments.size();
	}

	public List<ReadAlignment> search (FMIndex fMIndex, RawRead read) {
		return kmerBasedInexactSearchAlgorithm(fMIndex, read);
	}

	public List<ReadAlignment> exactSearch (FMIndex fMIndex, RawRead read) {
		return fMIndex.search(read.getSequenceString());
	}

	/**
	 * First approach to allow inexact search
	 * It iterates the Sequences of the genome and if there is at least MIN_ACCURACY percentage of the kmers
	 * it allow the alignment with the first an the last position of the kmers ocurrence
	 * @return 
	 */
	private List<ReadAlignment> kmerBasedInexactSearchAlgorithm (FMIndex fMIndex, RawRead read) 
	{
		String characters=read.getCharacters().toString();
		CharSequence[] kmers = KmersCounter.extractKmers(characters, SEARCH_KMER_LENGTH, true);

		List<ReadAlignment> finalAlignments =  new ArrayList<>();
		if(kmers==null) return finalAlignments;
		int kmersCount=((characters.length()/SEARCH_KMER_LENGTH)+1);
		

		HashMap<String, List<KmerAlignment>> seqHits = getSequenceHits(fMIndex, read,kmers);

		//Processing part
		KmerAlignmentComparator cmp = KmerAlignmentComparator.getInstance();
		Set<String> keys= seqHits.keySet();

		for (String sequenceName:keys)
		{
			List<KmerAlignment> alns = seqHits.get(sequenceName);
			
			Collections.sort(alns,cmp);
			
			Stack<KmerAlignment> stack = new Stack<KmerAlignment>();
			for (int i = 0; i < alns.size(); i++) 
			{
				KmerAlignment actual=alns.get(i);
								
				if(stack.isEmpty() || isKmerAlignmentConsistent(stack.peek(), actual))
				{
					stack.push(actual);
				}
				else 
				{
					insert(finalAlignments, kmersCount, sequenceName, stack,fMIndex,characters);
					//after save and clear the stack save the new alignment, could be good
					stack.push(actual);
				}
			}
			insert(finalAlignments, kmersCount, sequenceName, stack,fMIndex,characters);
			
		}
		return finalAlignments;
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

	
	private void insert(List<ReadAlignment> finalAlignments, int kmersCount, String sequenceName,Stack<KmerAlignment> stack, FMIndex fMIndex, String characters) 
	{
		double percent = (double) stack.size()/kmersCount;
		if(percent>=minProportionKmers)
		{
			KmerAlignment[] arr=new KmerAlignment[stack.size()];
			stack.toArray(arr);
			int first = arr[0].getReadAlignment().getFirst();
			int last = arr[arr.length-1].getReadAlignment().getLast();
			//Instead of just add the sequence we are going to use smith waterman
			CharSequence sequence = fMIndex.getSequence(sequenceName, first-1, last, arr[0].getReadAlignment().isNegativeStrand());
			String result = smithWatermanLocalAlingMent(characters, sequence.toString());
			
			finalAlignments.add(new ReadAlignment(sequenceName, first, first+result.length(), last-first, arr[0].getReadAlignment().getFlags()));
		}
		stack.clear();
	}

	/**
	 * It basically get the no overlapping kmers of  SEARCH_KMER_LENGTH length in the @param read 
	 * and find each alignment of each kmer using the @param fMIndex, saves the alignments in hashmap
	 * @param kmers
	 * @return HashMap with key SequenceName and value a List of alignments that has the kmer value.
	 */
	private HashMap<String, List<KmerAlignment>> getSequenceHits(FMIndex fMIndex, RawRead read,CharSequence[] kmers) {

		HashMap<String,List<KmerAlignment>> seqHits =  new HashMap<String,List<KmerAlignment>>();

		//Avoid overlaps
		for (int i = 0; i < kmers.length; i+=SEARCH_KMER_LENGTH) 
		{
			//Exit loop if kmers[i] is null
			if(kmers[i]==null) continue;

			String kmer =kmers[i].toString();

			//Where is located the kmer in exact way
			exactKmerSearch(fMIndex, i, kmer, seqHits);
		}
		
		if(read.getLength()%SEARCH_KMER_LENGTH!=0) {
			int l = kmers.length-1;
			exactKmerSearch(fMIndex, l, kmers[l].toString(), seqHits);
		}
		return seqHits;
	}
	private void exactKmerSearch(FMIndex fMIndex, int kmerNumber, String kmer, HashMap<String, List<KmerAlignment>> seqHits) {
		List<ReadAlignment> regions=fMIndex.search(kmer);

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
	
	
	
	
	private String smithWatermanLocalAlingMent(String reference, String sequence) 
	{
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
		int[] r={lowestMatrix.length-1,lowestMatrix[0].length-1};

		//The algorithm keeps going back until reach origin
		while(!(r[0]==0 && r[1]==0))
		{
			//We want the path with lowest cost
			try
			{
				//We check if the lower comes from above
				if( lowestMatrix [r[0]-1] [r[1]] < lowestMatrix [r[0]][r[1]-1] && lowestMatrix [r[0]-1] [r[1]] < lowestMatrix [r[0]-1] [r[1]-1])
				{
					int[] a = { r[0]-1,r[1]};

					//Add - to sequenceAlignment
					sequenceAlignment.push("-");

					//Push next character to stackReference
					referenceAlingment.push(stackReference.pop());


					// update r, the actual position
					r=a;
				}
				//We check if the lower comes from left
				else if( lowestMatrix [r[0]] [r[1]-1] < lowestMatrix [r[0]-1][r[1]] && lowestMatrix [r[0]] [r[1]-1] < lowestMatrix [r[0]-1][r[1]-1])
				{
					int[] a = { r[0],r[1]-1};

					//Add - to referenceAlingment
					referenceAlingment.push("-");

					//Push next character to sequenceAlignment
					sequenceAlignment.push(stackSequence.pop());

					// update r, the actual position
					r=a;
				}
				//Check diagonal (left upper)
				else
				{
					int[] a = { r[0]-1,r[1]-1};

					//Push next character to referenceAlingment
					referenceAlingment.push(stackReference.pop());

					//Push next character to sequenceAlignment
					sequenceAlignment.push(stackSequence.pop());

					// update r, the actual position
					r=a;
				}
			}
			catch (Exception e) 
			{
//				e.printStackTrace();
				// se trat� de llegar a posici�n negativa

				//si hay camino desde arriba
				if(r[0]-1>=0)
				{
					int[] a = { r[0]-1,r[1]};

					//Add - to sequenceAlignment
					sequenceAlignment.push("-");

					//Push next character to referenceAlingment
					referenceAlingment.push(stackReference.pop());

					// update r, the actual position
					r=a;
				}
				//si no debe haber camino por la izquierda
				else
				{
					int[] a = { r[0],r[1]-1};

					//Add - to referenceAlingment
					referenceAlingment.push("-");

					//Push next character to sequenceAlignment
					sequenceAlignment.push(stackSequence.pop());

					// update r, the actual position
					r=a;
				}
			}

		}
		//Words are in stacks just is needed to turn them around
		String p1="";
		String p2="";
		while(!referenceAlingment.isEmpty() && !sequenceAlignment.isEmpty())
		{
			p1+=referenceAlingment.pop();
			p2+=sequenceAlignment.pop();
		}

		//System.out.println(p1);
		//System.out.println(p2);
		return p2;
	}
}
