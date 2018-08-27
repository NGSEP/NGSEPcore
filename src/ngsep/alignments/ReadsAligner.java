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
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.alignments.io.ReadAlignmentFileWriter;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenomeFMIndex;
import ngsep.main.CommandsDescriptor;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
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
			//Load as DNAMaskedSequence to allow reverse complement
			reader.setSequenceType(DNAMaskedSequence.class);
			Iterator<RawRead> it = reader.iterator();
			while(it.hasNext()) {
				RawRead read = it.next();
				List<ReadAlignment> alns = alignRead(read);
				//System.out.println("Alignments for: "+read.getName()+" "+alns.size());
				for(ReadAlignment aln:alns) writer.write(aln);
				if(alns.size()==0) {
					ReadAlignment alnNoMap = new ReadAlignment(null, 0, 0, read.getLength(), ReadAlignment.FLAG_READ_UNMAPPED);
					alnNoMap.setReadName(read.getName());
					alnNoMap.setReadCharacters(read.getCharacters());
					alnNoMap.setQualityScores(read.getQualityScores());
					writer.write(alnNoMap);
				}
				int numAlns = alns.size();
				totalReads++;
				if(numAlns>0) readsAligned++;
				if(numAlns==1) uniqueAlignments++;
				if(totalReads%100000==0) log.info("Processed "+totalReads+" reads. Aligned: "+readsAligned);
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
		List<ReadAlignment> alignments = kmerBasedInexactSearchAlgorithm(read);
		for(ReadAlignment aln:alignments) aln.setReadName(read.getName());
		return filterAlignments(alignments);
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
	 * Only tries to align the given quey in the positive strand
	 * @return List<ReadAlignment>
	 */
	private List<ReadAlignment> kmerBasedInexactSearchAlgorithm (CharSequence query, String qualityScores) 
	{
		List<KmerWithStart> kmers = selectKmers(query);
		List<ReadAlignment> finalAlignments =  new ArrayList<>();
		if(kmers==null) return finalAlignments;
		//System.out.println("Query: "+query.toString()+" kmers: "+kmers.size());
		int kmersCount=kmers.size();
		List<ReadAlignment> initialKmerAlns = searchKmers (kmers);
		Collection<KmerAlignmentCluster> clusteredKmerAlns = clusterKmerAlignments(query, initialKmerAlns); 
		//System.out.println("Clusters: "+clusteredKmerAlns.size());
		

		for (KmerAlignmentCluster cluster:clusteredKmerAlns) {
			ReadAlignment readAln = createNewAlignmentFromConsistentKmers(cluster, kmersCount, query, qualityScores);
			if(readAln!=null) finalAlignments.add(readAln);
		}
		//System.out.println("Found "+finalAlignments.size()+" alignments for query: "+query.toString());
		return finalAlignments;
	}
	
	/**
	 * Selects the kmers that will be used to query the given sequence
	 * @param search sequence 
	 * @return List<KmerWithStart> List of kmers to search. Each k-mer includes the start position in the given sequence
	 */
	private List<KmerWithStart> selectKmers(CharSequence search) {
		List<KmerWithStart> kmers = new ArrayList<>();
		int n = search.length();
		int lastPos = 0;
		for (int i = 0; i+SEARCH_KMER_LENGTH <= n; i+=SEARCH_KMER_LENGTH) {
			CharSequence kmer = search.subSequence(i, i+SEARCH_KMER_LENGTH);
			if (DNASequence.isDNA(kmer.toString())) kmers.add(new KmerWithStart(kmer, i));
			lastPos = i;
		}
		if(n-SEARCH_KMER_LENGTH > lastPos) {
			CharSequence kmer = search.subSequence(n-SEARCH_KMER_LENGTH, n);
			if (DNASequence.isDNA(kmer.toString())) kmers.add(new KmerWithStart(kmer, n-SEARCH_KMER_LENGTH));
		}
		
		return kmers;
	}
	/**
	 * Searches the given kmers in the fmIndex 
	 * @param kmers to search
	 * @return List of alignments of each kmer. The read number of each alignment contains the kmer number.
	 */
	private List<ReadAlignment> searchKmers(List<KmerWithStart> kmers) {
		List<ReadAlignment> answer = new ArrayList<>();
		for (KmerWithStart kmer:kmers) {
			List<ReadAlignment> kmerAlns=fMIndex.search(kmer.getKmer().toString());
			for(ReadAlignment aln:kmerAlns) {
				aln.setReadNumber(kmer.getStart());
				answer.add(aln);
			}
		}
		return answer;
	}
	
	private Collection<KmerAlignmentCluster> clusterKmerAlignments(CharSequence query, List<ReadAlignment> initialKmerAlns) {
		List<KmerAlignmentCluster> clusters = new ArrayList<>();
		QualifiedSequenceList seqs = fMIndex.getSequencesMetadata();
		GenomicRegionSortedCollection<ReadAlignment> alnsCollection = new GenomicRegionSortedCollection<>(fMIndex.getSequencesMetadata());
		alnsCollection.addAll(initialKmerAlns);
		for(QualifiedSequence seq:seqs) {
			clusters.addAll(clusterSequenceKmerAlns(query, alnsCollection.getSequenceRegions(seq.getName())));
		}
		return clusters;
	}
	
	
	private Collection<KmerAlignmentCluster> clusterSequenceKmerAlns(CharSequence query, GenomicRegionSortedCollection<ReadAlignment> sequenceAlns) {
		Collection<KmerAlignmentCluster> answer = new ArrayList<>();
		//System.out.println("Alns to cluster: "+sequenceAlns.size());
		KmerAlignmentCluster cluster=null;
		for(ReadAlignment aln:sequenceAlns) {
			if(cluster==null || !cluster.addAlignment(aln)) {
				cluster = new KmerAlignmentCluster(query, aln);
				answer.add(cluster);
			}
		}
		return answer;
	}

	private ReadAlignment createNewAlignmentFromConsistentKmers(KmerAlignmentCluster cluster, int totalKmers, CharSequence query, String qualityScores) {
		int numDiffKmers = cluster.getNumDifferentKmers();
		double prop = (double) numDiffKmers/totalKmers;
		if(prop<minProportionKmers) return null;
		
		
		String sequenceName = cluster.getSequenceName();
		int first = cluster.getFirst();
		int last = cluster.getLast();
		//System.out.println("Reference region from k-mers: "+first+"-"+last+" all consistent: "+cluster.isAllConsistent()+" firstaln: "+cluster.isFirstAlnPresent()+" last aln: "+cluster.isLastAlnPresent());
		String cigar = query.length()+"M";
		double alnQual = 100.0;
		if (!cluster.isAllConsistent() || !cluster.isFirstAlnPresent() || !cluster.isLastAlnPresent()) {
			//Perform smith waterman
			if(!cluster.isFirstAlnPresent()) first -=10;
			first = Math.max(1, first);
			if(!cluster.isLastAlnPresent()) last+=10;
			last = Math.min(fMIndex.getReferenceLength(sequenceName), cluster.getLast());
			CharSequence refSeq = fMIndex.getSequence(sequenceName, first, last);
			if(refSeq == null) return null;
			//System.out.println(""+query);
			//System.out.println(""+refSeq);
			AlignmentResult result = smithWatermanLocalAlignment(query.toString(),refSeq);
			//TODO: Make better score
			if(result.getDistance()>0.5*query.length()) return null;
			
			//Last must be updated before first
			last = first+result.getSubjectLastIdx();
			first = first + result.getSubjectStartIdx();
			cigar = result.getCigarString();
			alnQual = 100.0* (query.length() - result.getDistance())/query.length();
		}
		ReadAlignment aln = new ReadAlignment(sequenceName, first, last, query.length(), 0);
		aln.setReadCharacters(query);
		aln.setQualityScores(qualityScores);
		//if(aln.getLast()-aln.getFirst()>aln.getReadLength()) System.out.println("Aln for read: "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+"-"+aln.getLast());
		aln.setCigarString(cigar);
		aln.setAlignmentQuality((short) Math.round(alnQual));
		return aln;
	}
	
	private AlignmentResult smithWatermanLocalAlignment(CharSequence query, CharSequence subject) {
		//Matrix for dynamic programming saves the lowest weight from 0,0 to i,j
		int[][] scores = new int[query.length()+1][subject.length()+1]; 
		int minLastRow = scores.length*scores[0].length;
		int minJLastRow = -1;
		for (int i = 0; i < scores.length; i++) {	
			scores[i][0] = i==0?0:scores[i-1][0]+1;
			for (int j = 1; j < scores[0].length; j++) 
			{
				if(i==0)
				{
					scores[i][j]=0;
				}
				//This case can be reached from upper (i-1,j)
				//From left (i,j-1)
				//or from upper left diagonal (i-1,j-1)
				else
				{
					int d = query.charAt(i-1)!=subject.charAt(j-1)?1:0;
					int min = scores[i-1][j-1]+d;
					int a = 1 + scores[i-1][j];
					if(min>a) min= a;
					int b = 1 + scores[i][j-1];
					if(min>b) min=b;
					scores[i][j]=min;
				}
				if(i==scores.length-1 && (j==1 || minLastRow>=scores[i][j])) {
					minLastRow = scores[i][j];
					minJLastRow = j; 
				}
			}
		}
		//At this point we have the minimum cost to reach the corner (A.length-1,A[0].length-1)
		//Now we have to go back and remember the decisions from the minimum of the last row
		int i = scores.length-1;
		int j= minJLastRow;
		//System.out.println("Subject length: "+subject.length()+" minJ last row: "+j+" minscore: "+minLastRow);
		AlignmentResult result = new AlignmentResult();
		result.setSubjectLastIdx(minJLastRow-1);
		result.setDistance(minLastRow);
		//The algorithm keeps going back until reach origin
		while(i>0 && j>0) {
			//We want the path with lowest cost
			int d = query.charAt(i-1)!=subject.charAt(j-1)?1:0;
			int sD = scores[i-1][j-1]+d;
			int sA = 1 + scores[i-1][j];
			if( sD == scores[i][j]) {
				result.addBacktrack(ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH));
				i--;
				j--;
			} else if( sA == scores[i][j]) {
				result.addBacktrack(ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION));
				i--;
			} else {
				result.addBacktrack(ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_DELETION));
				j--;
			}
		}
		result.setSubjectStartIdx(j);
		//System.out.println("Subject start: "+result.getSubjectStartIdx());
		while(i>0) {
			result.addBacktrack(ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION));
			i--;
		}
		return result;
	}
	private List<ReadAlignment> filterAlignments(List<ReadAlignment> alignments) {
		if (alignments.size()==0) return alignments;
		Collections.sort(alignments, new Comparator<ReadAlignment>() {
			@Override
			public int compare(ReadAlignment aln1, ReadAlignment aln2) {
				return aln2.getAlignmentQuality() - aln1.getAlignmentQuality();
			}
		});
		short bestQual = alignments.get(0).getAlignmentQuality();
		//TODO. Investigate alignment score
		double threshold = 0.8*bestQual;
		List<ReadAlignment> filteredAlignments = new ArrayList<>();
		int n = alignments.size();
		for (int i=0;i<n;i++) {
			ReadAlignment aln = alignments.get(i);
			if(aln.getAlignmentQuality()<threshold) continue;
			if(i>0) aln.setFlags(aln.getFlags()+ReadAlignment.FLAG_SECONDARY);
			filteredAlignments.add(aln);
		}
		if(filteredAlignments.size()>1) {
			for(ReadAlignment aln:filteredAlignments) aln.setAlignmentQuality((short) Math.round(0.1*aln.getAlignmentQuality()));
		}
		return filteredAlignments;
	}
}
class AlignmentResult {
	private int subjectStartIdx;
	private int subjectLastIdx;
	private int distance;
	private LinkedList<Character> path = new LinkedList<>();
	
	public void addBacktrack(char decision) {
		path.add(0, decision);
	}
	
	public String getCigarString () {
		StringBuilder cigar = new StringBuilder();
		int nextCount = 0;
		char next = 0;
		for(char c:path) {
			if(c!=next) {
				if(nextCount>0) {
					cigar.append(nextCount);
					cigar.append(next);
				}
				nextCount = 1;
				next = c;
			} else {
				nextCount++;
			}
		}
		if(nextCount>0) {
			cigar.append(nextCount);
			cigar.append(next);
		}
		return cigar.toString();
	}

	/**
	 * @return the subjectStartIdx
	 */
	public int getSubjectStartIdx() {
		return subjectStartIdx;
	}

	/**
	 * @param subjectStartIdx the subjectStartIdx to set
	 */
	public void setSubjectStartIdx(int subjectStartIdx) {
		this.subjectStartIdx = subjectStartIdx;
	}

	/**
	 * @return the subjectLastIdx
	 */
	public int getSubjectLastIdx() {
		return subjectLastIdx;
	}

	/**
	 * @param subjectLastIdx the subjectLastIdx to set
	 */
	public void setSubjectLastIdx(int subjectLastIdx) {
		this.subjectLastIdx = subjectLastIdx;
	}

	/**
	 * @return the distance
	 */
	public int getDistance() {
		return distance;
	}

	/**
	 * @param distance the distance to set
	 */
	public void setDistance(int distance) {
		this.distance = distance;
	}
	
	
}
class KmerWithStart {
	private CharSequence kmer;
	private int start;
	public KmerWithStart(CharSequence kmer, int start) {
		super();
		this.kmer = kmer;
		this.start = start;
	}
	/**
	 * @return the kmer
	 */
	public CharSequence getKmer() {
		return kmer;
	}
	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}	
}
class KmerAlignmentCluster implements GenomicRegion {
	private CharSequence query;
	private List<ReadAlignment> alns=new ArrayList<>();
	private String sequenceName;
	private int first;
	private int last;
	private Set<Integer> kmerNumbers = new HashSet<>();
	private boolean allConsistent = true;
	private boolean repeatedNumber = false;
	private boolean lastAlnPresent = false;
	
	public KmerAlignmentCluster(CharSequence query, ReadAlignment aln) {
		this.query = query;
		sequenceName = aln.getSequenceName();
		int kmerQueryStart = aln.getReadNumber();
		first = aln.getFirst() - kmerQueryStart;
		last = aln.getFirst()+(query.length()-kmerQueryStart-1);
		alns.add(aln);
		kmerNumbers.add(kmerQueryStart);
		lastAlnPresent = kmerQueryStart+aln.length()==query.length();
	}

	@Override
	public String getSequenceName() {
		return sequenceName;
	}

	@Override
	public int getFirst() {
		return first;
	}

	@Override
	public int getLast() {
		return last;
	}

	@Override
	public int length() {
		return last-first+1;
	}

	@Override
	public boolean isPositiveStrand() {
		return true;
	}

	@Override
	public boolean isNegativeStrand() {
		return false;
	}
	public boolean addAlignment(ReadAlignment aln) {
		int kmerQueryStart = aln.getReadNumber();
		int estFirst = aln.getFirst() - kmerQueryStart;
		int estLast = aln.getFirst()+(query.length()-kmerQueryStart-1);
		//System.out.println("Previous coords: "+first+"-"+last+" next cords: "+estFirst+"-"+estLast);
		if(first > estLast || last < estFirst) return false;
		if(first != estFirst) allConsistent = false;
		if(last != estLast) allConsistent = false;
		if(kmerNumbers.contains(kmerQueryStart)) repeatedNumber = true;
		else kmerNumbers.add(kmerQueryStart);
		if(kmerQueryStart+aln.length()==query.length()) lastAlnPresent=true;
		if(first>estFirst) first = estFirst;
		if(last<estLast) last = estLast;
		alns.add(aln);
		return true;	
	}

	/**
	 * @return the query
	 */
	public CharSequence getQuery() {
		return query;
	}

	/**
	 * @return the kmerNumbers
	 */
	public int getNumDifferentKmers() {
		return kmerNumbers.size();
	}

	/**
	 * @return the allConsistent
	 */
	public boolean isAllConsistent() {
		return allConsistent;
	}

	/**
	 * @return the repeatedNumber
	 */
	public boolean isRepeatedNumber() {
		return repeatedNumber;
	}

	public boolean isFirstAlnPresent() {
		return kmerNumbers.contains(0);
	}
	/**
	 * @return the lastAlnPresent
	 */
	public boolean isLastAlnPresent() {
		return lastAlnPresent;
	}
	
	
	
}
