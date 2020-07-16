package ngsep.gbs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.RawRead;

public class ReadCluster {
	public static final int MIN_ALLOWED_OVERLAP = 10;
	public static final double ALLOWED_OVERLAP_MISMATCH_PROPORTION = 0.1;
	public static final double CLUSTER_PERCENTAGE_TO_BUILD_ALIGNMENT = 0.1;
	private int clusterNumber;
	private String sequenceName;
	private int totalReads = 0;
	private List<RawRead> reads1 = new ArrayList<>();
	private List<RawRead> reads2 = null;
	private List<ReadAlignment> alignment = new ArrayList<ReadAlignment>();
	private List<String> sampleIds = new ArrayList<>();
	private String consensusSequence = null;
	private Integer breakPosition = null;
	
	public ReadCluster(int clusterNumber, boolean pairedEnd)      
	{                                                                 
		this.clusterNumber = clusterNumber;
		sequenceName = ""+clusterNumber;
		if (pairedEnd) reads2 = new ArrayList<RawRead>();
		
	}
	
	public String getSequenceName() {
		return sequenceName;
	}

	/**
	 * Builds the alignment of either single or paired-end reads
	 */
	public void buildAlignment () {
		alignment.clear();
		consensusSequence = null;
		//For single end
		if(reads2==null) {
			calculateConsensusSingleReads();
			for(int i=0;i<reads1.size();i++) {
				RawRead read = reads1.get(i);
				String sampleId = sampleIds.get(i);
				ReadAlignment readAlignment = buildReadAlignment(read, sampleId);
				alignment.add(readAlignment);
			}
		} else {
			calculateAlignmentPaired();	// calculates position where alignment starts
		}	
	}

	private ReadAlignment buildReadAlignment(RawRead read, String sampleId) {
		int readLength = read.getLength();
		ReadAlignment readAlignment = new ReadAlignment(sequenceName, 1, readLength, readLength, 0);
		readAlignment.setQualityScores(read.getQualityScores());
		readAlignment.setReadCharacters(read.getCharacters());
		readAlignment.setReadName(read.getName());
		readAlignment.setCigarString(""+readLength+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH));
		readAlignment.setReadGroup(sampleId);
		return readAlignment;
	}
	
	private void calculateAlignmentPaired() {
		ReadCluster rc1 = new ReadCluster(1, false);
		for(RawRead read: reads1) rc1.addSingleRead(read, "");
		rc1.buildAlignment();
		String consensus1 = rc1.getConsensusSequence();
		ReadCluster rc2 = new ReadCluster(2, false);
		for(RawRead read: reads2) rc2.addSingleRead(read, "");
		rc2.buildAlignment();
		String consensus2 = rc2.getConsensusSequence();
		String revC2 = DNAMaskedSequence.getReverseComplement(consensus2).toString();
		
		
		int pos = calculateAlignmentPos(consensus1, revC2);
		
		if (pos > 0) {
			consensusSequence = consensus1+revC2.substring(pos);
		} else {
			pos = calculateAlignmentPos(revC2, consensus1);
			if(pos>0) consensusSequence = consensus1.substring(0,pos);
			else {
				char [] nStr = new char[5];
				Arrays.fill(nStr, 'N');
				consensusSequence = consensus1+(new String(nStr))+revC2;
			}
		}
		int n = reads1.size();
		for(int i=0;i<n;i++) {
			RawRead r1 = reads1.get(i);
			RawRead r2 = reads2.get(i);
			String sampleId = sampleIds.get(i);
			alignment.addAll(alignPairedEndReadToConsensus(r1,r2,sampleId));			
		}
	}
	
	private List<ReadAlignment> alignPairedEndReadToConsensus(RawRead r1, RawRead r2, String sampleId) {
		List<ReadAlignment> answer = new ArrayList<ReadAlignment>(2);
		int cl = consensusSequence.length();
		int diff1 = cl-r1.getLength();
		if(diff1<=0) {
			//Consensus has been merged.
			DNAMaskedSequence seq1 = new DNAMaskedSequence(r1.getCharacters().subSequence(0, cl));
			String q1 = r1.getQualityScores().substring(0, cl);
			RawRead forward = new RawRead(r1.getName(), seq1, q1); 
			answer.add(buildReadAlignment(forward, sampleId));
			DNAMaskedSequence seq2 = new DNAMaskedSequence(r2.getCharacters());
			String q2 = r2.getQualityScores();
			if(seq2.length()>cl) {
				seq2 = (DNAMaskedSequence)seq2.subSequence(0, cl);
				q2 = q2.substring(0, cl);
			} else if(seq2.length()<cl){
				int numberOfN = cl-seq2.length();
				char [] nStr = new char[numberOfN];
				Arrays.fill(nStr, 'N');
				seq2.append(new String (nStr));
				char [] qStr = new char[numberOfN];
				Arrays.fill(qStr, '!');
				StringBuilder qb2 = new StringBuilder();
				qb2.append(q2);
				qb2.append(qStr);
				q2 = qb2.toString();
			}
			seq2 = seq2.getReverseComplement();
			q2 = reverseSequence(q2);
			RawRead reverse = new RawRead(r2.getName(), seq2, q2); 
			answer.add(buildReadAlignment(reverse, sampleId));
		} else {
			DNAMaskedSequence seq = new DNAMaskedSequence();
			StringBuilder qs = new StringBuilder();
			seq.append(r1.getCharacters());
			qs.append(r1.getQualityScores());
			if (diff1 <= r2.getLength()) {
				CharSequence rs2 = DNAMaskedSequence.getReverseComplement(r2.getCharacters().subSequence(0, diff1));
				seq.append(rs2);
				String qsr = reverseSequence(r2.getQualityScores().substring(0, diff1));
				qs.append(qsr);
			} else {
				int numberOfN = diff1- r2.getLength();
				char [] nStr = new char[numberOfN];
				Arrays.fill(nStr, 'N');
				seq.append(new String (nStr));
				seq.append(DNAMaskedSequence.getReverseComplement(r2.getCharacters()));
				char [] qStr = new char[numberOfN];
				Arrays.fill(qStr, '!');
				qs.append(qStr);
				qs.append(reverseSequence(r2.getQualityScores()));
			}
			RawRead combined = new RawRead(r2.getName(), seq, qs.toString());
			answer.add(buildReadAlignment(combined, sampleId));
		}
		return answer;
	}

	public List<ReadAlignment> getAlignedReads() {
		return alignment;
	}
	
	/**
	 * Returns the index on the second sequence where overlap ends. If no overlap is found, or the overlap is 
	 * smaller than min allowed overlap, returns -1. 
	 */
	private int calculateAlignmentPos(String seq1, String seq2) {
		int s1Length = seq1.length();
		int s2Length = seq2.length();
		int maxOverlapLength;
		
		if(s1Length <= s2Length) {
			maxOverlapLength = s1Length;
		} else {
			maxOverlapLength = s2Length;
		}
		for(int i=maxOverlapLength; i >= MIN_ALLOWED_OVERLAP; i--) {
			CharSequence suffix = seq1.substring(s1Length-i);
			CharSequence prefix = seq2.substring(0, i);
			if(checkOverlap(suffix, prefix)) {
				// Accept alignment
				return i;
			}
		}
		return -1;
	}
	
	/*
	 * This method takes a sequence and returns the reversed sequence
	 */
	private String reverseSequence(String seq) {
		StringBuilder reverse = new StringBuilder();
		reverse.append(seq);
		reverse.reverse();
		return reverse.toString();
	}
	
	/*
	 * Given a suffix and a prefix of the same length, returns the overlap score, where every match
	 * contributes +1 to the total score. 
	 */
	private boolean checkOverlap(CharSequence suffix, CharSequence prefix) {
		int numberMismatches = 0;
		int maxAllowedMismatches = (int)Math.round(ALLOWED_OVERLAP_MISMATCH_PROPORTION*suffix.length());
		for(int i = 0; i < suffix.length(); i++) {
			
			if(suffix.charAt(i) != prefix.charAt(i)) {
				numberMismatches++;
			}
			if(numberMismatches > maxAllowedMismatches) {
				return false;
			}
		}
		return true;
	}
	
	
	/**
	 * Calculates the average of hamming distances between each read and the consensus
	 * @return double averge of hamming distances
	 */
	public double getAverageHammingDistance() {
		if(alignment.size() == 0) {
			return 0;
		}
		double totHammingDist = 0;
		HammingSequenceDistanceMeasure measure = new HammingSequenceDistanceMeasure();
		for(ReadAlignment aln:alignment) {
			CharSequence sequence = aln.getReadCharacters();
			totHammingDist += measure.calculateDistanceDifferentLengths(sequence, consensusSequence);
		}
		return totHammingDist / alignment.size();
	}
	
	public void addSingleRead(RawRead read, String sampleId) {
		if(reads2!=null) throw new RuntimeException("Can not add single reads to a paired-end cluster");
		reads1.add(read);
		sampleIds.add(sampleId);
		totalReads++;
		consensusSequence = null;
	}
	
	public void addPairedEndRead(RawRead read1, RawRead read2, String sampleId) {
		if(reads2==null) throw new RuntimeException("Can not add paired-end reads to a single-end cluster");
		reads1.add(read1);
		reads2.add(read2);
		sampleIds.add(sampleId);
		totalReads++;
		consensusSequence = null;
	}
	
	private void calculateConsensusSingleReads() {
		
		int longestRead = 0;
		for(RawRead seq:reads1) {
			CharSequence sequence = seq.getCharacters();
			if(longestRead<sequence.length()) longestRead = sequence.length();
		}
		int numNucleotides = DNASequence.BASES_STRING.length();
		int[][] refSeqTable = new int[longestRead][numNucleotides+1];
		
		for(RawRead seq:reads1) {
			CharSequence sequence = seq.getCharacters();
			for(int i=0; i<sequence.length(); i++) {
				char c = sequence.charAt(i);
				int j = DNASequence.BASES_STRING.indexOf(c);
				if(j<0) j=numNucleotides;
				refSeqTable[i][j]++;
			}
		}
		char [] consensusC = new char[longestRead];
		int prevDepth = 0;
		Arrays.fill(consensusC, 'N');
		
		for(int i = 0; i < refSeqTable.length; i++) {
			int max = 0;
			int maxIdx = -1;
			int depth = 0;
			for(int j = 0; j < refSeqTable[i].length; j++) {
				int next = refSeqTable[i][j];
				if(next>0 && max <= next){
					maxIdx = j;
					max = next;
				}
				depth = depth + next;
			}
			
			//Report break in cluster structure when depth drops by 50%
			
			if(prevDepth > 2*depth) {
				this.breakPosition = i;
			}
			prevDepth = depth;
			
			if(maxIdx>=0 && maxIdx<numNucleotides) consensusC[i] = DNASequence.BASES_STRING.charAt(maxIdx);
		}
		consensusSequence = new String(consensusC);
	}
	
	/*
	 * Returns the position at which the depth of the cluster changes drastically (by more than 50%)
	 */
	public Integer getBreakPosition() {
		return breakPosition;
	}
	
	public String getConsensusSequence() {
		return consensusSequence;
	}
	
	public int getClusterNumber(){
		return clusterNumber;
	}
	
	public int getNumberOfTotalReads () {
		return totalReads;
	}

	/**
	 * @return the sampleIds
	 */
	public List<String> getSampleIds() {
		return sampleIds;
	}
	
	
}
