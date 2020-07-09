package ngsep.gbs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.RawRead;

public class ReadCluster {
	public static final int MIN_ALLOWED_OVERLAP = 10;
	public static final int ALLOWED_OVERLAP_MISMATCH = 2;
	public static final double CLUSTER_PERCENTAGE_TO_BUILD_ALIGNMENT = 0.1;
	private int clusterNumber;
	private int totalReads = 0;
	private List<RawRead> reads1 = new ArrayList<>();
	private List<RawRead> reads2 = null;
	private List<RawRead> alignment = new ArrayList<RawRead>();
	private List<String> sampleIds = new ArrayList<>();
	private String consensusSequence = null;
	private Integer breakPosition = null;
	
	public ReadCluster(int clusterNumber, boolean pairedEnd)      
	{                                                                 
		this.clusterNumber = clusterNumber;
		if (pairedEnd) reads2 = new ArrayList<RawRead>();
		
	}
	
	/**
	 * Builds the alignment of either single or paired-end reads
	 */
	public void buildAlignment () {
		alignment.clear();
		consensusSequence = null;
		//For single end
		if(reads2==null) {
			for(RawRead read:reads1) alignment.add(read);
		} else {
			calculateAlignmentPaired();	// calculates position where alignment starts
			// if no alignment was found join with "NNNNN"
			if(consensusSequence==null) {

				int alnLength = calculateMaxLength();
				for(int i=0; i< reads1.size();i++) {
					CharSequence sequence1 = reads1.get(i).getCharacters();
					CharSequence sequence2 = reads2.get(i).getCharacters();
					int numberOfN = alnLength - (sequence1.length() + sequence2.length());
					if (numberOfN == 0) numberOfN = 1;
					char [] nStr = new char[numberOfN];
					Arrays.fill(nStr, 'N');
					char [] qStr = new char[numberOfN];
					Arrays.fill(qStr, '!');
					CharSequence alignedReadChars = sequence1+(new String(nStr))+DNAMaskedSequence.getReverseComplement(sequence2);
					String alignedReadQS = reads1.get(i).getQualityScores()+(new String(qStr))+reverseSequence(reads2.get(i).getQualityScores());
					String alignedReadId = reads1.get(i).getName();
					RawRead alignedRead = new RawRead(alignedReadId, alignedReadChars, alignedReadQS);

					alignment.add(alignedRead);
				}
				calculateConsensus();
			}
		}	
	}
	
	/*
	 * Calculates the minimum length of the joined paired-end reads if there is no alignment
	 */
	private int calculateMaxLength() {
		int maxLength = 0;
		for(int i=0;i<reads1.size();i++) {
			int forwardLength = reads1.get(i).getLength();
			int reverseLength = reads2.get(i).getLength();
			if(maxLength <= forwardLength+reverseLength) {
				maxLength = forwardLength+reverseLength;
			}
		}
		return maxLength;
		
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
		
		if (pos == -1) {
			pos = calculateAlignmentPos(revC2, consensus1);
			if(pos==-1) return;
			consensusSequence = consensus1.substring(0,pos);
		} else {
			consensusSequence = consensus1+revC2.substring(pos);
		}
		int n = reads1.size();
		for(int i=0;i<n;i++) {
			RawRead r1 = reads1.get(i);
			RawRead r2 = reads2.get(i);
			RawRead aligned = alignPairedEndReadToConsensus(r1,r2);
			alignment.add(aligned);			
		}
	}
	
	private RawRead alignPairedEndReadToConsensus(RawRead r1, RawRead r2) {
		DNAMaskedSequence seq = new DNAMaskedSequence();
		StringBuilder qs = new StringBuilder();
		int cl = consensusSequence.length();
		int diff1 = cl-r1.getLength();
		if(diff1<=0) {
			seq.append(r1.getCharacters().subSequence(0, cl));
			qs.append(r1.getQualityScores().subSequence(0, cl));
		} else {
			seq.append(r1.getCharacters());
			qs.append(r1.getQualityScores());
			if (diff1 <= r2.getLength()) {
				CharSequence rs2 = DNAMaskedSequence.getReverseComplement(r2.getCharacters().subSequence(0, diff1));
				seq.append(rs2);
				CharSequence qsr = reverseSequence(r2.getQualityScores().subSequence(0, diff1));
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
		}
		RawRead answer = new RawRead(r1.getName(), seq, qs.toString());
		return answer;
	}

	public List<RawRead> getAlignedReads() {
		
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
	private CharSequence reverseSequence(CharSequence seq) {
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
		for(int i = 0; i < suffix.length(); i++) {
			
			if(suffix.charAt(i) != prefix.charAt(i)) {
				numberMismatches++;
			}
			if(numberMismatches > ALLOWED_OVERLAP_MISMATCH) {
				return false;
			}
		}
		return true;
	}
	
	
	//TODO: Contracts
	public double getAverageHammingDistance() {
		if(alignment.size() == 0) {
			return 0;
		}
		double totHammingDist = 0;
		HammingSequenceDistanceMeasure measure = new HammingSequenceDistanceMeasure();
		for(RawRead seq:alignment) {
			CharSequence sequence = seq.getCharacters();
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
	
	private void calculateConsensus() {
		
		int longestRead = 0;
		for(RawRead seq:alignment) {
			CharSequence sequence = seq.getCharacters();
			if(longestRead<sequence.length()) longestRead = sequence.length();
		}
		int numNucleotides = DNASequence.BASES_STRING.length();
		int[][] refSeqTable = new int[longestRead][numNucleotides+1];
		
		for(RawRead seq:alignment) {
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
		if(consensusSequence==null) {
			calculateConsensus();
		}
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
