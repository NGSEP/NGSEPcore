package ngsep.gbs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.RawRead;

public class ReadCluster {
	public static final int MIN_ALLOWED_OVERLAP = 10;
	public static final int ALLOWED_OVERLAP_MISMATCH = 1;
	public static final double CLUSTER_PERCENTAGE_TO_BUILD_ALIGNMENT = 0.1;
	private int alignmentLength = -1;
	private int alignmentPos;
	private double alignmentConfidence=0;		//alignment confidence below 0.5 might indicate something is off with the cluster
	private int clusterNumber;
	private int totalReads = 0;
	private List<RawRead> reads1 = new ArrayList<>();
	private List<RawRead> reads2 = null;
	private List<RawRead> alignment = new ArrayList<RawRead>();
	private List<String> sampleIds = new ArrayList<>();
	private String consensusSequence = null;
	private Integer breakPosition = null;
	private boolean oddAlignmentCluster = false;
	
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
		//For single end
		if(reads2==null) {
			for(RawRead read:reads1) alignment.add(read);
		} else {
			calculateAlignment();	// calculates position where alignment starts
			// if no alignment was found join with "NNNNN"
			if(this.alignmentPos == -1) {

				this.alignmentLength = calculateMaxLength();
				for(int i=0; i< reads1.size();i++) {
					int numberOfN = this.alignmentLength - (reads1.get(i).getCharacters().length() + reads2.get(i).getCharacters().length());
					if (numberOfN == 0) numberOfN = 1;
					String N = KmerPrefixReadsClusteringAlgorithm.PAIRED_END_READS_SEPARATOR.substring(0, numberOfN);
					String QS = KmerPrefixReadsClusteringAlgorithm.PAIRED_END_READS_QS.substring(0, numberOfN);
					CharSequence alignedReadChars = reads1.get(i).getCharacters()+N+DNAMaskedSequence.getReverseComplement(reads2.get(i).getCharacters());
					String alignedReadQS = reads1.get(i).getQualityScores()+QS+reverseSequence(reads2.get(i).getQualityScores());
					String alignedReadId = reads1.get(i).getName();
					RawRead alignedRead = new RawRead(alignedReadId, alignedReadChars, alignedReadQS);

					alignment.add(alignedRead);
				}
			} else {
				for(int i=0; i< reads1.size();i++) {
					CharSequence alignedReadChars = calculateAlignedReadSeq(reads1.get(i).getSequenceString(), reads2.get(i).getSequenceString());
					// Right now, the quality scores are not calculated. Instead, arbitrarily, the quality scores of the reverse read are used. 
					String alignedReadQS = calculateAlignedReadQS(reads1.get(i).getQualityScores(), reads2.get(i).getQualityScores());
					String alignedReadId = reads1.get(i).getName();
					RawRead alignedRead = new RawRead(alignedReadId, alignedReadChars, alignedReadQS);
					alignment.add(alignedRead);
				}
			}
			calculateConsensus();
		}
		
	}
	
	private String calculateAlignedReadQS(String forwardQS, String reverseQS) {
		String alignedReadQS;
		if(this.alignmentLength == -1) {
			throw new RuntimeException("Alignment position has not been calculated for this cluster. CalculateAlignment() must be run first");
		}
		alignedReadQS = forwardQS.substring(0, this.alignmentPos) + reverseSequence(reverseQS);
		return alignedReadQS;
	}
	
	private CharSequence calculateAlignedReadSeq(String forward, String reverse) {
		CharSequence alignedReadSeq;
		if(this.alignmentLength == -1) {
			throw new RuntimeException("Alignment position has not been calculated for this cluster. CalculateAlignment() must be run first");
		}
		alignedReadSeq = forward.substring(0, this.alignmentPos) + DNAMaskedSequence.getReverseComplement(reverse);
		return alignedReadSeq;
	}
	
	/*
	 * Returns true if oddities where found when internally aligning the cluster
	 * see calculateAlignment()
	 */
	public boolean getClusterStatus() {
		return oddAlignmentCluster;
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
	
	/*
	 * This method calculates the consensus position where the forward and the reverse reads start overlapping.
	 * If no consensus is found, the cluster is flagged as an odd cluster
	 */
	public void calculateAlignment() {
		Random random = new Random();
		int numberOfReadsToAlign = (int)(reads1.size() * CLUSTER_PERCENTAGE_TO_BUILD_ALIGNMENT);
		if(numberOfReadsToAlign < KmerPrefixReadsClusteringAlgorithm.getMinClusterDepth()) numberOfReadsToAlign = KmerPrefixReadsClusteringAlgorithm.getMinClusterDepth();
		
		int[] alignmentsPos = new int[numberOfReadsToAlign];
		int[] reverseLengths = new int[numberOfReadsToAlign];
		for(int i = 0; i < numberOfReadsToAlign; i++) {
			int chosenRead = random.nextInt(reads1.size()); 
			int alignmentPos = calculateAlignmentPos(reads1.get(chosenRead), reads2.get(chosenRead));
			reverseLengths[i] = reads2.get(chosenRead).getLength();
			alignmentsPos[i] = alignmentPos;
		}
		
		int consensusReverseLength = reverseLengths[0];
		int consensusAlignPos = alignmentsPos[0];
		int alignPosConfidence = 1;
		int altReverseLength = 0;
		int altAlignPosConfidence = 0;
		int alternativeAlgnPos = -2;
		// allows for one disagreement in alignment, otherwise flags it as an odd cluster.
		for(int i=1; i<numberOfReadsToAlign;i++) {
			if(alignmentsPos[i]==consensusAlignPos) {
				alignPosConfidence++;
			} else if (alternativeAlgnPos==-2) {
				alternativeAlgnPos = alignmentsPos[i];
				altReverseLength = reverseLengths[i];
				altAlignPosConfidence = 1;
			} else if (alternativeAlgnPos == alignmentsPos[i]) {
				altAlignPosConfidence++;
			} else {
				this.alignmentPos = -1;
				oddAlignmentCluster = true;
				break;
			}
		}
		
		if(this.alignmentPos != -1) {
			if(alignPosConfidence <= altAlignPosConfidence) {
				this.alignmentPos = alternativeAlgnPos;
				this.alignmentLength = this.alignmentPos + altReverseLength;
				this.alignmentConfidence = (double) altAlignPosConfidence / numberOfReadsToAlign;
			} else {
				this.alignmentPos = consensusAlignPos;
				this.alignmentLength = this.alignmentPos + consensusReverseLength;
				this.alignmentConfidence = (double) alignPosConfidence / numberOfReadsToAlign;
			}
		} else {
			// Even if no alignment was found, we want to calculate the alignment confidence
			if(alignPosConfidence <= altAlignPosConfidence) {
				this.alignmentConfidence = (double) altAlignPosConfidence / numberOfReadsToAlign;
			} else {
				this.alignmentConfidence = (double) alignPosConfidence / numberOfReadsToAlign;
			}
		}
	}
	
	public List<RawRead> getAlignedReads() {
		
		return alignment;
	}
	private String calculateQualityScores(int i) {
		return null;
	}
	
	/*
	 * Returns the index on the forward read where overlaps starts. If no overlap is found, or the overlap is 
	 * smaller than min allowed overlap, returns -1. 
	 */
	private int calculateAlignmentPos(RawRead forward, RawRead reverse) {
		int alignment = -1;
		int forwardLength = forward.getLength();
		int reverseLength = reverse.getLength();
		int minOverlapLength;
//		CharSequence reverseCharacters = reverseSequence(reverse.getCharacters());
		CharSequence reverseCharacters = DNAMaskedSequence.getReverseComplement(reverse.getCharacters());
		CharSequence reverseQualityScore = reverseSequence(reverse.getQualityScores());
		if(forwardLength <= reverseLength) {
			minOverlapLength = forwardLength;
		} else {
			minOverlapLength = reverseLength;
		}
		int overlapScore;
		for(int i=minOverlapLength; i > MIN_ALLOWED_OVERLAP; i--) {
			CharSequence suffix = forward.getCharacters().subSequence(forwardLength-i, forwardLength);
			CharSequence prefix = reverseCharacters.subSequence(0, i);
			if(suffix.length() != prefix.length()) {
				throw new RuntimeException("Suffix and prefix must be same length");
			}
			overlapScore = calculateOverlapScore(suffix, prefix);
			if(suffix.length() - overlapScore <= ALLOWED_OVERLAP_MISMATCH) {
				// Accept alignment
				alignment = forwardLength-i;
				return alignment;
			}
		}
		return alignment;
	}
	
	/*
	 * This method takes a sequence and returns the reversed sequence
	 */
	private CharSequence reverseSequence(CharSequence seq) {
		char[] reverseArray = new char[seq.length()];
		for(int i=0; i < seq.length(); i++) {
			reverseArray[i] = seq.charAt(seq.length() - i - 1);
		}
		CharSequence reverseSequence = new String(reverseArray);
		return reverseSequence;
	}
	
	/*
	 * Given a suffix and a prefix of the same length, returns the overlap score, where every match
	 * contributes +1 to the total score. 
	 */
	private int calculateOverlapScore(CharSequence suffix, CharSequence prefix) {
		int overlapScore = 0;
		int numberMismatches = 0;
		for(int i = 0; i < suffix.length(); i++) {
			if(numberMismatches > 1) {
				return overlapScore;
			}
			if(suffix.charAt(i) == prefix.charAt(i)) {
				overlapScore++;
			} else {
				numberMismatches++;
			}
		}
		
		return overlapScore;
	}
	
	private List<RawRead> mergePairedEndReads() {
		List<RawRead> mergedReads = null;
		
		return mergedReads;
	}
	
	private List<RawRead> adjustReadsLength(List<RawRead> reads) {
		int maxLength = 0;
		for (RawRead read : reads) {
			int l = read.getSequenceString().length();
			if (maxLength < l) maxLength = l;
		}
		
		List<RawRead> newList = new ArrayList<RawRead>();
		for (RawRead read: reads) {
			String s = read.getSequenceString();
			String q = read.getQualityScores();
			int insertPos = s.indexOf(KmerPrefixReadsClusteringAlgorithm.PAIRED_END_READS_SEPARATOR);
			char scoreChar = KmerPrefixReadsClusteringAlgorithm.PAIRED_END_READS_QS.charAt(0);
			char seqChar = KmerPrefixReadsClusteringAlgorithm.PAIRED_END_READS_SEPARATOR.charAt(0);
			
			String seq = "";
			String score = "";
			for(int k = s.length(); k < maxLength; k++) {
				seq += seqChar;
				score += scoreChar;
			}
			
			s = String.format("%s%s%s", s.substring(0, insertPos), seq, s.substring(insertPos));
			q = String.format("%s%s%s", q.substring(0, insertPos), score, q.substring(insertPos));
			RawRead newRead = new RawRead(read.getName(), s, q);
			newList.add(newRead);
		}
		
		return newList;
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
		int[][] refSeqTable = new int[longestRead][DNASequence.BASES_ARRAY.length];
		
		for(RawRead seq:alignment) {
			CharSequence sequence = seq.getCharacters();
			for(int i=0; i<sequence.length(); i++) {
				char c = sequence.charAt(i);
				if(!DNASequence.isInAlphabeth(c)) {
					continue;
				}
				int j = DNASequence.BASES_STRING.indexOf(c);
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
			if(depth > 0) {
				if((prevDepth / depth) > 2) {
					this.breakPosition = i;
				}
			}
			prevDepth = depth;
			
			if(maxIdx>=0) consensusC[i] = DNASequence.BASES_STRING.charAt(maxIdx);
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
