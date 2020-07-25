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
	public static final String MIDDLE_N_SEQUENCE_PAIRED_END = "NNN";
	private int clusterNumber;
	private String sequenceName;
	private int totalReads = 0;
	private List<RawRead> reads1 = new ArrayList<>();
	private List<RawRead> reads2 = null;
	private List<ReadAlignment> alignment = new ArrayList<ReadAlignment>();
	private List<String> sampleIds = new ArrayList<>();
	private String consensusSequence = null;
	private Integer breakPosition1 = null;
	private Integer breakPosition2 = null;
	
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
				ReadAlignment readAlignment = buildUngappedReadAlignment(read, sampleId, false);
				alignment.add(readAlignment);
			}
		} else {
			calculateAlignmentPaired();	// calculates position where alignment starts
		}	
	}

	private ReadAlignment buildUngappedReadAlignment(RawRead read, String sampleId, boolean reverse) {
		String matchCode = ""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		String skipCode = ""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_SKIPFROMREAD);
		int readLength = read.getLength();
		int flag = 0;
		CharSequence seq = read.getCharacters();
		String qs =read.getQualityScores();
		int cl = consensusSequence.length();
		int alnLength = Math.min(readLength, cl);
		int first = 1;
		String cigar = ""+alnLength+matchCode;
		if(reverse) {
			flag=16;
			seq = DNAMaskedSequence.getReverseComplement(seq);
			qs = reverseSequence(qs);
			first = cl-alnLength+1;
			if(readLength>alnLength) cigar = ""+(readLength-alnLength)+skipCode+cigar;
		} else if (readLength>alnLength) cigar += (readLength-alnLength)+skipCode;
		ReadAlignment readAlignment = new ReadAlignment(sequenceName, first, first+alnLength-1, readLength, flag);
		readAlignment.setQualityScores(qs.toString());
		readAlignment.setReadCharacters(seq);
		readAlignment.setReadName(read.getName());
		readAlignment.setCigarString(cigar);
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
				consensusSequence = consensus1+MIDDLE_N_SEQUENCE_PAIRED_END+revC2;
			}
		}
		int n = reads1.size();
		for(int i=0;i<n;i++) {
			RawRead r1 = reads1.get(i);
			RawRead r2 = reads2.get(i);
			String sampleId = sampleIds.get(i);
			alignment.add(buildUngappedReadAlignment(r1, sampleId, false));
			alignment.add(buildUngappedReadAlignment(r2, sampleId, true));			
		}
		if(rc1.breakPosition1!=null && rc1.breakPosition1<consensusSequence.length()) breakPosition1 = rc1.getBreakPosition1();
		if(rc2.getBreakPosition1()!=null && rc2.breakPosition1<consensusSequence.length()) breakPosition2 = consensusSequence.length()-rc2.getBreakPosition1()-1;
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
	 * @return double average of hamming distances
	 */
	public double getAverageHammingDistance() {
		if(alignment.size() == 0) {
			return 0;
		}
		double totHammingDist = 0;
		HammingSequenceDistanceMeasure measure = new HammingSequenceDistanceMeasure();
		for(ReadAlignment aln:alignment) {
			CharSequence sequence = aln.getReadCharacters();
			if(aln.getFirst()==1) {
				totHammingDist += measure.calculateDistanceDifferentLengths(sequence, consensusSequence);
			} else {
				//TODO: Take into account start position
			}	
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
				this.breakPosition1 = i;
			}
			prevDepth = depth;
			
			if(maxIdx>=0 && maxIdx<numNucleotides) consensusC[i] = DNASequence.BASES_STRING.charAt(maxIdx);
		}
		consensusSequence = new String(consensusC);
	}
	
	/**
	 * Returns the position from the consensus start at which the depth of the cluster changes drastically (by more than 50%)
	 * @return Integer
	 */
	public Integer getBreakPosition1() {
		return breakPosition1;
	}
	/**
	 * For paired end consensus, returns the position from the consensus end at which the depth of the cluster changes drastically (by more than 50%)
	 * @return Integer
	 */
	public Integer getBreakPosition2() {
		return breakPosition2;
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
