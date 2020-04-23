package ngsep.gbs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ngsep.sequences.DNASequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.RawRead;

public class ReadCluster {
	private int clusterNumber;
	private int totalReads = 0;
	private List<RawRead> reads1 = new ArrayList<>();
	private List<RawRead> reads2 = null;
	private List<CharSequence> alignment = new ArrayList<CharSequence>();
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
		//For single end
		if(reads2==null) {
			for(RawRead read:reads1) alignment.add(read.getCharacters());
		} else {
			//TODO: Paired end reads
			//Build naive initial version using the algorithm concatenating Ns
			//Create draft consensus
			calculateConsensus();
			//Check with consensus if ends could be merged
			//if so, Merge reads and produce final alignment and consensus
		}
		
	}
	public List<RawRead> getAlignedReads() {
		
		return null;
	}
	private String calculateQualityScores(int i) {
		return null;
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
		for(CharSequence sequence:alignment) {
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
		for(CharSequence sequence:alignment) {
			if(longestRead<sequence.length()) longestRead = sequence.length();
		}
		int[][] refSeqTable = new int[longestRead][DNASequence.BASES_ARRAY.length];
		
		for(CharSequence sequence:alignment) {
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
	public List<CharSequence> getAlignment() {
		return alignment;
	}
	
	
}
