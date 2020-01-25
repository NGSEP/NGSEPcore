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
	private List<RawRead> reads = new ArrayList<>();
	private List<String> sampleIds = new ArrayList<>();
	private String consensusSequence = null;
	private Integer breakPosition = null;
	
	public ReadCluster(int clusterNumber)      
	{                                                                 
		this.clusterNumber = clusterNumber;
		
	}
	
	public double getAverageHammingDistance() {
		if(reads.size() == 0) {
			return 0;
		}
		double totHammingDist = 0;
		HammingSequenceDistanceMeasure measure = new HammingSequenceDistanceMeasure();
		for(RawRead read:reads) {
			totHammingDist += measure.calculateDistanceDifferentLengths(read.getSequenceString(), getConsensusSequence());
		}
		return totHammingDist / reads.size();
	}
	
	public void addRead(RawRead read, String sampleId) {
		reads.add(read);
		sampleIds.add(sampleId);
		totalReads++;
		consensusSequence = null;
	}
	
	private void calculateConsensus() {
		int longestRead = 0;
		for(RawRead read:reads) {
			String s = read.getSequenceString();
			if(longestRead<s.length()) longestRead = s.length();
		}
		int[][] refSeqTable = new int[longestRead][DNASequence.BASES_ARRAY.length];
		
		for(RawRead read:reads) {
			String s = read.getSequenceString();
			for(int i=0; i<s.length(); i++) {
				if(!DNASequence.isInAlphabeth(s.charAt(i))) {
					continue;
				}
				int j = DNASequence.BASES_STRING.indexOf(s.charAt(i));
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
	
	public List<RawRead> getReads() {
		return reads;
	}

	/**
	 * @return the sampleIds
	 */
	public List<String> getSampleIds() {
		return sampleIds;
	}
	
	
}
