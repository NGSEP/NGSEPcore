package ngsep.gbs;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.sequences.DNASequence;
import ngsep.sequences.RawRead;
import ngsep.sequencing.ReadsDemultiplex;

public class ReadCluster {
	private Logger log = Logger.getLogger(ReadsDemultiplex.class.getName());
	private int clusterNumber;
	private int totalReads = 0;
	private int longestRead = 0;
	private int[][] refSeqTable;
	private char[] refCharSeq;
	private boolean refSeqDone;
	private List<RawRead> reads = new ArrayList<>();
	private List<String> sampleIds = new ArrayList<>();
	String refSeq = "";
	
	public ReadCluster(int clusterNumber)      
	{                                                                 
		this.clusterNumber = clusterNumber;
		
	}
	
	public double getAvgHammingDist() {
		if(reads.size() == 0) {
			return 0;
		}
		double totHammingDist = 0;
		for(RawRead read:reads) {
			totHammingDist += hammingDist(read.getSequenceString(), getRefSeq());
		}
		return totHammingDist / reads.size();
	}
	
	public void addRead(RawRead read, String sampleId) {
		reads.add(read);
		sampleIds.add(sampleId);
		if(read.getLength() >= longestRead) {
			longestRead = read.getLength();
		}
		totalReads++;
	}
	
	private void calcRefSeq() {
		this.refSeqTable = new int[this.longestRead][DNASequence.BASES_ARRAY.length];
		this.refCharSeq = new char[this.longestRead];
		for(RawRead read:reads) {
			String s = read.getSequenceString();
			if(longestRead<s.length()) longestRead = s.length();
			for(int i=0; i<s.length(); i++) {
				if(!DNASequence.isInAlphabeth(s.charAt(i))) {
					continue;
				}
				int j = DNASequence.BASES_STRING.indexOf(s.charAt(i));
				refSeqTable[i][j]++;
			}
		}
		for(int i = 0; i < this.longestRead; i++) {
			int max = 0;
			for(int j = 0; j < DNASequence.BASES_STRING.length(); j++) {
				int next = refSeqTable[i][j];
				if((max <= next) && (next != 0)){
					this.refCharSeq[i] = DNASequence.BASES_STRING.charAt(j);
					max = next;
				}
			} this.refSeq += this.refCharSeq[i];
		}
		refSeqDone = true;
	}
	
	public String getRefSeq() {
		if(refSeqDone) {
			return refSeq;
		} else {
			calcRefSeq();
			return refSeq;
		}

	}

	private int hammingDist(String str1, String str2) {
	    int count = 0;
	    for(int i = 0; i < str1.length(); i++)
	    {
	        if (str1.charAt(i) != str2.charAt(i))
	            count++;
	    }
	    return count;
	}
	
	public int getClusterNumber(){
		return clusterNumber;
	}
	
	public int longestRead() {
		return longestRead;
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
