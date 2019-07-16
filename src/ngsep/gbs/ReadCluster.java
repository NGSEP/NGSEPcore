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
	private char[] refSeq;
	private List<RawRead> reads = new ArrayList<>();
	private List<String> sampleIds = new ArrayList<>();
	
	public ReadCluster(int clusterNumber)      
	{                                                                 
		this.clusterNumber = clusterNumber;
		
	}
	
	public void addRead(RawRead read, String sampleId) {
		reads.add(read);
		sampleIds.add(sampleId);
		if(read.getLength() >= longestRead) {
			longestRead = read.getLength();
		}
		totalReads++;
	}
	
	public String getRefSeq() {
		String refSeq = "";
		this.refSeqTable = new int[this.longestRead][DNASequence.BASES_ARRAY.length];
		this.refSeq = new char[this.longestRead];
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
					this.refSeq[i] = DNASequence.BASES_STRING.charAt(j);
					max = next;
				}
			} refSeq += this.refSeq[i];
		}
		return refSeq;
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
	
}
