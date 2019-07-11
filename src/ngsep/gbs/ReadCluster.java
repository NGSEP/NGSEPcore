package ngsep.gbs;

import ngsep.sequences.DNASequence;

public class ReadCluster {
	
	int clusterNumber;
	int totalReads = 0;
	int longestRead = 0;
	int[][] refSeqTable;
	char[] refSeq;
	boolean complete;
	
	public ReadCluster(int clusterNumber, int longestRead)      
	{                                                                 
		this.clusterNumber = clusterNumber;
		this.longestRead = longestRead;
		this.refSeqTable = new int[this.longestRead][DNASequence.BASES_ARRAY.length];
		this.refSeq = new char[this.longestRead];
	}
	
	public int[][] addRead(String s) {
		totalReads++;
		for(int i=0; i<s.length(); i++) {
			if(!DNASequence.isInAlphabeth(s.charAt(i))) {
				continue;
			}
			int j = DNASequence.BASES_STRING.indexOf(s.charAt(i));
			refSeqTable[i][j]++;
		}
		return refSeqTable;
	}
	
	public void calcRefSeq() {
		for(int i = 0; i < this.longestRead; i++) {
			int max = 0;
			for(int j = 0; j < DNASequence.BASES_STRING.length(); j++) {
				int next = refSeqTable[i][j];
				if((max <= next) && (next != 0)){
					refSeq[i] = DNASequence.BASES_STRING.charAt(j);
					max = next;
				}
			}
		}
	}
	
	public char[] getRefSeq() {
		return refSeq;
	}
	
	public int getClusterNumber(){
		return clusterNumber;
	}
	public boolean clusterComplete() {
		return complete;
	}
	
	public int longestRead() {
		return longestRead;
	}
	
	public int getNumberOfTotalReads () {
		return totalReads;
	}
	
}
