package ngsep.sequences;

import java.io.PrintStream;
import java.util.Arrays;

import ngsep.math.Distribution;

public class ShortArrayDNAKmersMapImpl implements KmersMap {

	private byte kmerLength;
	private short [] kmerCounts;
	private int size = 0;
	private static final DNASequence dummySequence = new DNASequence();
	
	public ShortArrayDNAKmersMapImpl( byte kmerLength) {
		if(kmerLength>15) throw new IllegalArgumentException("The kmer length must be less or equal than 15");
		this.kmerLength = kmerLength;
		kmerCounts = new short [(int)Math.pow(4, kmerLength)];
		Arrays.fill(kmerCounts, (short)0);
	}
	
	@Override
	public int size() {
		return size;
	}

	@Override
	public int getCount(CharSequence kmer) {
		if(kmer.length()!=kmerLength) throw new IllegalArgumentException("Unexpected length for query: "+kmer+" expected: "+kmerLength);
		long code = AbstractLimitedSequence.getHash(kmer, 0, kmer.length(), dummySequence);
		return getCount(code);
	}
	
	public int getCount(long kmerCode) {
		int hash = getHashCode(kmerCode);
		return kmerCounts[hash];
	}

	public void setCount(CharSequence kmer, int count) {
		if(kmer.length()!=kmerLength) throw new IllegalArgumentException("Unexpected length for query: "+kmer+" expected: "+kmerLength);
		long code = AbstractLimitedSequence.getHash(kmer, 0, kmer.length(), dummySequence);
		int hash = getHashCode(code);
		if(count>Short.MAX_VALUE) count = Short.MAX_VALUE;
		if(kmerCounts[hash]==0 && count>0) size++;
		kmerCounts[hash] = (short) count;
	}

	@Override
	public void addOcurrance(CharSequence kmer) {
		if(kmer.length()!=kmerLength) throw new IllegalArgumentException("Unexpected length for query: "+kmer+" expected: "+kmerLength);
		long code = AbstractLimitedSequence.getHash(kmer, 0, kmer.length(), dummySequence);
		addCodeOccurance(code);
		
	}
	
	public void addCodeOccurance(long code) {
		// TODO Auto-generated method stub
		int hash = getHashCode(code);
		if(kmerCounts[hash]==0) size++;
		if(kmerCounts[hash]<Short.MAX_VALUE) (kmerCounts[hash])++;
	}
	
	private int getHashCode (long code) {
		if(kmerLength<=15) return (int)code;
		return (int) (code % 1000000000);
	}

	@Override
	public void filterKmers(int minAbundance) {
		for(int i=0;i<kmerCounts.length;i++) {
			if(kmerCounts[i]>0 && kmerCounts[i]<minAbundance) {
				kmerCounts[i]=0;
				size--;
			}
		}
		
	}
	
	

	@Override
	public Distribution calculateAbundancesDistribution() {
		Distribution dist = new Distribution(1, Short.MAX_VALUE, 1);
		for(int i=0;i<kmerCounts.length;i++) {
			if(kmerCounts[i]>0) dist.processDatapoint(kmerCounts[i]);
		}
		return dist;
	}

	@Override
	public void save(PrintStream out) {
		for(int i=0;i<kmerCounts.length;i++) {
			if(kmerCounts[i]>0) {
				char [] sequence = AbstractLimitedSequence.getSequence(i, kmerLength, dummySequence);
				out.println(new String(sequence)+"\t"+kmerCounts[i]);
			}
		}
		
	}

}
