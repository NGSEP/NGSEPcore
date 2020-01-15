package ngsep.sequences;

import java.io.PrintStream;
import java.util.Arrays;

import ngsep.math.Distribution;

public class ByteArrayKmersMapImpl implements KmersMap {

	byte kmerLength;
	private byte [] kmerCounts;
	private int size = 0;
	private static final DNASequence dummySequence = new DNASequence();
	
	public ByteArrayKmersMapImpl( byte kmerLength) {
		if(kmerLength>15) throw new IllegalArgumentException("The kmer length must be less or equal than 15");
		this.kmerLength = kmerLength;
		kmerCounts = new byte [(int)Math.pow(4, kmerLength)];
		Arrays.fill(kmerCounts, (byte)0);
	}
	
	@Override
	public int size() {
		return size;
	}

	@Override
	public int getCount(CharSequence kmer) {
		if(kmer.length()!=kmerLength) throw new IllegalArgumentException("Unexpected length for query: "+kmer+" expected: "+kmerLength);
		int hash = (int) AbstractLimitedSequence.getHash(kmer, 0, kmer.length(), dummySequence);
		return kmerCounts[hash];
	}

	public void setCount(CharSequence kmer, int count) {
		if(kmer.length()!=kmerLength) throw new IllegalArgumentException("Unexpected length for query: "+kmer+" expected: "+kmerLength);
		int hash = (int) AbstractLimitedSequence.getHash(kmer, 0, kmer.length(), dummySequence);
		if(count>Byte.MAX_VALUE) count = Byte.MAX_VALUE;
		if(kmerCounts[hash]==0 && count>0) size++;
		kmerCounts[hash] = (byte) count;
	}

	@Override
	public void addOcurrance(CharSequence kmer) {
		if(kmer.length()!=kmerLength) throw new IllegalArgumentException("Unexpected length for query: "+kmer+" expected: "+kmerLength);
		int hash = (int) AbstractLimitedSequence.getHash(kmer, 0, kmer.length(), dummySequence);
		if(kmerCounts[hash]==0) size++;
		if(kmerCounts[hash]<Byte.MAX_VALUE) (kmerCounts[hash])++;
		
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
		Distribution dist = new Distribution(1,Byte.MAX_VALUE, 1);
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
