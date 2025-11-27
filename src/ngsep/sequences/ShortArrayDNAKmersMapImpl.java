package ngsep.sequences;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.math.Distribution;

public class ShortArrayDNAKmersMapImpl implements KmersMap {

	private byte kmerLength;
	private short [] kmerCounts;
	private int size = 0;
	
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
		long code = DNASequence.EMPTY_DNA_SEQUENCE.getLongCode(kmer, 0, kmer.length());
		return getCount(code);
	}
	
	public int getCount(long kmerCode) {
		int hash = (int)kmerCode;
		if(hash<0 || hash>=kmerCounts.length) return 0;
		return kmerCounts[hash];
	}

	public void setCount(CharSequence kmer, int count) {
		if(kmer.length()!=kmerLength) throw new IllegalArgumentException("Unexpected length for query: "+kmer+" expected: "+kmerLength);
		//Casting is possible because the kmer length is less than 16
		int hash = (int)DNASequence.EMPTY_DNA_SEQUENCE.getLongCode(kmer, 0, kmer.length());
		if(count>Short.MAX_VALUE) count = Short.MAX_VALUE;
		if(kmerCounts[hash]==0 && count>0) size++;
		kmerCounts[hash] = (short) count;
	}

	@Override
	public void addOcurrance(CharSequence kmer) {
		if(kmer.length()!=kmerLength) throw new IllegalArgumentException("Unexpected length for query: "+kmer+" expected: "+kmerLength);
		long code = DNASequence.EMPTY_DNA_SEQUENCE.getLongCode(kmer, 0, kmer.length());
		addCodeOccurance(code);
		
	}
	
	public void addCodeOccurance(long code) {
		// TODO Auto-generated method stub
		int hash = (int)code;
		//synchronized (kmerCounts) {
			if(kmerCounts[hash]==0) size++;
			if(kmerCounts[hash]<Short.MAX_VALUE) (kmerCounts[hash])++;
		//}
		
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
				char [] sequence = DNASequence.EMPTY_DNA_SEQUENCE.getSequenceFromCode(i, kmerLength);
				out.println(new String(sequence)+"\t"+kmerCounts[i]);
			}
		}
		
	}

	@Override
	public List<CharSequence> getKmersWithCount(int count) {
		List<CharSequence> answer = new ArrayList<CharSequence>();
		for(int i=0;i<kmerCounts.length;i++) {
			if(kmerCounts[i]==count) {
				char [] sequence = DNASequence.EMPTY_DNA_SEQUENCE.getSequenceFromCode(i, kmerLength);
				answer.add(new DNAShortKmer(new String(sequence)));
			}
		}
		return answer;
	}
	
	public Map<Integer,Short> getKmerCodesWithCount(int minCount, int maxCount) {
		Map<Integer,Short> answer = new HashMap<Integer,Short>();
		for(int i=0;i<kmerCounts.length;i++) {
			if(kmerCounts[i]>=minCount && kmerCounts[i]<=maxCount) {
				answer.put(i,kmerCounts[i]);
			}
		}
		return answer;
	}

}
