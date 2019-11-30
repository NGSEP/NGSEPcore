package ngsep.alignments;

import java.util.ArrayList;
import java.util.List;

import ngsep.sequences.DNASequence;

public class KmerWithStart {
	private CharSequence kmer;
	private int start;
	public KmerWithStart(CharSequence kmer, int start) {
		super();
		this.kmer = kmer;
		this.start = start;
	}
	/**
	 * @return the kmer
	 */
	public CharSequence getKmer() {
		return kmer;
	}
	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}
	/**
	 * Selects the kmers that will be used to query the given sequence
	 * @param search sequence 
	 * @return List<KmerWithStart> List of kmers to search. Each k-mer includes the start position in the given sequence
	 */
	public static List<KmerWithStart> selectKmers(CharSequence search, int kmerLength, int kmerOffset) {
		List<KmerWithStart> kmers = new ArrayList<>();
		int n = search.length();
		int lastPos = 0;
		for (int i = 0; i+kmerLength <= n; i+=kmerOffset) {
			CharSequence kmer = search.subSequence(i, i+kmerLength);
			if (DNASequence.isDNA(kmer.toString())) kmers.add(new KmerWithStart(kmer, i));
			lastPos = i;
		}
		if(n-kmerLength > lastPos) {
			CharSequence kmer = search.subSequence(n-kmerLength, n);
			if (DNASequence.isDNA(kmer.toString())) kmers.add(new KmerWithStart(kmer, n-kmerLength));
		}

		return kmers;
	}
}