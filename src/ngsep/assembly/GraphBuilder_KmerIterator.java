package ngsep.assembly;

import java.util.Iterator;
import java.util.AbstractMap.SimpleEntry;
import java.util.Map.Entry;

import ngsep.sequences.DNAMaskedSequence;

public class GraphBuilder_KmerIterator {
	private static final double LN2 = 0.693147806;
	private static final double LN100000 = 11.51292546;

	public int SEARCH_KMER_LENGTH;
	public int SEARCH_KMER_DISTANCE;
	public int REAL_KMER_DISTANCE;
	public int MAX_KMER_DES;

	public GraphBuilder_KmerIterator(double Rate_of_changes, double Rate_of_cuts, double Rate_of_cover) {
		double rate_of_error = Rate_of_changes + Rate_of_cuts - Rate_of_changes * Rate_of_cuts;
		//SEARCH_KMER_LENGTH = (int) (LN2 / rate_of_error);
		SEARCH_KMER_LENGTH = (int) (1 / rate_of_error);
		SEARCH_KMER_DISTANCE = (int) (SEARCH_KMER_LENGTH * ((1 / Rate_of_cover) - 1));
		MAX_KMER_DES = (int) (Rate_of_cuts * (LN100000 / rate_of_error));
		REAL_KMER_DISTANCE = SEARCH_KMER_DISTANCE + SEARCH_KMER_LENGTH;

		System.out.println("		SEARCH_KMER_LENGTH: " + SEARCH_KMER_LENGTH);
		System.out.println("		SEARCH_KMER_DISTANCE: " + SEARCH_KMER_DISTANCE);
		System.out.println("		MAX_KMER_DES: " + MAX_KMER_DES);
	}

	public Iterable<Entry<Integer, String>> positiveStrand(CharSequence seq) {
		String str = seq.toString();
		
		return new Iterable<Entry<Integer, String>>() {
			public Iterator<Entry<Integer, String>> iterator() {
				return new Iterator<Entry<Integer, String>>() {
					int i = 0;

					public boolean hasNext() {
						return i < str.length() - SEARCH_KMER_LENGTH;
					}

					public Entry<Integer, String> next() {
						String str2 = str.substring(i, i + SEARCH_KMER_LENGTH);
						Entry<Integer, String> ans = new SimpleEntry<>(i, str2);
						i += REAL_KMER_DISTANCE;
						return ans;
					}
				};
			};
		};
	}
	
	public Iterable<Entry<Integer, String>> negativeStrand(CharSequence seq) {
		String str = DNAMaskedSequence.getReverseComplement(seq).toString();
		return new Iterable<Entry<Integer, String>>() {
			public Iterator<Entry<Integer, String>> iterator() {
				return new Iterator<Entry<Integer, String>>() {
					int i = 0;

					public boolean hasNext() {
						return i < str.length() - SEARCH_KMER_LENGTH;
					}

					public Entry<Integer, String> next() {
						String str2 = str.substring(i, i + SEARCH_KMER_LENGTH);
						Entry<Integer, String> ans = new SimpleEntry<>(i, str2);
						i += REAL_KMER_DISTANCE;
						return ans;
					}
				};
			};
		};
	}

}