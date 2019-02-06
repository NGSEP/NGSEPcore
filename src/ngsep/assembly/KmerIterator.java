package ngsep.assembly;

import java.util.AbstractMap.SimpleEntry;
import java.util.Iterator;
import java.util.Map.Entry;

@FunctionalInterface
public interface KmerIterator {
	public static final int SEARCH_KMER_LENGTH = 15;
	public static final int SEARCH_KMER_DISTANCE = 0;
	public static final int MAX_KMER_DES = 6;

	public Iterable<Entry<Integer, String>> iter(CharSequence seq);

	public static KmerIterator NORMAL = (CharSequence seq) -> {
		String str = seq.toString();
		return new Iterable<Entry<Integer, String>>() {
			@Override
			public Iterator<Entry<Integer, String>> iterator() {
				return new Iterator<Entry<Integer, String>>() {
					int i = 0;

					@Override
					public boolean hasNext() {
						return i < str.length() - SEARCH_KMER_LENGTH;
					}

					@Override
					public Entry<Integer, String> next() {
						String ans = str.substring(i, i + SEARCH_KMER_LENGTH);
						i += SEARCH_KMER_LENGTH + SEARCH_KMER_DISTANCE;
						return new SimpleEntry<>(i, ans);
					}
				};
			}
		};
	};

	public static KmerIterator REVERSE = (CharSequence seq) -> {
		String str = (new StringBuilder(seq)).reverse().toString();
		return new Iterable<Entry<Integer, String>>() {
			@Override
			public Iterator<Entry<Integer, String>> iterator() {
				return new Iterator<Entry<Integer, String>>() {
					int i = 0;

					@Override
					public boolean hasNext() {
						return i < str.length() - SEARCH_KMER_LENGTH;
					}

					@Override
					public Entry<Integer, String> next() {
						String ans = str.substring(i, i + SEARCH_KMER_LENGTH);
						i += SEARCH_KMER_LENGTH + SEARCH_KMER_DISTANCE;
						return new SimpleEntry<>(i, ans);
					}
				};
			}
		};
	};

}
