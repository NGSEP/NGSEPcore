package ngsep.assembly;

import java.util.Iterator;

public class KmerEmmbeddedIterator {
	public static final int NUMBER_OF_KMERS = 10;
	public static final int SEARCH_KMER_LENGTH = 15;
	public static final double PORCENTAJE_OF_KMERS = 0.1;
	public static final double PORCENTAJE_OF_TOLERANCE = 0.3;
	public static final double FACTOR = PORCENTAJE_OF_KMERS / (double) SEARCH_KMER_LENGTH;
	private final String s;
	private final int number;

	public KmerEmmbeddedIterator(CharSequence s) {
		this.s = s.toString();
		number = (int) Math.round(s.length() * FACTOR);
	}

	public Iterable<String> firts() {
		return new Iterable<String>() {
			@Override
			public Iterator<String> iterator() {
				return new Iterator<String>() {
					private int i = 0;
					private int c = 0;

					@Override
					public boolean hasNext() {
						return c < number;
					}

					@Override
					public String next() {
						String ans = s.substring(i, i + SEARCH_KMER_LENGTH - 1);
						i += SEARCH_KMER_LENGTH;
						c++;
						return ans;
					}
				};
			}
		};
	}

	public Iterable<String> lasts() {
		return new Iterable<String>() {
			@Override
			public Iterator<String> iterator() {
				return new Iterator<String>() {
					private int i = s.length();
					private int c = 0;

					@Override
					public boolean hasNext() {
						return c < number;
					}

					@Override
					public String next() {
						String ans = s.substring(i - SEARCH_KMER_LENGTH - 1, i);
						i -= SEARCH_KMER_LENGTH;
						c++;
						return ans;
					}
				};
			}
		};
	}
	
	public int getNumber() {
		return number;
	}
}
