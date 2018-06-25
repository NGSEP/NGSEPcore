package ngsep.assembly;

import java.util.Iterator;

public class KmerEmmbeddedIterator {
    public final static int NUMBER_OF_KMERS = 4;
    public final static int SEARCH_KMER_LENGTH = 15;
    private final String s;

    public KmerEmmbeddedIterator(CharSequence s) {
	this.s = s.toString();
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
			return c < NUMBER_OF_KMERS;
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
			return c < NUMBER_OF_KMERS;
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
}
