package ngsep.sequences;

public class DNAShortKmer implements CharSequence, Comparable<DNAShortKmer> {
	private long index;
	private byte length;
	private static final DNASequence EMPTYDNASEQ = new DNASequence();
	
	public DNAShortKmer(CharSequence kmerSeq) {
		length = (byte) kmerSeq.length();
		if(length>31) throw new IllegalArgumentException("The maximum k-mer size for this class is 31. Input k-mer: "+kmerSeq);
		index = AbstractLimitedSequence.getHash(kmerSeq, 0, length,EMPTYDNASEQ);
	}
	
	@Override
	public char charAt(int i) {
		char [] characters = AbstractLimitedSequence.getSequence(index, length, EMPTYDNASEQ);
		return characters[i];
	}
	@Override
	public int length() {
		return length;
	}
	@Override
	public CharSequence subSequence(int start, int end) {
		String s = new String (AbstractLimitedSequence.getSequence(index, length, EMPTYDNASEQ));
		return s.subSequence(start, end);
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o) {
		if(!(o instanceof DNAShortKmer)) return false;
		DNAShortKmer kmer2 = (DNAShortKmer) o;
 		return length== kmer2.length && index==kmer2.index; 
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return (int) (index%100000000);
	}

	@Override
	public int compareTo(DNAShortKmer kmer2) {
		if(length!=kmer2.length) return length - kmer2.length;
		if(index == kmer2.index) return 0;
		return (index<kmer2.index)?-1:1;
	}
}
