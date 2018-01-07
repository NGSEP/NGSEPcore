package ngsep.sequences;

public class DNAShortKmer implements CharSequence {
	private int index1;
	private int index2=0;
	private byte l1;
	private byte l2=0;
	private static final String ALPHABET = DNASequence.BASES_STRING;
	private static final int ALP_LENGTH = 4;
	
	public DNAShortKmer(CharSequence key) {
		int k = key.length();
		if(k>15) {
			l1 = 15;
			l2 = (byte) (k-15); 
		} else {
			l1 = (byte) k;
		}
		index1 = getHash(key, 0, l1);
		if(l2>0)index2 = getHash(key, l1, l1+l2);
	}
	
	@Override
	public char charAt(int arg0) {
		throw new RuntimeException("Unimplemented method");
	}
	@Override
	public int length() {
		// TODO Auto-generated method stub
		return l1+l2;
	}
	@Override
	public CharSequence subSequence(int arg0, int arg1) {
		throw new RuntimeException("Unimplemented method");
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o) {
		if(!(o instanceof DNAShortKmer)) return false;
		DNAShortKmer kmer2 = (DNAShortKmer) o;
 		return l1== kmer2.l1 && index1==kmer2.index1 && l2==kmer2.l2 && index2==kmer2.index2; 
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		return index1%100000000;
	}
	
	/**
	 * Returns the number corresponding with a suitable size substring of the given sequence 
	 * @param seq Sequence to calculate hash
	 * @param start Zero based first position
	 * @param end Zero based last position
	 * @return int Number representing the substring of seq between start (included) and end (not included)
	 */
	private int getHash(CharSequence seq, int start, int end) {
		int number =0;
		for(int i=start;i<end;i++) {
			number*=ALP_LENGTH;
			int index = ALPHABET.indexOf(seq.charAt(i));
			if(index <0) {
				throw new IllegalArgumentException("Character "+seq.charAt(i)+" not supported by sequence of type "+getClass().getName());
			}
			number+=index;
			if(number<0) throw new RuntimeException("Encoding reached a negative number for sequence: "+seq+" between "+start+" and "+end);
		}
		return number;
	}
	/**
	 * Gets the sequence corresponding with the given hash and the given size
	 * @param number Hash number to decode
	 * @param size Size of the sequence to return
	 * @return char[] Decoded String as a char array
	 */
	private char [] getSequence(int number, int size) {
		char [] answer = new char[size];
		for(int i=0;i<size;i++) {
			int nextDigit = (int)(number%ALP_LENGTH);
			int index = size-i-1;
			answer[index] = ALPHABET.charAt(nextDigit);
			number = number/ALP_LENGTH;
		}
		return answer;
	}
}
