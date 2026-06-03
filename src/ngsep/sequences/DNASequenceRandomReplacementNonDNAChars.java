package ngsep.sequences;

import java.util.Random;

public class DNASequenceRandomReplacementNonDNAChars extends DNASequence {

	private static final long serialVersionUID = 1L;
	private static Random random = new Random();
	
	public DNASequenceRandomReplacementNonDNAChars (CharSequence sequence) {
		super(sequence);
	}
	/**
	 * Returns the number corresponding with a suitable size substring of the given sequence 
	 * @param seq Sequence to calculate hash
	 * @param start Zero based start position
	 * @param end Zero based end position
	 * @param targetSeq AbstractLimitedSequence having the target alphabet
	 * @return long Positive number representing the substring of seq between start (included) and end (not included)
	 */
	@Override
	public long getLongCode(CharSequence seq, int start, int end) {
		try {
			return super.getLongCode(seq, start, end);
		} catch (IllegalArgumentException e) {
			StringBuilder segment = new StringBuilder(seq.subSequence(start, end).toString().toUpperCase());
			for(int i=0;i<segment.length();i++) {
				if(!isInAlphabet(segment.charAt(i))) {
					int idxR = random.nextInt(4);
					segment.setCharAt(i, DNASequence.BASES_STRING.charAt(idxR));
				}
			}
			return super.getLongCode(segment.toString(), 0, segment.length());
		}
	}
	public static void main(String[] args) {
		String seq = "XAGCCACGGTACattACAGCACCCAAAACACCGCTCGTACNTACSCTCCCACYCACAK";
		DNASequenceRandomReplacementNonDNAChars instance = new DNASequenceRandomReplacementNonDNAChars(seq);
		System.out.println(seq);
		System.out.println(instance.toString());
	}
}
