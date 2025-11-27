package ngsep.sequences;

import java.util.Random;

public class AminoacidSequence extends AbstractLimitedSequence {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	/**
	 * Extended aminoacids to facilitate the indexOf operation using ASCII code
	 */
	public static final String AMINOACIDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	
	public static final AminoacidSequence EMPTY_AA_SEQUENCE = new AminoacidSequence();
	/**
	 * Creates a new Aminoacid sequence with the given characters
	 * @param sequence CharSequence with the aminoacids for the new sequence
	 */
	public AminoacidSequence(CharSequence sequence) {
		this.setSequence(sequence);
	}
	/**
	 * Creates a new empty DNA Sequence
	 */
	public AminoacidSequence() {
	}
	@Override
	public String getAlphabet() {
		return AMINOACIDS;
	}
	@Override
	public int getAlphabetSize() {
		return AMINOACIDS.length();
	}
	@Override
	protected int getBitsPerCharacter() {
		return 5;
	}
	
	@Override
	protected int getDefaultIndex() {
		return getAminoacidIndex('X');
	}
	private static int getAminoacidIndex (char aa) {
		int i = aa - 65;
		if(i<0 || i>=AMINOACIDS.length()) {
			return -1;
		}
		return i;
	}
	/**
	 * Test main for methods of AbstractLimitedSequence
	 * @param args
	 */
	public static void main(String[] args) {
		StringBuilder randomSequence = new StringBuilder();
		Random r = new Random();
		int length = r.nextInt(10000)+100000;
		for(int i=0;i<length;i++) {
			int bpI = r.nextInt(20);
			randomSequence.append(AminoacidSequence.AMINOACIDS.charAt(bpI));
		}
		AminoacidSequence aaSeq = new AminoacidSequence(randomSequence);
		if(!randomSequence.toString().equals(aaSeq.toString())) throw new RuntimeException("Sequences not equal");
		long time1 = System.currentTimeMillis();
		long code = AminoacidSequence.EMPTY_AA_SEQUENCE.getLongCode(randomSequence, 0, 6);
		for(int i=6;i<randomSequence.length();i++) {
			code = AminoacidSequence.EMPTY_AA_SEQUENCE.getNextCode(code, 6, randomSequence.charAt(i));
			if(code!= (code & 0x3FFFFFFF)) System.out.println("Code "+code+" larger than 2E30");
			char [] seq2 = AminoacidSequence.EMPTY_AA_SEQUENCE.getSequenceFromCode(code, 6);
			String k1 = new String (seq2);
			String k2 = randomSequence.substring(i-5, i+1);
			if(!k1.equals(k2)) System.out.println("Error encoding / decoding next character. Expected: "+k2+" given: "+k1); 
		}
		long time2 = System.currentTimeMillis();
		System.out.println("Time kmers encoding and decoding: "+ (time2 - time1));
	}
}
