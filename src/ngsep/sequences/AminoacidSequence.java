package ngsep.sequences;

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
}
