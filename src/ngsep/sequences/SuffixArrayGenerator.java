package ngsep.sequences;

public interface SuffixArrayGenerator {
	/**
	 * @return the suffix array. The first position is sequence.length
	 */
	public int[] getSuffixArray();
}
