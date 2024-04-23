package ngsep.alignments;

public interface PairwiseAligner {
	/**
	 * Calculates a pairwise alignment between two sequences
	 * @param sequence1 to align
	 * @param sequence2 to align
	 * @return PairwiseAlignment objet with the aligned sequences and context information
	 */
	public PairwiseAlignment calculateAlignment (CharSequence sequence1, CharSequence sequence2);
}
