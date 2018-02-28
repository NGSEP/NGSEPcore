package ngsep.sequences;

public interface SequenceDistanceMeasure {
	/**
	 * Calculates a measure of distance between two sequences.
	 * The specific measure and its properties is implementation dependent
	 * @param seq1 first sequence
	 * @param seq2 second sequence
	 * @return double distance between the two sequences
	 */
	public double calculateDistance(CharSequence seq1, CharSequence seq2);
	/**
	 * Calculates a measure of distance between two sequences normalized by the length of the sequences.
	 * The specific measure and its properties is implementation dependent
	 * @param seq1 first sequence
	 * @param seq2 second sequence
	 * @return double distance between the two sequences normalized by the length of the sequences
	 */
	public double calculateNormalizedDistance(CharSequence seq1, CharSequence seq2);
}
