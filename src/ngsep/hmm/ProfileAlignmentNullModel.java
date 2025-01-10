package ngsep.hmm;

public interface ProfileAlignmentNullModel {
	/**
	 * Calculates the score of the given sequence
	 * @param sequence to calculate score
	 * @return Double log of the probability of the sequence under the null model
	 */
	public Double calculateScore(String sequence);
}
