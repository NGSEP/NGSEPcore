package ngsep.clustering;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class BestStarMultipleSequenceAlignmentAlgorithm implements MultipleSequenceAlignmentAlgorithm {

	@Override
	public QualifiedSequenceList calculateMultipleSequenceAlignment(QualifiedSequenceList sequences) {
		// TODO Auto-generated method stub
		return null;
	}

	/**
	 * Calculates the distance matrix for a collection of sequences, putting in each entry the calculation
	 * of the minimum edit distance.
	 * @param sequences - The list of sequences
	 * @return The associated distance matrix
	 */
	public double [][] calculatePairwiseEditDistanceMatrix (QualifiedSequenceList sequences) {
		int n = sequences.size();
		// Distance matrix between the sequences
		double [][] D = new double[n][n];
		int i = 0;
		int j = 0;

		for (QualifiedSequence seq1: sequences){
			for (QualifiedSequence seq2: sequences){
				if (i == j){
					D[i][j] = 0.0;
				}else{
					D[i][j] = calculateMinEditDistance(seq1.getCharacters(), seq2.getCharacters());
				}
				j++;
			}
			j = 0;
		}

		return D;
	}


	/**
	 * Calculates the edit distance for two strings of characters
	 * @param s1 - First string
	 * @param s2 - Second string
	 * @return The edit distance matrix
	 */
	private double[][] calculateEditDistanceMatrix(CharSequence s1, CharSequence s2){

		double[][] E = new double[s1.length() + 1][s2.length() + 1];

		for (int i = 0; i < s1.length() + 1; i++) {
			E[i][0] = i;
		}

		for (int i = 0; i < s2.length() + 1; i++) {
			E[0][i] = i;
		}

		for (int i = 1; i < s1.length() + 1; i++) {
			for (int j = 1; j < s2.length() + 1; j++) {
				if (s1.charAt(i - 1) == s2.charAt(j - 1)){
					E[i][j] = E[i - 1][j - 1];
				}else {
					E[i][j] =
							Math.min(E[i][j - 1] + 1, Math.min(E[i - 1][j] + 1, E[i - 1][j - 1] + 1));
				}
			}
		}

		return E;
	}

	/**
	 * Calculates the minimum edit distance between two strings of characters.
	 * @param s1 - First string
	 * @param s2 - Second string
	 * @return minimum edit distance
	 */
	private double calculateMinEditDistance(CharSequence s1, CharSequence s2){
		return calculateEditDistanceMatrix(s1, s2)[s1.length()][s2.length()];
	}

	/**
	 * Gets the index of the sequences that minimizes the sum of its row.
	 * @param D - Distance Matrix
	 * @return The index of the center
	 */
	private int getCenterOfStar (double[][] D){
		int min = -1;
		double minSum = Double.MAX_VALUE;
		for (int i = 0; i < D.length; i++) {
			double rowSum = sumOverRow(D, i);
			if (rowSum < minSum ){
				minSum = rowSum;
				min = i;
			}
		}
		return min;
	}

	/**
	 * Calculates the sum over the i'th row of a distance matrix.
	 * @param D - Distance Matrix
	 * @return sum of the elements of the i'th row
	 */
	private double sumOverRow (double[][] D, int i){
		double sum = 0.0;
		for (int j = 0; j < D.length; j++) {
			sum += D[i][j];
		}
		return sum;
	}

}
