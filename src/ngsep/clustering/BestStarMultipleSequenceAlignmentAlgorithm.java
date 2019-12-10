package ngsep.clustering;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.SimpleEditDistanceMeasure;

public class BestStarMultipleSequenceAlignmentAlgorithm implements MultipleSequenceAlignmentAlgorithm {

	/**
	 * Attribute that calculates the minimum edit distance and the pairwise alignment
	 */
	private SimpleEditDistanceMeasure editDistanceMeasure;

	public BestStarMultipleSequenceAlignmentAlgorithm(){
		editDistanceMeasure = new SimpleEditDistanceMeasure();
	}

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
					D[i][j] = editDistanceMeasure.calculateDistance(seq1.getCharacters(), seq2.getCharacters());
				}
				j++;
			}
			i++;
			j = 0;
		}

		return D;
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
