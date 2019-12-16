package ngsep.clustering;

import ngsep.sequences.LimitedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.SimpleEditDistanceMeasure;

import java.util.ArrayList;
import java.util.List;

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
		QualifiedSequenceList alignedSequences = new QualifiedSequenceList();

		double[][] D = calculatePairwiseEditDistanceMatrix(sequences);
		int minSeqIndex = getCenterOfStar(D);

		// get all pairwise alignments
		List<QualifiedSequenceList> pairwiseAlignments = getAllPairwiseAlignments(sequences, minSeqIndex);

		// add the first pairwise alignment
		alignedSequences.add(pairwiseAlignments.get(0).get(0));
		alignedSequences.add(pairwiseAlignments.get(0).get(1));

		// align pairwise alignments
		for (int i = 1; i < pairwiseAlignments.size(); i++) {

			QualifiedSequence prevCenter = pairwiseAlignments.get(0).get(0);
			QualifiedSequence nextCenter = pairwiseAlignments.get(i).get(0);
			QualifiedSequence pairedSeq = pairwiseAlignments.get(i).get(1);

			// get a new center in accordance to the next pairwise alignment
			CharSequence newCenter = replaceCenter(prevCenter.getCharacters(), nextCenter.getCharacters());

			// Add the sequences to the list
			alignedSequences.get(0).setCharacters(newCenter);
			alignedSequences.add(pairedSeq);

			// If the next center has a gap, force this gap into the rest of the sequences (except for the new paired sequence)
			String nextCenterChars = nextCenter.getCharacters().toString();
			for (int j = 0; j < nextCenterChars.length(); j++) {
				if (nextCenterChars.charAt(j) == LimitedSequence.GAP_CHARACTER){
					for (int k = 1; k < alignedSequences.size() - 1; k++) {
						QualifiedSequence seq = alignedSequences.get(k);
						StringBuilder seqChars = new StringBuilder(seq.getCharacters());
						StringBuilder newSeqChars = seqChars.insert(j, LimitedSequence.GAP_CHARACTER);
						seq.setCharacters(newSeqChars);
					}
				}
			}

			// If there are sequences that do not have the same size as the center, add gaps at the end
			for (int j = 1; j < alignedSequences.size(); j++) {
				QualifiedSequence seq = alignedSequences.get(j);

				if (seq.getLength() < newCenter.length()){
					StringBuilder seqChars = new StringBuilder(seq.getCharacters());
					for (int k = seqChars.length(); k < newCenter.length(); k++) {
						seqChars.append(LimitedSequence.GAP_CHARACTER);
					}
					seq.setCharacters(seqChars);
				}
			}

		}

		return alignedSequences;
	}

	/**
	 * TODO: Make corrections
	 * Gets a new center for the alignment list. if the next center (result from the next pairwise alignment) has
	 * greater or equal size than the previous center, this is the new center, otherwise, the new center is equal
	 * to the next center with the missing characters from the previous center appended to it (such that the length
	 * of the new center is equal to the maximum lenght between the previous and the next).
	 * @param prev - Previous center
	 * @param next - Next center
	 * @return New center
	 */
	private CharSequence replaceCenter(CharSequence prev, CharSequence next){
		StringBuilder newCenter = new StringBuilder(prev);

		for (int i = 0; i < next.length(); i++) {
			if (i > next.length() - 1){
				newCenter.append(next.charAt(i));
			} else if (next.charAt(i) == LimitedSequence.GAP_CHARACTER){
				newCenter.insert(i, LimitedSequence.GAP_CHARACTER);
			}
		}

		return newCenter.toString();
	}

	/**
	 * Gets a list of all pairwise alignments of the sequences with the center of the star given
	 * by its index.
	 * @param sequences - Sequences to align
	 * @param centerIndex - Position of the center of the star in the list of sequences
	 * @return All pairwise alignments
	 */
	private List<QualifiedSequenceList> getAllPairwiseAlignments(QualifiedSequenceList sequences, int centerIndex){
		QualifiedSequence center = sequences.get(centerIndex);
		List<QualifiedSequenceList> pairwiseAlignments = new ArrayList<>();

		for (int i = 0; i < sequences.size(); i++) {
			if (i != centerIndex){
				QualifiedSequence seq = sequences.get(i);
				List<CharSequence> alignedChars = editDistanceMeasure.pairwiseAlignment(center.getCharacters(), seq.getCharacters());
				List<QualifiedSequence> alignedSeqs = new ArrayList<>();
				alignedSeqs.add(new QualifiedSequence(center.getName(), alignedChars.get(0)));
				alignedSeqs.add(new QualifiedSequence(seq.getName(), alignedChars.get(1)));
				pairwiseAlignments.add(new QualifiedSequenceList(alignedSeqs));
			}
		}

		return pairwiseAlignments;
	}

	/**
	 * Calculates the distance matrix for a collection of sequences, putting in each entry the calculation
	 * of the minimum edit distance.
	 * @param sequences - The list of sequences
	 * @return The associated distance matrix
	 */
	private double [][] calculatePairwiseEditDistanceMatrix (QualifiedSequenceList sequences) {
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
