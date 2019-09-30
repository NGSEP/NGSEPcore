package ngsep.clustering;

import ngsep.sequences.QualifiedSequenceList;

public class BestStarMultipleSequenceAlignmentAlgorithm implements MultipleSequenceAlignmentAlgorithm {

	@Override
	public QualifiedSequenceList calculateMultipleSequenceAlignment(QualifiedSequenceList sequences) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public double [][] calculatePairwiseEditDistanceMatrix (QualifiedSequenceList sequences) {
		int n = sequences.size();
		double [][] distances = new double[n][n];
		//TODO: Implement
		return distances;
	}

}
