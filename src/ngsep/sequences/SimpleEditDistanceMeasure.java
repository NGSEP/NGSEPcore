package ngsep.sequences;

public class SimpleEditDistanceMeasure implements SequenceDistanceMeasure {

	@Override
	public double calculateDistance(CharSequence seq1, CharSequence seq2) {
		int[][] scores = new int[seq1.length()+1][seq2.length()+1];
		int indelDistance = 1;
		int mismatchDistance = 1;
		for (int i = 0; i < scores.length; i++) {	
			for (int j = 0; j < scores[i].length; j++) {
				if(i==0 && j==0) {
					scores[i][j]=0;
				} else if (i==0) {
					scores[i][j]=scores[i][j-1]+indelDistance;
				} else if (j==0) {
					scores[i][j]=scores[i-1][j]+indelDistance;
				} else {
					
					scores[i][j]=scores[i-1][j-1];
					boolean equal = seq1.charAt(i-1)==seq2.charAt(j-1);
					if(!equal) scores[i][j]+=mismatchDistance;
					if(scores[i][j-1]+indelDistance < scores[i][j]) scores[i][j] = scores[i][j-1]+indelDistance;
					if(scores[i-1][j]+indelDistance < scores[i][j]) scores[i][j] = scores[i-1][j]+indelDistance;
				}
			}
		}
		return scores[seq1.length()][seq2.length()];
	}

	@Override
	public double calculateNormalizedDistance(CharSequence seq1, CharSequence seq2) {
		// TODO Implement
		return 0;
	}

}
