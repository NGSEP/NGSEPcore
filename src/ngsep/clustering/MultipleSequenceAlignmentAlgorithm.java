package ngsep.clustering;

import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.QualifiedSequenceList;

public interface MultipleSequenceAlignmentAlgorithm {
	public QualifiedSequenceList calculateMultipleSequenceAlignment(QualifiedSequenceList sequences);
	
	public static double calculateSumOfPairsScore(QualifiedSequenceList alignedSequences){
        int score = 0;
        int n = alignedSequences.size();
        HammingSequenceDistanceMeasure hamming = new HammingSequenceDistanceMeasure();
        for (int i = 0; i < n; i++) {
        	CharSequence seq1 = alignedSequences.get(i).getCharacters();
            for (int j = i+1; j < n; j++) {
            	CharSequence seq2 = alignedSequences.get(j).getCharacters();
                score += hamming.calculateDistance(seq1, seq2);
            }
        }
        return score;
    }
}
