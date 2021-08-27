/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.sequences;


/**
 * 
 * @author Jorge Duitama
 *
 */
public class SimpleEditDistanceMeasure implements SequenceDistanceMeasure {

	private int indelDistance = 1;
	private int mismatchDistance = 1;
	@Override
	public double calculateDistance(CharSequence seq1, CharSequence seq2) {
		int [][] scores = calculateScoresMatrix(seq1, seq2);
		return scores[seq1.length()][seq2.length()];
	}
	private int[][] calculateScoresMatrix(CharSequence seq1, CharSequence seq2) {
		int[][] scores = new int[seq1.length()+1][seq2.length()+1];
		
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
		return scores;
	}
	public String [] calculateAlignment(CharSequence seq1, CharSequence seq2) {
		StringBuilder alignedSequence1 = new StringBuilder();
		StringBuilder alignedSequence2 = new StringBuilder();
		int [][] scores = calculateScoresMatrix(seq1, seq2);
		int i = scores.length - 1;
        int j = scores[0].length - 1;

        while (i != 0 || j != 0){
            if(i == 0 && j > 0){
                alignedSequence1.append(LimitedSequence.GAP_CHARACTER);
                alignedSequence2.append(seq2.charAt(j-1));
                j--;
            } else if (j == 0 && i > 0){
                alignedSequence2.append(LimitedSequence.GAP_CHARACTER);
                alignedSequence1.append(seq1.charAt(i-1));
                i--;
            } else {
                int diagonal = scores[i - 1][j - 1];
                boolean equal = seq1.charAt(i-1)==seq2.charAt(j-1);
				if(!equal) diagonal+=mismatchDistance;
                int left = scores[i][j - 1]+indelDistance;
                int up = scores[i - 1][j]+indelDistance;

                int min = Math.min(up, Math.min(left, diagonal));
                

                if(diagonal == min){
                    alignedSequence1.append(seq1.charAt(i-1));
                    alignedSequence2.append(seq2.charAt(j-1));
                    i--;
                    j--;
                } else if (up == min) {
                	alignedSequence1.append(seq1.charAt(i-1));
                	alignedSequence2.append(LimitedSequence.GAP_CHARACTER);
                    i--;
                } else if (left == min) {
                	alignedSequence1.append(LimitedSequence.GAP_CHARACTER);
                	alignedSequence2.append(seq2.charAt(j-1));
                    j--;
                }
            }
        }
		String [] answer = {alignedSequence1.reverse().toString(),alignedSequence2.reverse().toString()};
		return answer;
	}

	@Override
	public double calculateNormalizedDistance(CharSequence seq1, CharSequence seq2) {
		// TODO Implement
		return 0;
	}

}
