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

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ngsep.math.CountsRankHelper;

/**
 * Hamming distance measure implementation 
 * @author Jorge Duitama
 */
public class HammingSequenceDistanceMeasure implements SequenceDistanceMeasure {

	/**
	 * Calculates the hamming distance between the two sequences
	 * PRE: The two sequences have the same length
	 * @param seq1 first sequence
	 * @param seq2 second sequence
	 * @return double hamming distance between the two sequences
	 */
	@Override
	public double calculateDistance(CharSequence seq1, CharSequence seq2) {
		assert seq1.length() == seq2.length();
		int answer = 0;
		int l = seq1.length();
		for(int i=0;i<l;i++) {
			if(seq1.charAt(i)!=seq2.charAt(i)) answer++;
		}
		return answer;
	}
	/**
	 * Calculates the hamming distance between the two sequences up to the length of the smallest sequence
	 * @param seq1 first sequence
	 * @param seq2 second sequence
	 * @return double hamming distance between the two sequences
	 */
	public double calculateDistanceDifferentLengths(CharSequence seq1, CharSequence seq2) {
		int answer = 0;
		int l1 = seq1.length();
		int l2 = seq2.length();
		for(int i=0;i<l1 && i<l2;i++) {
			if(seq1.charAt(i)!=seq2.charAt(i)) answer++;
		}
		return answer;
	}
	
	/**
	 * Calculates the hamming distance between the two sequences divided by the length of the sequences
	 * PRE: The two sequences have the same length
	 * @param seq1 first sequence
	 * @param seq2 second sequence
	 * @return double hamming distance between the two sequences divided by the common length of the sequences.
	 * Returns zero if the length of the sequences is zero 
	 */
	@Override
	public double calculateNormalizedDistance(CharSequence seq1, CharSequence seq2) {
		if(seq1.length()==0) return 0;
		return calculateDistance(seq1, seq2)/seq1.length();
	}
	
	/**
	 * Makes a consensus of the given sequences
	 * PRE: The list is not empty and all sequences have the same length
	 * @param sequences to calculate consensus
	 * @return String representing the consensus of the given sequences
	 */
	public static String makeHammingConsensus(List<? extends CharSequence> sequences) {
		assert sequences.size()>0;
		int l = sequences.get(0).length();
		for(CharSequence seq:sequences) assert seq.length()==l;
		StringBuilder sb = new StringBuilder();
		for(int i=0;i<l;i++) {
			CountsRankHelper<Character> countsChar = new CountsRankHelper<>();
			for(CharSequence seq:sequences) {
				countsChar.add(seq.charAt(i));
			}
			sb.append(countsChar.selectBest(1).keySet().iterator().next());
		}
		return sb.toString();
	}
	public static double[] calculateMinorRelativeFrequencies(List<? extends CharSequence> sequences) {
		assert sequences.size()>0;
		int l = sequences.get(0).length();
		for(CharSequence seq:sequences) assert seq.length()==l;
		double [] answer = new double[l];
		Arrays.fill(answer, 0);
		for(int i=0;i<l;i++) {
			CountsRankHelper<Character> countsChar = new CountsRankHelper<>();
			for(CharSequence seq:sequences) {
				countsChar.add(seq.charAt(i));
			}
			if(countsChar.getNumDifferent()>=2) {
				Map<Character, Integer> twoBest = countsChar.selectBest(2);
				Iterator<Integer> it = twoBest.values().iterator();
				double c1 = it.next();
				double c2 = it.next();
				assert c2 <= c1;
				answer[i] = c2 / (c1+c2); 
			}
		}
		return answer;
	}
}
