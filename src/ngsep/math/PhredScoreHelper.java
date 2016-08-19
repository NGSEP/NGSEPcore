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
package ngsep.math;
/**
 * Helper class to transfer probabilities in phred scale forth and back
 * @author Jorge Duitama
 */
public class PhredScoreHelper {
	/**
	 * Calculates the phred score given the probability
	 * @param p Probability to calculate
	 * @return byte score=-10*Math.log10(p). 255 if the probability is zero or the score is higher than 255
	 */
	public static short calculatePhredScore(double p) {
		if(p==0) {
			return 255;
		}
		double score = -10*Math.log10(p);
		if(score > 255) {
			return 255;
		}
		return (short)Math.round(score);
	}
	/**
	 * Calculates the probability related with the given score
	 * @param phredScore Score in Phred scale
	 * @return double probability=Math.pow(10.0, (-phredScore/10.0)). Zero if the score is equal to 255
	 */
	public static double calculateProbability(short phredScore) {
		if(phredScore >= 255) {
			return 0; 
		}
		return Math.pow(10.0, -((double)phredScore/10.0));
	}
}
