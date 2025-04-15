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

import java.util.HashMap;

/**
 * @author Nicolas Rozo Fajardo
 */

public class ShannonEntropyCalculator implements EntropyCalculator {
	
	// Matrix to stores precomputed log inverse probability values
	private static final double[][] LOG_CACHE = new double[100][100];
	// Value of log_10(2) to avoid redundant calculations
	private static final double LOG2_BASE10 = Math.log10(2d);
	// Maximum entropy value for an alphabet of size n
	private final double maxEntropy;

	// For every possible combination of count and n between 1 and 100, the value 
    // of log_2(n / count) is stored. This helps avoid redundant calculations and 
    // improves performance
    static{
        for(int i = 0; i < 100; i++) {
			int n = i + 1;
			for(int j = 0; j <= i; j++) {
				double count = j + 1;
				LOG_CACHE[i][j] = calculateTerm(count, n);
				//System.out.println("Next entropy at "+i+" "+j+":"+LOG_CACHE[i][j] );
			}
		}
    }

	/**
     * Class constructor
     * 
     * @param alphabetSize Size of the alphabet to be used in the calculation.
     * Set to 0 if the information is unavailable or if computing the maximum entropy is 
     * not relevant
     */
	public ShannonEntropyCalculator(int alphabetSize) {
		maxEntropy = (alphabetSize == 0) ? 0 : (Math.log10(alphabetSize) / LOG2_BASE10);
	}

	/**
     * Method that calculates the log_2 of the inverse probability of a character occurring in a 
     * sequence, based on the number of times it appears and the total length of the sequence
     * 
     * @param count Number of appearances or occurrences of the character in the sequence
     * @param n Length of the sequence
     * @return Log_2 of the  inverse probability of occurrence of the character within the sequence
     */
    public static double calculateTerm(double count, int n) {
    	if(count==0) return 0;
        double probability = count / n;
		double inverse = 1 / probability;
		return (probability * Math.log10(inverse) / LOG2_BASE10);
    }

	/**
     * Method that calculates the shannon entropy for a given sequence
     * 
     * @param sequence Input sequence used to compute the shannon entropy
     * @return Shannon entropy for the given sequence
     */
    public double calculateEntropy(CharSequence sequence) {
		HashMap<Character, Integer> charFrequencies = new HashMap<Character, Integer>();
		int n = sequence.length();
		if (n == 0) return 0;
		// Step to compute the ocurrence frequence for each character in the sequence
		for(int i = 0; i < n; i++) {
			charFrequencies.compute(sequence.charAt(i), (k,v) -> (v == null) ? 1 : v+1);
		}
		double entropy = 0d;
		// Step to compute the log_2 inverse probability of each character in the sequence
		// and compute the total entropy
		for(int count : charFrequencies.values()) { 
			if (n < 100) entropy += LOG_CACHE[n-1][count-1];
			else entropy += calculateTerm(count, n);
		}
		return entropy;
	}

	/**
     * Method that retrive the maximum entropy for the given alphabet
     * 
     * @return Maximum entropy for the alphabet
     */
	public double getMaxEntropy() {
        return this.maxEntropy;
    }
}
