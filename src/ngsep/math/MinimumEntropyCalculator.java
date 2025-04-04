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

public class MinimumEntropyCalculator implements EntropyCalculator {
    
    // Value of log_10(2) to avoid redundant calculations
    private static final double LOG2_BASE10 = Math.log10(2);
    // Maximum entropy value for an alphabet of size n
    private final double maxEntropy;

    /**
     * Class constructor
     * 
     * @param alphabetSize Size of the alphabet to be used in the calculation.
     * Set to 0 if the information is unavailable or if computing the maximum entropy is 
     * not relevant
     */
    public MinimumEntropyCalculator(int alphabetSize) {
        maxEntropy = (alphabetSize == 0) ? 0 : (Math.log10(alphabetSize) / LOG2_BASE10);
    }

    /**
     * Method that calculates the minimum entropy for a given sequence
     * 
     * @param sequence Input sequence used to compute the minimum entropy
     * @return Minimum entropy for the given sequence
     */
    public double calculateEntropy(CharSequence sequence) {
        HashMap<Character, Integer> charFrequencies = new HashMap<Character, Integer>();
        int n = sequence.length();
        if (n == 0) return 0;
        // Step to compute the ocurrence frequence for each character in the sequence
        for(int i = 0; i < n; i++) {
            charFrequencies.compute(sequence.charAt(i), (k,v) -> (v == null)? 1 : v + 1);
        }
        double max = 0d;
        // Step to find the maximum probability
        for(int count : charFrequencies.values()) {
            double probability = (double) count / n;
            max = Math.max(max, probability);
        }
         // Computes the minimum entropy for the given sequence
        double entropy = -1d * (Math.log10(max) / LOG2_BASE10);
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
