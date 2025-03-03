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
public class CollisionEntropyCalculator implements EntropyCalculator {
    
    private static double[][] SQUARES_CACHE = new double[100][100]; 
    private static double LOG2_BASE10 = Math.log10(2);

    static{
        for(int i = 0; i < 100; i++) {
            int n = i + 1;
            for(int j = 0; j <= i; j++) {
                double count = j + 1;
                SQUARES_CACHE[i][j] = calculateTerm(count, n);
            }
        }
    }

    public static double calculateTerm(double count, int n) {
        double probability = count / n;
        double square = Math.pow(probability, 2d);
        return square;   
    }

    public double calculateEntropy(CharSequence sequence) {
        HashMap<Character, Integer> charFrequencies = new HashMap<Character, Integer>();
        int n = sequence.length();
        if (n == 0) return 0;
        for(int i = 0; i < n; i++) {
            charFrequencies.compute(sequence.charAt(i), (k,v) -> (v == null)? 1 : v + 1);
        }
        double sum = 0d;
        for(int count : charFrequencies.values()) {
			if (n < 100) sum += SQUARES_CACHE[n-1][count-1];
			else sum += calculateTerm(count, n);
		}
        double entropy = -1d * (Math.log10(sum) / LOG2_BASE10);
        return entropy;
    }
}
