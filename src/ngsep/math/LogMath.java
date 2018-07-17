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
 * Class with static methods performing basic math operations that receive and
 * return logarithms of the values to operate. Minus infinitum is represented
 * as a null object 
 * @author Jorge Duitama
 *
 */
public class LogMath {
	public static final double MAXLOGDIFF=20;
	/**
	 * Null aware sum of probabilities, also scalable to small values.
	 * The sum is calculated as p+q = log(p)+log(1+exp(log(q)-log(p)))
	 * @param log1 10-based logarithm of the first probability to add
	 * @param log2 10-based logarithm of the second probability to add
	 * @return Double logarithm of the sum of the probabilities. Null if both parameters are null (0+0=0)
	 */
	public static Double logSum (Double log1, Double log2) {
		if(log2==null) return log1;
		if(log1==null) return log2;
		if(log1-log2>MAXLOGDIFF) return log1;
		if(log2-log1>MAXLOGDIFF) return log2;
		return log1 + Math.log10(1+Math.pow(10.0, log2-log1));
	}
	
	/**
	 * Null aware product of two probabilities
	 * @param log1 Log of the first probability
	 * @param log2 Log of the second probability
	 * @return Double sum of the two logarithms. Null if either parameter is null (0*x=0)
	 */
	public static Double logProduct (Double log1, Double log2) {
		if(log1==null || log2==null) return null;
		return log1+log2;
	}
	/**
	 * Null aware 10 power 
	 * @param exponent
	 * @return double power(10,exponent). Zero if exponent is null
	 */
	public static double power10(Double exponent) {
		if(exponent==null) return 0;
		return Math.pow(10.0, exponent);
	}
	/**
	 * Takes the 10-base logarithm of the given value
	 * @param prob Value to take the logarithm
	 * @return Double log10(value). Null if value is less or equal than zero
	 */
	public static Double log10 (double value) {
		if(value > 0) return Math.log10(value);
		else return null;
	}

	public static void normalizeLogs(Double[] logProbs) {
		Double total = null;
		int n = logProbs.length;
		if(n==0) throw new IllegalArgumentException("Array of logarithms must have at least one entry");
		for(int j=0;j<n;j++)  total = LogMath.logSum(total, logProbs[j]);
		for(int j=0;j<n;j++)  logProbs[j] = LogMath.logProduct(logProbs[j],-total);
		
	}
}
