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

import java.util.Random;

/**
 * Implements the Fisher exact test
 * @author Jorge Duitama
 *
 */
public class FisherExactTest {
	private static double [] logFactorials;
	private static boolean quick = true;
	
	public static void main(String[] args) throws Exception {
		Random r = new Random();
		int limit=100000;
		for(int i=0;i<100;i++) {
			int a = r.nextInt(limit);
			int b = r.nextInt(limit);
			int c = r.nextInt(limit);
			int d = r.nextInt(limit);
			quick = false;
			long time = System.currentTimeMillis();
			double p1 = calculatePValue(a, b, c, d);
			time = System.currentTimeMillis()-time;
			System.out.println("Calculated p-value for "+a+" "+b+" "+c+" "+d+" in "+time+" milliseconds. Result: "+p1);
			time = System.currentTimeMillis();
			quick = true;
			double p2 = calculatePValue(a, b, c, d);
			time = System.currentTimeMillis()-time;
			double diff = p1 - p2;
			System.out.println("Calculated p-value on quick mode for "+a+" "+b+" "+c+" "+d+" in "+time+" milliseconds. Result: "+p2+" difference: "+diff);
			if(diff>p1/100) throw new Exception("Quick method not too good"); 
		}
	}
	
	/**
	 * Calculates the p-value of the 2x2 table defined by the given values of a, b, c and d
	 * PRE: a, b, c, and d are non negative
	 * @param a First value of the table.
	 * @param b Second value of the table.
	 * @param c Third value of the table.
	 * @param d Fourth value of the table.
	 * @return double p-value of the table according to the Fisher exact test
	 */
	public static double calculatePValue(int a, int b, int c, int d) {
		//Put smaller value in the top left or bottom right
		if(a>b) {
			a = a + b;
			b = a - b;
			a = a - b;
			c = c + d;
			d = c - d;
			c = c - d;
		}
		if(a>c) {
			a = a + c;
			c = a - c;
			a = a - c;
			b = b + d;
			d = b - d;
			b = b - d;
		}
		int e = Math.min(a,d);
		double answer = 0;
		while (a>=0 && d>=0) {
			double p = calculateExactValue(a, b, c, d);
			if(quick && e>=10 && answer>100*e*p) {
				//Further calculations will not increment the two most significant digits
				break;
			}
			answer+=p;
			a--;
			b++;
			c++;
			d--;
			e++;
		}
		return answer;
	}
	private static void initLogFactorials(int n) {
		if(n<10000) n=10000;
		logFactorials = new double [n+1];
		logFactorials[0] = logFactorials[1] = 0;
		for(int i=2;i<=n;i++) {
			logFactorials[i] = logFactorials[i-1]+LogMath.log10(i);
			//System.out.println("Log factorials ["+i+"] : "+logFactorials[i]);
		}
	}
	
	/**
	 * Calculates the exact probability of the 2x2 table defined by the given values of a, b, c and d
	 * PRE: a, b, c, and d are non negative
	 * @param a First value of the table.
	 * @param b Second value of the table.
	 * @param c Third value of the table.
	 * @param d Fourth value of the table.
	 * @return double p-value of the table according to the Fisher exact test
	 */
	public static double calculateExactValue(int a, int b, int c, int d) {
		int n = a + b + c + d;
		if(logFactorials==null || logFactorials.length<=n) {
			initLogFactorials(n);
		}
		double answer = logFactorials[a+b];
		answer+=logFactorials[c+d];
		answer+=logFactorials[a+c];
		answer+=logFactorials[b+d];
		answer-=logFactorials[a];
		answer-=logFactorials[b];
		answer-=logFactorials[c];
		answer-=logFactorials[d];
		answer-=logFactorials[n];
		return LogMath.power10(answer);
	}
}
