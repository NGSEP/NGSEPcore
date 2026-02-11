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
	private static int CACHE_SIZE = 1000000;
	private static boolean quick = true;
	
	static {
		initLogFactorials(CACHE_SIZE);
	}
	
	public static void main(String[] args) throws Exception {
		Random r = new Random();
		//Test for 2x2 tables
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
		
		//Test for genotypes
		for(int i=0;i<100;i++) {
			int n11 = r.nextInt(limit);
			int n12 = r.nextInt(limit);
			int n22 = r.nextInt(limit);
			quick = false;
			long time = System.currentTimeMillis();
			double p1 = calculateExactTestGenotypeCounts(n11, n12, n22);
			time = System.currentTimeMillis()-time;
			System.out.println("Calculated fisher test for genotype counts "+n11+" "+n12+" "+n22+" in "+time+" milliseconds. Result: "+p1);
			time = System.currentTimeMillis();
			quick = true;
			double p2 = calculateExactTestGenotypeCounts(n11, n12, n22);
			time = System.currentTimeMillis()-time;
			double diff = p1 - p2;
			System.out.println("Calculated fisher test for genotype counts on quick mode for genotype counts"+n11+" "+n12+" "+n22+" in "+time+" milliseconds. Result: "+p2+" difference: "+diff);
			if(diff>p1/100) throw new Exception("Quick method not too good"); 
		}
		//calculateExactTestGenotypeCounts(14262, 59894, 96565);
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
		if(n<CACHE_SIZE) n=CACHE_SIZE;
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
	
	public static double calculateLogCondProbHet (int n12, int N, int n1 ) {
		int n2 = 2*N-n1;
		int n11 = (n1-n12)/2;
		int n22 = N-n11-n12;
		if(logFactorials==null || logFactorials.length<=2*N) {
			initLogFactorials(2*N);
		}
		double answer = n12*LogMath.LOG2;
		answer+=logFactorials[N];
		answer+=logFactorials[n1];
		answer+=logFactorials[n2];
		answer-=logFactorials[2*N];
		answer-=logFactorials[n11];
		answer-=logFactorials[n12];
		answer-=logFactorials[n22];
		return answer;
	}
	public static double calculateExactTestGenotypeCounts (int n11, int n12, int n22 ) {
		if(n12+n22 == 0) return 1;
		if(n11+n12 == 0) return 1;
		int n1 = 2*n11+n12;
		int N = n11+n12+n22;
		int N2 = 2*N;
		int n2 = N2-n1;
		int nMin = Math.min(n1, n2);
		double p1 = ((double)n1)/N2;
		double expectedCountHet = 2.0*p1*(1-p1)*N;
		double next = calculateLogCondProbHet(n12, N, n1);
		//System.out.println("Value for "+n12+" het "+next+" expected  count: "+expectedCountHet);
		double answer = next;
		if(n12 < expectedCountHet) {
			for(int x=n12;x>1;x-=2) {
				next = next + LogMath.log10(x) + LogMath.log10(x-1);
				next = next - LogMath.log10((n1-x)/2+1) - LogMath.log10((n2-x)/2+1);
				next -= LogMath.LOG4;
				if(quick && answer-next>5) break;
				
				answer=LogMath.logSum(answer, next);
				//if(x>59000)System.out.println("Next Value for "+(x-2)+" het "+next+" answer "+answer+" "+LogMath.power10(answer));
			}
		} else {
			for(int x=n12+2;x<=nMin;x+=2) {
				next += LogMath.LOG4;
				next = next + LogMath.log10((n1-x)/2+1) + LogMath.log10((n2-x)/2+1);
				next = next - LogMath.log10(x) - LogMath.log10(x-1);
				if(quick && answer-next>5) break;
				answer=LogMath.logSum(answer, next);
			}
		}
		return LogMath.power10(answer);
	}
	
}
