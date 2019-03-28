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
import java.util.List;
import java.util.Map;
import JSci.maths.statistics.FDistribution;

/**
 * 
 * @author Andrea Parra
 *
 */
public class FactorialAnova {
	private double sumSquareTotal;
	private double treatmentSumSquares;
	private double errorSumSquares;
	private double r2;
	private double treatmentMeanSquares;
	private double errorMeanSquares;
	private double fStatistic;
	private FDistribution fdist;
	private double pValue;
	private List<Double> y;
	private List<String> x;

	
	private double average;
	private int numberOfLevels;
	private int numberOfSamples;
	private Map<String, Integer> xCountMap; //Describes the number of samples per step
	private Map<String, Double> averageMap; //Describes the average per step
	private Map<String, Double> ySumMap; //Describes the average per step
	private Map<String, Double> varianceMap; //Describes the variance per step
	
	
	
	/**
	 * Pre: |x|=|y|
	 * @param x
	 * @param y
	 */
	
	//-----------------------------------------------------------------
	// 							CONSTRUCTOR
	// ----------------------------------------------------------------
			
	public FactorialAnova (List<String> _x, List<Double> _y) {
		
		y = _y;
		x = _x;
		
		numberOfSamples = _x.size();
		average = calculateAverage();
		xCountMap = generateXCountMap();
		ySumMap = generateYSumMap();
		numberOfLevels = this.xCountMap.keySet().size();
		averageMap = generateAverageMap(this.ySumMap, this.xCountMap);
		varianceMap = generateVarianceMap(this.averageMap, this.xCountMap);
		errorSumSquares = calculateSE(); //change names to real names
		treatmentSumSquares = calculateSR();
		sumSquareTotal = this.treatmentSumSquares + this.errorSumSquares;
		r2 = this.treatmentSumSquares / this.sumSquareTotal;	
		fdist = new FDistribution(this.numberOfLevels-1,this.numberOfSamples-1);
		errorMeanSquares = this.errorSumSquares / (this.numberOfSamples - 1);
		treatmentMeanSquares = this.treatmentSumSquares / (this.numberOfLevels - 1);
		pValue = 1-fdist.cumulative(this.treatmentMeanSquares/this.errorMeanSquares);
			
	}

	public double getSumSquareTotal() {
		return sumSquareTotal;
	}

	public double getTreatmentSumSquares() {
		return treatmentSumSquares;
	}


	public double getErrorSumSquares() {
		return errorSumSquares;
	}

	public double getR2() {
		return r2;
	}

	public double getTreatmentMeanSquares() {
		return treatmentMeanSquares;
	}

	public double getErrorMeanSquares() {
		return errorMeanSquares;
	}


	public double getpValue() {
		return pValue;
	}

	public FDistribution getpFdist() {
		return fdist;
	}


	/**
	 * calculates the overall average of the dataset
	 * @return average (double)
	 */
	private double calculateAverage() {
		double total = 0.0;
		for(double i : this.y) {
			total = total + i;
		}
		double avg = total / this.numberOfSamples;
		return avg;
	}
	
	
	/**
	 * Creates a Hashmap with the keys being the alleles for a given sample: (i.e R, RA, A for biallelic SNPs)
	 * The value is the number of such variations in the sample pool.  
	 * @return a HashMap containing the number of each allelic variation.
	 */
	private Map<String, Integer> generateXCountMap() {
		
		Map<String, Integer> map = new HashMap<String, Integer>();
		for(String str : this.x) {
			map.putIfAbsent(str, 0);
		}
	
		for(String allele : this.x) {
			Integer count = map.get(allele);
			map.put(allele, count + 1);
		}
		return map;
	}
	
	/**
	 * Creates a Hashmap with the keys being the alleles for a given sample: (i.e R, RA, A for biallelic SNPs)
	 * The value is the sum of the fenotypes for such variation. 
	 * @return: a HashMap containing the sum of the fenotypes corresponding each allelic variation 
	 */
	private Map<String, Double> generateYSumMap() {
		Map<String, Double> map = new HashMap<String, Double>();
		
		// Initialize HashMap
		for(String str : this.x) {
			map.putIfAbsent(str, 0.0);
		}
		
		double _y;
		String _x;
		double sum;
		for(int i=0;i<this.numberOfSamples;i++) {
			_x = this.x.get(i);
			_y = this.y.get(i);
			sum = map.get(_x);
			map.put(_x, sum + _y);
		}
		return map;
	}
	
	/**
	 * Creates a Hashmap with the keys being the alleles for a given sample: (i.e R, RA, A for biallelic SNPs)
	 * The value is the average of the fenotypes for such variation. 
	 * @return: a HashMap containing the average of the fenotypes corresponding each allelic variation 
	 */
	private Map<String, Double> generateAverageMap(Map<String, Double> sumMap, Map<String, Integer> countMap) {
		Map<String, Double> map = new HashMap<String, Double>();
		
		// Initialize HashMap
		for(String str : this.x) {
			map.putIfAbsent(str, 0.0);
		}
		
		double sum;
		int count;
		double average;
		for(String key : sumMap.keySet()) {
			sum = sumMap.get(key);
			count = countMap.get(key);
			average = sum / count;
			map.put(key, average);
		}
		return map;
	}
	
	/**
	 * Creates a Hashmap with the keys being the alleles for a given sample: (i.e R, RA, A for biallelic SNPs)
	 * The value is the variance of the fenotypes for such variation. 
	 * TODO: check if the variance to be used is population or sample variance. 
	 * @return: a HashMap containing the variance of the fenotypes corresponding each allelic variation 
	 */
	private Map<String, Double> generateVarianceMap(Map<String, Double> avgMap, Map<String, Integer> countMap){

		Map<String, Double> map = new HashMap<String, Double>();
		double _y;
		String _x;
		double tmp;
		
		// Initialize HashMap
		for(String str : this.x) {
			map.putIfAbsent(str, 0.0);
		}
		
		for(int i=0; i<this.numberOfSamples;i++) {
			tmp = 0.0;
			_x = this.x.get(i);
			_y = this.y.get(i);
			tmp = map.get(_x) + (Math.pow((_y - avgMap.get(_x)), 2) / countMap.get(_x));
			map.put(_x, tmp);
		}
		return map;
	}
	
	/**
	 * Calculates SR. SR = SUM(J[i]*(y_bar[i] - y_bar)^2)
	 * @return SR (double)
	 */
	private double calculateSR() {
		double sum = 0.0;
		for(String str : this.averageMap.keySet()) {
			sum = sum + (this.xCountMap.get(str) * Math.pow((this.averageMap.get(str) - this.average), 2));
		}
		double SR = sum;
		return SR;
	}
	
	/**
	 * Calculates SE, SE = SUM((J[i] - 1) * s[i]^2)
	 * @return SE (double)
	 */
	private double calculateSE() {
		double sum = 0.0;
		for(String str : this.averageMap.keySet()) {
			sum = sum + ((this.xCountMap.get(str) - 1) * Math.pow(this.varianceMap.get(str), 2));
		}
		double SE = sum;
		return SE;
	
	}
}
