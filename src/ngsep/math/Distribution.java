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

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ngsep.main.io.ParseUtils;

/**
 * This class is used to accumulate values and infer parameters of a
 * distribution from which the given values could be drawn. Useful to build
 * histograms
 * @author Jorge Duitama
 */
public class Distribution {
	private double sum=0;
	private double sumSquare=0;
	private double count=0;
	private int maxIdx = -1;
	private double distribution [];
	private double cumulative[];
	private boolean cumulativeUpdated = false;
	private double minValueData = Integer.MAX_VALUE;
	private double maxValueData = Integer.MIN_VALUE;
	private List<Double> outliersLess = new ArrayList<Double>();
	private List<Double> outliersMore = new ArrayList<Double>();
	private double minValueDistribution;
	private double maxValueDistribution;
	private double binLength;
	
	public Distribution(double minValueDistribution, double maxValueDistribution, double binLength) {
		super();
		this.minValueDistribution = minValueDistribution;
		this.maxValueDistribution = maxValueDistribution;
		this.binLength = binLength;
		int nBins = (int)((maxValueDistribution-minValueDistribution)/binLength)+1;
		if(nBins<=0) throw new IllegalArgumentException("The given distribution parameters: min "+minValueDistribution+" max "+ maxValueDistribution+" bin length" + binLength+" produce an empty distribution");
		distribution = new double[nBins];
		Arrays.fill(distribution, 0.0);
	}
	public void processDatapoint(double value) {
		processDatapoint(1, value);
	}
	public void processDatapoint(double weigth, double value) {
		double valueW = weigth*value;
		sum+=valueW;
		sumSquare+=(valueW*valueW);
		count+=weigth;
		
		if(value < minValueData ) minValueData = value;
		if(value > maxValueData ) maxValueData = value;
		if(value>=minValueDistribution && value <= maxValueDistribution) {
			int bin = (int)((value-minValueDistribution)/binLength);
			distribution[bin]+=weigth;
			if(maxIdx==-1 || distribution[maxIdx] < distribution[bin] ) maxIdx = bin;
		} else if(value < minValueDistribution) {
			outliersLess.add(value);
		} else if(value > maxValueDistribution) {
			outliersMore.add(value);
		}
		cumulativeUpdated = false;
	}
	public double [] getCumulative () {
		if (cumulativeUpdated) return cumulative;
		cumulative = new double[distribution.length];
		for(int i=0;i<distribution.length;i++) {
			if (i==0) cumulative[i] = outliersLess.size()+distribution[0];
			else cumulative[i] = cumulative[i-1]+distribution[i];
		}
		cumulativeUpdated = true;
		return cumulative;
	}
	public double getSum() {
		return sum;
	}
	public double getSumSquare() {
		return sumSquare;
	}
	public double getCount() {
		return count;
	}
	public double[] getDistribution() {
		return distribution;
	}
	public double getMinValueData() {
		return minValueData;
	}
	public double getMaxValueData() {
		return maxValueData;
	}
	public List<Double> getOutliers() {
		List<Double> outliers = new ArrayList<>();
		outliers.addAll(outliersLess);
		outliers.addAll(outliersMore);
		return outliers;
	}
	public double getMinValueDistribution() {
		return minValueDistribution;
	}
	public double getMaxValueDistribution() {
		return maxValueDistribution;
	}
	public double getBinLength() {
		return binLength;
	}
	public double getAverage() {
		if(count ==0) return 0;
		return sum/count;
	}
	public double getVariance() {
		return (sumSquare-sum*sum/count)/(count-1);
	}
	public double getMaximumFrequency() {
		if(maxIdx==-1) return distribution[0];
		return distribution[maxIdx];
	}
	public double getMaximumBinStart () {
		if(maxIdx==-1) return minValueDistribution;
		return minValueDistribution+maxIdx*binLength;
	}
	public int getBinIndex (int value) {
		return (value-(int)minValueDistribution)/(int)binLength;
	}
	public int getBinIndex (double value) {
		return (int)((value-minValueDistribution)/binLength);
	}
	/**
	 * Calculates the distribution count for the bin where the given value is located
	 * @param value to calculate. If it is an outlier returns the size of the corresponding list of outliers
	 * @return Count of elements falling in the same distribution bin as the given value
	 */
	public double getDistributionCount (double value) {
		int idx = getBinIndex(value);
		if(idx<0) return outliersLess.size();
		if(idx>=distribution.length) return outliersMore.size();
		return distribution[idx];
	}
	/**
	 * Calculates the cumulative count for the given value
	 * @param value to calculate. If it is an outlier returns the size of the left outliers in one case and the complete count in the other end
	 * @return double Count of elements with bin locations less or equal than the bin where the given value is located
	 */
	public double getCumulativeCount (double value) {
		double [] cumulative = getCumulative();
		int idx = getBinIndex(value);
		if(idx<0) return outliersLess.size();
		if(idx>=distribution.length) return count;
		return cumulative[idx];
	}
	/**
	 * Calculates the percentage of datapoints with values equal or more extreme than the given value
	 * @param value to calculate. It must be between the limits of the distribution
	 * @return double Empirical p-value of the given value based on the distribution
	 */
	public double getEmpiricalPvalue (double value) {
		if(count==0) return 1;
		double distCount = getDistributionCount(value);
		double cumulativeCount = getCumulativeCount(value);
		double testCount = Math.max(1, Math.min(cumulativeCount, count-cumulativeCount+distCount));
		return testCount/count;
	}
	/**
	 * Calculates the distribution count for the bin where the given value is located
	 * @param value to calculate. If it is an outlier returns the size of the corresponding list of outliers
	 * @return Count of elements falling in the same distribution bin as the given value
	 */
	public double getDistributionCount (int value) {
		int idx = getBinIndex(value);
		if(idx<0) return outliersLess.size();
		if(idx>=distribution.length) return outliersMore.size();
		return distribution[idx];
	}
	/**
	 * Calculates the cumulative count for the given value
	 * @param value to calculate. If it is an outlier returns the size of the left outliers in one case and the complete count in the other end
	 * @return double Count of elements with bin locations less or equal than the bin where the given value is located
	 */
	public double getCumulativeCount (int value) {
		double [] cumulative = getCumulative();
		int idx = getBinIndex(value);
		if(idx<0) return outliersLess.size();
		if(idx>=distribution.length) return count;
		return cumulative[idx];
	}
	
	/**
	 * Calculates the percentage of datapoints with values equal or more extreme than the given value
	 * @param value to calculate. It must be between the limits of the distribution
	 * @return double Empirical p-value of the given value based on the distribution
	 */
	public double getEmpiricalPvalue (int value) {
		if(count==0) return 1;
		double distCount = getDistributionCount(value);
		double cumulativeCount = getCumulativeCount(value);
		double testCount = Math.max(1, Math.min(cumulativeCount, count-cumulativeCount+distCount));
		return testCount/count;
	}
	
	/**
	 * Return the most popular value within the subset of the distribution defined by the given limits
	 * @param leftValue
	 * @param rightValue
	 * @return
	 */
	public double getLocalMode(double leftValue, double rightValue) {
		int i = getBinIndex(leftValue);
		if(i<0) i=0;
		int endBin = getBinIndex(rightValue)+1;
		if(endBin>distribution.length) endBin = distribution.length;
		int maxBinIdx = i;
		double maxValue = distribution[i];
		while(i<endBin) {
			double value = distribution[i];
			if(value > maxValue) {
				maxValue = value;
				maxBinIdx = i;
			}
			i++;
		}
		return minValueDistribution+(maxBinIdx*binLength);
	}
	/**
	 * Return the less popular value within the subset of the distribution defined by the given limits
	 * @param leftValue
	 * @param rightValue
	 * @return double
	 */
	public double getLocalMinimum(double leftValue, double rightValue) {
		int i = getBinIndex(leftValue);
		if(i<0) i=0;
		int endBin = getBinIndex(rightValue)+1;
		if(endBin>distribution.length) endBin = distribution.length;
		int minBinIdx = i;
		double minValue = distribution[i];
		while(i<endBin) {
			double value = distribution[i];
			if(value < minValue) {
				minValue = value;
				minBinIdx = i;
			}
			i++;
		}
		return minValueDistribution+(minBinIdx*binLength);
	}
	public double getEstimatedStandardDeviationPeak(double peakValue) {
		int peakBin = (int)((peakValue-minValueDistribution)/binLength);
		if(peakBin<0) peakBin = 0;
		if(peakBin >= distribution.length) peakBin = distribution.length-1;
		double peakFrequency = distribution[peakBin];
		//Forward estimation
		double estimationF=maxValueDistribution-peakValue;
		for(int i=peakBin+1;i<distribution.length;i++) {
			if(0.13*peakFrequency>distribution[i]) {
				estimationF = (minValueDistribution+i*binLength-peakValue)/2;
				break;
			}
		}
		//Backward estimation
		double estimationB=peakValue-minValueDistribution;
		for(int i=peakBin-1;i>=0;i--) {
			if(0.13*peakFrequency>distribution[i]) {
				estimationB = (peakValue - (minValueDistribution+i*binLength))/2;
				break;
			}
		}
		return (estimationF+estimationB)/2;
	}
	public void printDistribution(PrintStream out) {
		printDistribution(out,false,maxValueDistribution);
	}
	public void printDistributionInt(PrintStream out) {
		printDistribution(out,true,maxValueDistribution);
	}
	
	public void printDistribution(PrintStream out, boolean integerCounts, double maxValue) {
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
		if(outliersLess.size()>0) {
			out.print("Less\t");
			if(integerCounts) out.println(outliersLess.size());
			else out.println(fmt.format(outliersLess.size()));
		}
		int maxIdx = getBinIndex(maxValue);
		for(int i=0;i<distribution.length && i<=maxIdx;i++) {
			if(integerCounts) {
				int binMinimum = (int)(minValueDistribution+i*binLength);
				out.println(""+binMinimum+"\t"+(int)Math.round(distribution[i]));
			} else {
				double binMinimum = minValueDistribution+i*binLength;
				out.println(""+fmt.format(binMinimum)+"\t"+fmt.format(distribution[i]));
			}
			
		}
		double moreCount = outliersMore.size();
		for(int i=maxIdx+1;i<distribution.length;i++) {
			moreCount+=distribution[i];
		}
		if(moreCount>0) {
			out.print("More\t");
			if(integerCounts) out.println((int)Math.round(moreCount));
			else out.println(fmt.format(moreCount));
		}
		if(count>Integer.MAX_VALUE) out.println("Count\t"+(long)Math.round(count));
		else out.println("Count\t"+(int)Math.round(count));
		
		if(integerCounts) {
			if(sum>Integer.MAX_VALUE) out.println("Sum\t"+(long)Math.round(sum));
			else out.println("Sum\t"+(int)Math.round(sum));
		} else {
			out.println("Sum\t"+fmt.format(sum));
		}
		
		if(count>0) out.println("Average\t"+fmt.format(getAverage()));
		if(count>1) {
			double variance = getVariance();
			out.println("Variance\t"+fmt.format(variance));
			out.println("STDev\t"+fmt.format(Math.sqrt(variance)));
		}
	}
	public void reset() {
		sum=0;
		sumSquare=0;
		count=0;
		maxIdx = -1;
		Arrays.fill(distribution, 0);
		minValueData = Integer.MAX_VALUE;
		maxValueData = Integer.MIN_VALUE;
		outliersLess.clear();
		outliersMore.clear();	
	}

}
