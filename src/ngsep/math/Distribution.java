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
	public double getLocalMode(double leftValue, double rightValue) {
		int i = (int)((leftValue-minValueDistribution)/binLength);
		if(i<0) i=0;
		int endBin = (int)((rightValue-minValueDistribution)/binLength)+1;
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
	public void printDistribution(PrintStream out,double maxValue) {
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
		int valueBin = (int)((maxValue-minValueDistribution)/binLength);
		for(int i=0;i<distribution.length && i<=valueBin;i++) {
			double binMinimum = minValueDistribution+i*binLength;
			out.println(""+fmt.format(binMinimum)+"\t"+fmt.format(distribution[i]));
		}
	}
	public void printDistribution(PrintStream out) {
		printDistribution(out,false);
	}
	public void printDistributionInt(PrintStream out) {
		printDistribution(out,true);
	}
	
	public void printDistribution(PrintStream out, boolean integerCounts) {
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
		if(outliersLess.size()>0) {
			out.print("Less\t");
			if(integerCounts) out.println(outliersLess.size());
			else out.println(fmt.format(outliersLess.size()));
		}
		for(int i=0;i<distribution.length;i++) {
			if(integerCounts) {
				int binMinimum = (int)(minValueDistribution+i*binLength);
				out.println(""+binMinimum+"\t"+(int)distribution[i]);
			} else {
				double binMinimum = minValueDistribution+i*binLength;
				out.println(""+fmt.format(binMinimum)+"\t"+fmt.format(distribution[i]));
			}
			
		}
		if(outliersMore.size()>0) {
			out.print("More\t");
			if(integerCounts) out.println(outliersMore.size());
			else out.println(fmt.format(outliersMore.size()));
		}
		out.println("Count\t"+(int)count);
		if(integerCounts) {
			out.println("Sum\t"+(int)sum);
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
