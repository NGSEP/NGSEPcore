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
package ngsep.variants;

import java.util.Arrays;
import java.util.List;

import JSci.maths.statistics.ChiSqrDistribution;

public class DiversityStatistics {
	private static ChiSqrDistribution chiSquare = new ChiSqrDistribution(1);
	private int numSamplesGenotyped;
	private double [] alleleFrequencies;
	private int [] alleleCounts;
	private int numCalledAlleles;
	private double expectedHeterozygosity;
	private double observedHeterozygosity;
	private double fStatistic;
	private double maf=0;
	private byte mafIndex = -1;
	//Index of the wild type allele
	private byte wtIndex = -1;
	//Statistics for biallelic samples
	
	private double chiSquareValue = 0;
	private double chiSquarePValue = 1;
	
	
	public int getNumSamplesGenotyped() {
		return numSamplesGenotyped;
	}
	public void setNumSamplesGenotyped(int numSamplesGenotyped) {
		this.numSamplesGenotyped = numSamplesGenotyped;
	}
	
	public int[] getAlleleCounts() {
		return alleleCounts;
	}
	public void setAlleleCounts(int[] alleleCounts) {
		this.alleleCounts = alleleCounts;
	}
	public double[] getAlleleFrequencies() {
		return alleleFrequencies;
	}
	public void setAlleleFrequencies(double[] alleleFrequencies) {
		this.alleleFrequencies = alleleFrequencies;
	}
	public double getExpectedHeterozygosity() {
		return expectedHeterozygosity;
	}
	public void setExpectedHeterozygosity(double expectedHeterozygosity) {
		this.expectedHeterozygosity = expectedHeterozygosity;
	}
	public double getObservedHeterozygosity() {
		return observedHeterozygosity;
	}
	public void setObservedHeterozygosity(double observedHeterozygosity) {
		this.observedHeterozygosity = observedHeterozygosity;
	}
	
	public double getfStatistic() {
		return fStatistic;
	}
	public void setfStatistic(double fStatistic) {
		this.fStatistic = fStatistic;
	}
	public double getMaf() {
		return maf;
	}
	public void setMaf(double maf) {
		this.maf = maf;
	}
	public byte getMafIndex() {
		return mafIndex;
	}
	public void setMafIndex(byte mafIndex) {
		this.mafIndex = mafIndex;
	}
	public byte getWtIndex() {
		return wtIndex;
	}
	public void setWtIndex(byte wtIndex) {
		this.wtIndex = wtIndex;
	}
	public double getChiSquareValue() {
		return chiSquareValue;
	}
	public void setChiSquareValue(double chiSquareValue) {
		this.chiSquareValue = chiSquareValue;
	}
	public double getChiSquarePValue() {
		return chiSquarePValue;
	}
	public void setChiSquarePValue(double chiSquarePValue) {
		this.chiSquarePValue = chiSquarePValue;
	}
	public int getNumCalledAlleles() {
		return numCalledAlleles;
	}
	public void setNumCalledAlleles(int numCalledAlleles) {
		this.numCalledAlleles = numCalledAlleles;
	}
	public boolean isPolymorphic() {
		return numCalledAlleles>1;
	}
	public static DiversityStatistics calculateDiversityStatistics(List<CalledGenomicVariant> genotypeCalls, boolean assumeAlwaysDiploid) {
		String [] alleles = genotypeCalls.get(0).getAlleles();
		int [] counts = new int[alleles.length];
		Arrays.fill(counts, 0);
		double sum = 0;
		int numGenotyped = 0;
		
		int numHeterozygous = 0;
		int [] numHomozygous = new int[alleles.length];
		for(CalledGenomicVariant calledVar:genotypeCalls) {
			if(calledVar.isUndecided()) continue;
			numGenotyped++;
			if(calledVar.isHeterozygous()) numHeterozygous++;
			byte [] indexesCalledAlleles = calledVar.getIndexesCalledAlleles();
			if(calledVar.isHomozygous()) numHomozygous[indexesCalledAlleles[0]]++;
			if(assumeAlwaysDiploid) {
				if(indexesCalledAlleles.length==1) counts[indexesCalledAlleles[0]]+=2;
				else {
					counts[indexesCalledAlleles[0]]++;
					counts[indexesCalledAlleles[1]]++;
				}
				sum+=2;
			} else {
				short [] allelesCN = calledVar.getAllelesCopyNumber();
				for(int i=0;i<indexesCalledAlleles.length;i++) {
					int j=indexesCalledAlleles[i];
					if(j>=0 && j<counts.length) {
						counts[j]+=allelesCN[j];
						sum+=allelesCN[j];
					}
					
				}
			}
		}
		double afs [] = new double[alleles.length];
		double expectedHet = 1;
		int numCalledAlleles = 0;
		int minAC = 0;
		byte minACIdx = -1;
		int maxAC = 0;
		byte maxACIdx = -1;
		for(int i=0;i<alleles.length;i++) {
			afs[i] = 0;
			if(sum>0) afs[i] = ((double)counts[i])/sum;
			if(counts[i]>0) {
				numCalledAlleles++;
				if(minAC == 0 || minAC>counts[i]) {
					minAC = counts[i];
					minACIdx = (byte) i;
				}
				if(maxACIdx == -1 || maxAC<counts[i]) {
					maxAC = counts[i];
					maxACIdx = (byte) i;
				}
			}
			expectedHet-=(afs[i]*afs[i]);
		}
		
		double observedHet = 0;
		if(numGenotyped > 0) observedHet = (double)numHeterozygous/numGenotyped;
		DiversityStatistics answer = new DiversityStatistics();
		answer.setAlleleCounts(counts);
		answer.setAlleleFrequencies(afs);
		answer.setNumCalledAlleles(numCalledAlleles);
		answer.setNumSamplesGenotyped(numGenotyped);
		answer.setExpectedHeterozygosity(expectedHet);
		answer.setObservedHeterozygosity(observedHet);
		if(numCalledAlleles<2) {
			//Monomorphic
			answer.setMaf(0);
			answer.setMafIndex((byte) -1);
			if(numCalledAlleles==1) answer.setWtIndex(maxACIdx);
		} else {
			answer.setMaf((double)minAC/sum);
			answer.setMafIndex(minACIdx);
			answer.setWtIndex(maxACIdx);
		}
		if(expectedHet>0.0001)answer.setfStatistic(1-observedHet/expectedHet);
		else answer.setfStatistic(0);
		if(alleles.length==2 && numCalledAlleles==2 && numGenotyped>3) {
			int expCountHomo0 = (int)Math.round(afs[0]*afs[0]*numGenotyped);
			if(expCountHomo0==0) expCountHomo0 = 1;
			else if (expCountHomo0>numGenotyped-2) expCountHomo0=numGenotyped-2;
			int expCountHomo1 = (int)Math.round(afs[1]*afs[1]*numGenotyped);
			if(expCountHomo1==0) expCountHomo1 = 1;
			else if (expCountHomo0+expCountHomo1>numGenotyped-1) expCountHomo1 = numGenotyped-expCountHomo0-1;
			int expCountHetero = numGenotyped - expCountHomo0 - expCountHomo1;
			double chiValue = Math.pow(numHomozygous[0]-expCountHomo0,2)/expCountHomo0;
			chiValue += Math.pow(numHomozygous[1]-expCountHomo1,2)/expCountHomo1;
			chiValue += Math.pow(numHeterozygous-expCountHetero, 2)/expCountHetero;
			answer.setChiSquareValue(chiValue);
			answer.setChiSquarePValue(1-chiSquare.cumulative(chiValue));
		}
		return answer;
	}
	
}
