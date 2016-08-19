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
package ngsep.discovery;

import java.util.Arrays;
import java.util.List;

import JSci.maths.ExtraMath;
import JSci.maths.SpecialMath;
import ngsep.math.PhredScoreHelper;
import ngsep.sequences.DNASequence;


/**
 * Class that stores allele counts and genotype probabilities and offer methods to calculate 
 * posterior probabilities from conditionals 
 * @author Jorge Duitama
 */
public class CountsHelper {
	private int totalCount=0;
	private int lowBaseQualityCount = 0;
	private int [] counts;
	private double [][] logConditionalProbs;
	private double pAllele1=0.5;
	private short maxBaseQS=255;
	private List<String> alleles;
	/**
	 * Creates a default counts helper
	 */
	public CountsHelper () {
		setAlleles(DNASequence.BASES_ARRAY);
	}
	public CountsHelper (String [] alleles) {
		setAlleles(alleles);
	}
	/**
	 * Starts all counts to zero
	 */
	public void startCounts () {
		totalCount=0;
		lowBaseQualityCount=0;
		for(int i=0;i<logConditionalProbs.length;i++) {
			counts[i] = 0;
			for(int j=0;j<logConditionalProbs[0].length;j++) {
				logConditionalProbs[i][j] = 0;
			}
		}
	}
	/**
	 * Updates counts and conditional probabilities for the given allele call
	 * @param allele New allele call to count
	 * @param charQualityScore Quality score of the allele call in Phred+33 scale
	 * @param mappingQuality Quality score of the mapping in Phred scale
	 */
	public void updateCounts (String allele, char charQualityScore, short mappingQuality) {
		totalCount++;
		short qualScore = (short)(charQualityScore-33);
		if(qualScore<=2) {
			lowBaseQualityCount++;
			return;
		} else if (qualScore>maxBaseQS) {
			qualScore = maxBaseQS;
		}
		double readProb = 1;
		double errorProb = PhredScoreHelper.calculateProbability(qualScore);
		double successProb = (1 - errorProb);
		errorProb = errorProb/(alleles.size()-1);
		if(mappingQuality<30) {
			readProb = (1 - Math.pow(10.0, (-mappingQuality/10.0)));
		}
		int index = alleles.indexOf(allele);
		if(index>=0 && successProb > errorProb) {
			counts[index]++;
			double term = pAllele1*(1-alleles.size()*errorProb);
			for(int i=0;i<logConditionalProbs.length;i++) {
				double errorCont = readProb*Math.log10(errorProb);
				if(i==index) {
					logConditionalProbs[i][i] += readProb*Math.log10(successProb); 
				} else {
					logConditionalProbs[i][i] += errorCont;
				}
				for(int j=0;j<logConditionalProbs[0].length;j++) {
					if(i!=j) {
						if(i==index) {
							logConditionalProbs[i][j] += readProb*Math.log10(errorProb+term); 
						} else if (j==index) {
							logConditionalProbs[i][j] += readProb*Math.log10(successProb-term);
						} else {
							logConditionalProbs[i][j] += errorCont;
						}
					}
				}
			}
		}
	}
	public void addAllelicImbalanceFactor(double alpha, double beta) {
		double totalCount = 0;
		for(int i=0;i<counts.length;i++) {
			totalCount+=counts[i];
		}
		if(totalCount==0) return;
		for(int i=0;i<logConditionalProbs.length;i++) {
			for(int j=0;j<logConditionalProbs[0].length;j++) {
				if(i!=j) {
					double pB = 0.000001;
					if(counts[i] > 0 || counts[j] > 0) {
						double beta1 = SpecialMath.beta(alpha+counts[i], beta+counts[j]);
						//System.out.println("Beta1: "+beta1);
						double beta2 = SpecialMath.beta(alpha, beta);
						//System.out.println("Beta2: "+beta2);
						double combinatorial = ExtraMath.binomial((double)(counts[i]+counts[j]), (double)counts[i]);
						pB = combinatorial*beta1/beta2;
					}
					if(pB <0.000001) pB = 0.000001;
					double pLog = Math.log10(pB);
					logConditionalProbs[i][j] += pLog;
					logConditionalProbs[j][i] += pLog;
				}
			}
		}
	}
	/**
	 * 
	 * @return double[][] Conditional probability of each possible genotype  
	 */
	public double[][] getLogConditionalProbs() {
		return logConditionalProbs;
	}

	/**
	 * Calculates the posterior probability of each genotype
	 * @param hetRate Prior heterozygosity rate 
	 * @return double [][] Squared matrix of probabilities.
	 */
	public double [][] getPosteriorProbabilities (double hetRate) {
		int nAlleles = alleles.size();
		//Calculate prior probabilities. Takes into accont alleles order while defining events
		int heteroGenotypes = nAlleles*(nAlleles-1);
		double logPriorHetero = Math.log10(hetRate/heteroGenotypes);
		double logPriorHomo = Math.log10((1-hetRate)/nAlleles);
		
		double [] eventsArray = new double [nAlleles*nAlleles];
		double [][] posteriorProb = new double [nAlleles][nAlleles];
		int indCond=0;
		for(int i=0;i<nAlleles;i++) {
			eventsArray[indCond] = logConditionalProbs[i][i]+logPriorHomo;
			indCond++;
			for(int j=0;j<nAlleles;j++) {
				if(i!=j) {
					eventsArray[indCond] = logConditionalProbs[i][j]+logPriorHetero;
					indCond++;
				}
			}
		}
		calculatePosteriorProbabilities(eventsArray);
		indCond=0;
		for(int i=0;i<nAlleles;i++) {
			posteriorProb[i][i] = eventsArray[indCond];
			indCond++;
			for(int j=0;j<nAlleles;j++) {
				if(i!=j) {
					posteriorProb[i][j] = eventsArray[indCond];
					indCond++;
				}
			}
		}
		return posteriorProb;
	}
	/**
	 * @param eventsArray conditional times prior with logaritmich scale.
	 * This array at the end stores the posterior probabilities
	 */
	private void calculatePosteriorProbabilities (double [] eventsArray) {
		//Calculate max prob
		double logMax = 1;
		for(int i=0;i<eventsArray.length;i++) {
			if(logMax > 0 || logMax < eventsArray[i]) {
				logMax = eventsArray[i];
			}
		}
		//Calculate total probability and contribution of each possible event
		double totalProb = 0;
		for(int i=0;i<eventsArray.length;i++) { 
			eventsArray[i] -= logMax;
			if(eventsArray[i] < -20) {
				eventsArray[i] = 0.0;
			} else {
				eventsArray[i] = Math.pow(10.0, eventsArray[i]);
			}
			totalProb+=eventsArray[i];
		}
		//Normalize
		for(int i=0;i<eventsArray.length;i++) { 
			eventsArray[i] = eventsArray[i]/totalProb;
		}
	}
	/**
	 * Prints the given probabilities. 
	 * @param probs Probabilities matrix to print
	 * @param logs If true, it assumes that the given matrix have log probabilities and then
	 * calculates exp(10,probs[i][j]) before printing
	 */
	public void printProbs (double [][] probs,boolean logs) {
		for(int i=0;i<probs.length;i++) { 
			for(int j=0;j<probs[0].length;j++) {
				if(logs) {
					System.out.print(" "+Math.pow(10.0, probs[i][j]));
				} else {
					System.out.print(" "+probs[i][j]);
				}
				
			}
			System.out.println();
		}
	}
	/**
	 * @return double [] Counts for each possible allele
	 */
	public int[] getCounts() {
		return counts;
	}
	/**
	 * Gets the count for the given allele
	 * @param allele Allele to look for
	 * @return int Count for the given allele
	 */
	public int getCount(String allele) {
		int index = alleles.indexOf(allele);
		if(index>=0) {
			return counts[index];
		}
		return 0;
	}
	/**
	 * @return double Total count of alleles
	 */
	public int getTotalCount() {
		return totalCount;
	}
	public String [] getAlleles() {
		return alleles.toArray(new String [0]);
	}
	public void setAlleles(String [] alleles) {
		this.alleles = Arrays.asList(alleles);
		int nAlleles = alleles.length;
		counts = new int [nAlleles];
		logConditionalProbs = new double [nAlleles][nAlleles];
		startCounts();
	}
	
	public short getMaxBaseQS() {
		return maxBaseQS;
	}
	public void setMaxBaseQS(short maxBaseQS) {
		this.maxBaseQS = maxBaseQS;
	}
	
	public double getpAllele1() {
		return pAllele1;
	}
	public void setpAllele1(double pAllele1) {
		this.pAllele1 = pAllele1;
	}
	public int getLowBaseQualityCount() {
		return lowBaseQualityCount;
	}
	
}
