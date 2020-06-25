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
import java.util.Collections;
import java.util.List;

import JSci.maths.ExtraMath;
import JSci.maths.SpecialMath;
import ngsep.math.FisherExactTest;
import ngsep.math.LogMath;
import ngsep.math.PhredScoreHelper;
import ngsep.sequences.DNASequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;


/**
 * Class that stores allele counts and genotype probabilities and offer methods to calculate 
 * posterior probabilities from conditionals 
 * @author Jorge Duitama
 */
public class CountsHelper {
	
	public static final double DEF_HETEROZYGOSITY_RATE_DIPLOID = 0.001;
	public static final double DEF_HETEROZYGOSITY_RATE_HAPLOID = 0.000001;
	private static final byte DEF_MIN_BASE_QS = 3;
	public static final byte DEF_MAX_BASE_QS = 100;
	private static final double DEF_LOG_ERROR_PROB_INDEL = Math.log10(0.0001);
	
	private int totalCount=0;
	private int lowBaseQualityCount = 0;
	private int [] counts;
	private int [][] countsStrand;
	private double [][] logConditionalProbs;
	private byte maxBaseQS = DEF_MAX_BASE_QS;
	
	private List<String> alleles;
	private static double [][][] logProbCache;
	
	private boolean verbose = false;
	
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
	 * Creates a CountsHelper object with counts of the given alleles according to the given calls
	 * @param alleles Array of possible alleles to take into account
	 * @param calls Allele calls to count
	 * @param maxBaseQS maximum base quality score. Larger quality scores are equalized to this value
	 * @return CountsHelper object with counts and probabilities to call variants
	 */
	public static CountsHelper calculateCounts(String [] alleles, List<PileupAlleleCall> calls, byte maxBaseQS) {
		CountsHelper answer = new CountsHelper(alleles);
		if(maxBaseQS>0) answer.setMaxBaseQS(maxBaseQS);
		for(PileupAlleleCall call: calls) {
			answer.updateCounts(call.getAlleleString(), call.getQualityScores(), call.isNegativeStrand());
		}
		return answer;
	}
	/**
	 * Calculates counts to call SNVs for the given pileup
	 * @param pileup with alignments spanning a given position
	 * @param maxBaseQS maximum base quality score. Larger quality scores are equalized to this value
	 * @param readGroups to return alignments. If null, all alignments of this pileup are processed
	 * @return CountsHelper object with counts and probabilities to call SNVs
	 */
	public static CountsHelper calculateCountsSNV (List<PileupAlleleCall> calls, byte maxBaseQS) {
		CountsHelper answer = new CountsHelper();
		if(maxBaseQS>0) answer.setMaxBaseQS(maxBaseQS);
		for(PileupAlleleCall call:calls ) {
			byte q = (byte)(Math.min(DEF_MAX_BASE_QS, call.getQualityScores().charAt(0)-33));
			answer.updateCounts(call.getAlleleString().substring(0,1), q, call.isNegativeStrand());
		}
		return answer;
	}
	
	/**
	 * @return the verbose
	 */
	public boolean isVerbose() {
		return verbose;
	}
	/**
	 * @param verbose the verbose to set
	 */
	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}
	public void setAlleles(String [] alleles) {
		this.alleles = Arrays.asList(alleles);
		int nAlleles = alleles.length;
		counts = new int [nAlleles];
		countsStrand = new int [nAlleles][2];
		logConditionalProbs = new double [nAlleles][nAlleles];
		updateProbabilitiesCache(nAlleles);
		startCounts();
	}
	private void updateProbabilitiesCache(int n) {
		int m = DEF_MAX_BASE_QS+1;
		if(n<=GenomicVariant.MAX_NUM_ALLELES)n=GenomicVariant.MAX_NUM_ALLELES+1;
		if(logProbCache!=null && logProbCache[0].length>=n) return;
		logProbCache = new double [m][n][3];
		for(byte i=DEF_MIN_BASE_QS;i<m;i++) {
			double errorProb = PhredScoreHelper.calculateProbability(i);
			double successProb = (1 - errorProb);
			logProbCache[i][0][0] = Math.log10(successProb);
			logProbCache[i][0][2] = Math.log10(errorProb);
			for(int j=2;j<n;j++) {
				double epa = errorProb/(j-1);
				logProbCache[i][j][2] = Math.log10(epa);
				double term = 0.5*(1-j*epa);
				logProbCache[i][j][0] = Math.log10(successProb-term);
				logProbCache[i][j][1] = Math.log10(epa+term);
			}
		}
	}
	/**
	 * Starts all counts to zero
	 */
	public void startCounts () {
		totalCount=0;
		lowBaseQualityCount=0;
		for(int i=0;i<logConditionalProbs.length;i++) {
			counts[i] = 0;
			countsStrand [i][0] = countsStrand [i][1] = 0;
			for(int j=0;j<logConditionalProbs[0].length;j++) {
				logConditionalProbs[i][j] = 0;
			}
		}
	}
	/**
	 * Updates counts and conditional probabilities for the given allele call
	 * @param allele New allele call to count
	 * @param qualScore Quality score of the allele call in Phred scale
	 * @param negativeStrand True if the allele comes from a read aligned to the negative strand
	 */
	public void updateCounts (String allele, byte qualScore, boolean negativeStrand) {
		totalCount++;
		if(qualScore<=DEF_MIN_BASE_QS) {
			lowBaseQualityCount++;
			return;
		} else if (qualScore>maxBaseQS) {
			qualScore = maxBaseQS;
		}
		int index = alleles.indexOf(allele);
		if(index>=0) {
			//Update raw count
			counts[index]++;
			//Update strand counts
			if(negativeStrand) countsStrand[index][0]++;
			else countsStrand[index][1]++;
			
			int n = alleles.size();
			//Update probabilities
			for(int i=0;i<logConditionalProbs.length;i++) {
				if(i==index) {
					logConditionalProbs[i][i] += logProbCache[qualScore][0][0]; 
				} else {
					logConditionalProbs[i][i] += logProbCache[qualScore][n][2];
				}
				for(int j=0;j<logConditionalProbs[i].length;j++) {
					if(i!=j) {
						if(i==index) {
							logConditionalProbs[i][j] += logProbCache[qualScore][n][1];
						} else if (j==index) {
							logConditionalProbs[i][j] += logProbCache[qualScore][n][0];
						} else {
							logConditionalProbs[i][j] += logProbCache[qualScore][n][2];
						}
					}
				}
			}
		}
	}
	public void updateCounts(String call, String qualityScores, boolean negativeStrand) {
		totalCount++;
		int index = alleles.indexOf(call);
		if(index>=0) {
			//Update raw count
			counts[index]++;
			//Update strand counts
			if(negativeStrand) countsStrand[index][0]++;
			else countsStrand[index][1]++;
		}
			
		int n = alleles.size();
		double [] conditionals = new double [n];
		for(int i=0;i<logConditionalProbs.length;i++) {
			String alleleI = alleles.get(i);
			if(alleleI.length()==call.length()) {
				conditionals[i] = calculateConditional(alleleI,call,qualityScores); 
			} else {
				conditionals[i] = DEF_LOG_ERROR_PROB_INDEL;
			}
			if(verbose) System.out.println("Allele: "+alleleI+" call: "+call+ " log cond: "+conditionals[i]);
		}
		//Update probabilities
		for(int i=0;i<logConditionalProbs.length;i++) {
			logConditionalProbs[i][i] += conditionals[i];
			for(int j=i+1;j<logConditionalProbs[i].length;j++) {
				double average = LogMath.logSum(conditionals[i], conditionals[j])-0.3;
				logConditionalProbs[i][j]+=average;
				logConditionalProbs[j][i]+=average;
			}
		}
		if(verbose) printProbs(logConditionalProbs, true);
	}
	/**
	 * PRE: allele.length()==call.length()==qualityScores.length()
	 * @param allele
	 * @param call
	 * @param qualityScores phred+33 format
	 * @return double log conditional probability of the given call with the given scores given the allele
	 */
	private double calculateConditional(String allele, String call, String qualityScores) {
		double logCond = 0;
		for(int i=0;i<allele.length();i++) {
			char c = call.charAt(i);
			byte qualScore = (byte)Math.min(maxBaseQS, (qualityScores.charAt(i)-33));
			if(allele.charAt(i)==c) logCond+=logProbCache[qualScore][0][0];
			else logCond+=logProbCache[qualScore][4][2];
		}
		return logCond;
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
	
	public List<String> getAllelesList() {
		return Collections.unmodifiableList(alleles);
	}
	
	public byte getMaxBaseQS() {
		return maxBaseQS;
	}
	public void setMaxBaseQS(byte maxBaseQS) {
		this.maxBaseQS = maxBaseQS;
	}
	
	public int getLowBaseQualityCount() {
		return lowBaseQualityCount;
	}
	
	public double getPValueStrandBiasFisher (int i1, int i2) {
		if(i1<0 || i1>=countsStrand.length) throw new IllegalArgumentException("Invalid first allele index: "+i1);
		if(i2<0 || i2>=countsStrand.length) throw new IllegalArgumentException("Invalid second allele index: "+i2);
		int a = countsStrand[i1][0];
		int b = countsStrand[i2][0];
		int c = countsStrand[i1][1];
		int d = countsStrand[i2][1];
		return FisherExactTest.calculatePValue(a, b, c, d);
	}
	
	public byte getScoreStrandBiasFisher (int i1, int i2) {
		double pvalueSB = getPValueStrandBiasFisher(i1, i2);
		return (byte) Math.min(CalledGenomicVariant.MAX_STRAND_BIAS_SCORE, PhredScoreHelper.calculatePhredScore(pvalueSB));
	}
	
	
}
