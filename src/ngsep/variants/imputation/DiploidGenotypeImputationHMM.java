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
package ngsep.variants.imputation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.hmm.RecombinationHMM;
import ngsep.math.LogMath;
import ngsep.math.NumberArrays;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledSNV;

public class DiploidGenotypeImputationHMM extends RecombinationHMM {
	
	//private int startsBaumWelch = HaplotypeClustersHMM.DEF_STARTS_BAUM_WELCH;
	private int startsBaumWelch = 1;
	
	private Double [][] forwardLogs=new Double[0][0];
	private Double [][] backwardLogs=new Double[0][0];
	private HaplotypeClustersHMM haploidBaseHMM;
	public DiploidGenotypeImputationHMM(HaplotypeClustersHMM baseHMM, List<? extends HaplotypePairHMMState> states, int numMarkers, List<Integer> positions) {
		super(states, numMarkers, positions);
		haploidBaseHMM = baseHMM;
		
	}

	public static DiploidGenotypeImputationHMM createHMM (Map<String, List<CalledSNV>> genotypes, List<String> parentIds, int k, boolean inbreds) {
		List<CalledSNV> snvs = genotypes.values().iterator().next();
		List<Integer> positions = new ArrayList<Integer>();
		for(CalledSNV snv:snvs) positions.add(snv.getFirst());
		int m = snvs.size();
		System.out.println("Number of variants: "+m);
		List<HaplotypeClusterHMMState> statesHaploid = HaplotypeClusterHMMState.createEmptyStates(m, k);
		System.out.println("Created "+statesHaploid.size()+" empty states");
		if(inbreds) HaplotypeClusterHMMState.fillHMMStatesInbred(statesHaploid, genotypes, parentIds);
		else HaplotypeClusterHMMState.fillHMMStatesPhasedDiploid(statesHaploid, genotypes, parentIds);
		System.out.println("Filled states with parents");
		List<HaplotypePairHMMState> pairStates = createHaplotypePairStates(statesHaploid);
		System.out.println("Created states for pairs of samples");
		HaplotypeClustersHMM baseHMM = new HaplotypeClustersHMM(statesHaploid, m, positions);
		System.out.println("Created haploid HMM");
		DiploidGenotypeImputationHMM diploidHMM = new DiploidGenotypeImputationHMM(baseHMM, pairStates, m, positions);
		System.out.println("Created diploid HMM");
		return diploidHMM;
	}

	private static List<HaplotypePairHMMState> createHaplotypePairStates(List<HaplotypeClusterHMMState> haploidStates) {
		List<HaplotypePairHMMState> pairStates = new ArrayList<>();
		int k = haploidStates.size();
		double logStart = LogMath.log10(1.0/(k*k));
		for (int i=0;i<haploidStates.size();i++) {
			HaplotypeClusterHMMState state1 = haploidStates.get(i);
			for (int j=0;j<haploidStates.size();j++) {
				HaplotypeClusterHMMState state2 = haploidStates.get(j);
				HaplotypePairHMMState pairState = new HaplotypePairHMMState(i,state1, j,state2);
				pairState.setLogStart(logStart);
				pairStates.add(pairState);
			}
		}
		return pairStates;
	}

	public double getAvgCMPerKbp() {
		return haploidBaseHMM.getAvgCMPerKbp();
	}

	public void setAvgCMPerKbp(double avgCMPerKbp) {
		super.setAvgCMPerKbp(avgCMPerKbp);
		haploidBaseHMM.setAvgCMPerKbp(avgCMPerKbp);
	}

	public boolean isFixedTransitions() {
		return haploidBaseHMM.isFixedTransitions();
	}

	public void setFixedTransitions(boolean fixedTransitions) {
		super.setFixedTransitions(fixedTransitions);
		haploidBaseHMM.setFixedTransitions(fixedTransitions);
	}
	
	public void imputeGenotypes (Map<String, List<CalledSNV>> genotypes) {
		List<String> sampleIds = new ArrayList<String>();
		sampleIds.addAll(genotypes.keySet());
		int n = sampleIds.size();
		
		int m = genotypes.values().iterator().next().size();
		if(m!=getSteps()) throw new IllegalArgumentException("Number of variants: "+m+" in the set of genotypes does not coincide with steps of the HMM: "+getSteps());
		double [][][] sumGenotypeProbs = new double [n][m][3];
		double [][][] nextGenotypeProbs = new double [n][m][3];
		//TODO: Check memory and use for diploid imputation
		//double [][][] sumStateProbs = new double [n][m][getNumStates()];
		//Double [][][] nextStateLogProbs = new Double [n][m][getNumStates()];
		for(int i=0;i<n;i++) {
			for(int j=0;j<m;j++) {
				Arrays.fill(sumGenotypeProbs[i][j], 0.0);
				//Arrays.fill(sumStateProbs[i][j], 0.0);
			}
		}
		//if(getReferenceHaplotypes() == null) setReferenceHaplotypes(makeHaplotypesWithHomozygous(genotypes));
		for(int h=0;h<startsBaumWelch;h++) {
			getLog().info("Training and sampling iteration: "+h);
			train();
			getLog().info("Model trained");
			for(int i=0;i<n;i++) {
				String sampleId = sampleIds.get(i);
				List<CalledSNV> genotypesSample = genotypes.get(sampleId);
				calculateGenotypePosteriors(genotypesSample, nextGenotypeProbs[i]);
				//calculatePosteriors(genotypesSample, nextStateLogProbs[i]);
				NumberArrays.accumulate(sumGenotypeProbs[i],nextGenotypeProbs[i]);
				//NumberArrays.accumulate(sumStateProbs[i],nextStateProbs[i]);
				getLog().info("Calculated posteriors for sample: "+sampleId);
				
			}
		}
		//Map<String,List<Integer>> assignments = new TreeMap<String, List<Integer>>();
		for(int i=0;i<n;i++) {
			String sampleId = sampleIds.get(i);
			getLog().info("Choosing best genotypes for sample: "+sampleId);
			imputeGenotypes(genotypes.get(sampleId),sumGenotypeProbs[i],startsBaumWelch);
			//assignments.put(sampleId, calculateFinalAssignments(sumStateProbs[i]));
		}
	}
	private void train() {
		//TODO: Decide when to train
		//haploidBaseHMM.train();
		haploidBaseHMM.setRandomTransitions();
		//TODO: Update transitions this HMM
		getLog().info("Trained haploid model ");
		int n = getSteps();
		int kD = getNumStates();
		Double [][] logTransitionsStep = new Double [kD][kD];
		for(int step=0;step<n-1;step++) {
			for(int i = 0;i<kD; i++) {
				HaplotypePairHMMState statePair1 = (HaplotypePairHMMState)getState(i);
				for(int j = 0;j<kD; j++) {
					HaplotypePairHMMState statePair2 = (HaplotypePairHMMState)getState(j);
					Double t1 = haploidBaseHMM.getTransition(statePair1.getIndex1(), statePair2.getIndex1(), step);
					if(statePair1.getIndex1()!=statePair2.getIndex1() && t1 > -1) {
						getLog().info("WARN: Abnormally high transition between: "+statePair1.getState1().getId()+" and "+statePair2.getState1().getId()+" at step: "+step+" value: "+t1);
					}
					Double t2 = haploidBaseHMM.getTransition(statePair1.getIndex2(), statePair2.getIndex2(), step);
					if(statePair1.getIndex2()!=statePair2.getIndex2() && t2 > -1) {
						//getLog().info("WARN: Abnormally high transition between: "+statePair1.getState2().getId()+" and "+statePair2.getState2().getId()+" at step: "+step+" value: "+t2);
					}
					logTransitionsStep[i][j] = LogMath.logProduct(t1, t2);
				}
			}
			setTransitions(logTransitionsStep, step);
		}		
	}

	public void calculateGenotypePosteriors(List<CalledSNV> genotypes, double[][] genotypePosteriors) {
		int m = genotypes.size();
		int k = getNumStates();
		initArrays(k, m);
		calculateForward(genotypes, forwardLogs);
		calculateBackward(genotypes, backwardLogs);
		for(int i=0;i<m;i++) {
			Double log0 = null;
			Double log1 = null;
			Double log2 = null;
			for(int j=0;j<k;j++) {
				Double fTimesB = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i][j]);
				log0 = LogMath.logSum(log0, LogMath.logProduct(fTimesB, getEmission(j, CalledSNV.GENOTYPE_HOMOREF, i)));
				log1 = LogMath.logSum(log1, LogMath.logProduct(fTimesB, getEmission(j, CalledSNV.GENOTYPE_HETERO, i)));
				log2 = LogMath.logSum(log2, LogMath.logProduct(fTimesB, getEmission(j, CalledSNV.GENOTYPE_HOMOALT, i)));
			}
			//Normalize and raise to calculate final probabilities of genotypes
			Double logSum = LogMath.logSum(log0, log1);
			logSum = LogMath.logSum(logSum, log2);
			double prob0 = LogMath.power10(LogMath.logProduct(log0, -logSum));
			double prob1 = LogMath.power10(LogMath.logProduct(log1, -logSum));
			double prob2 = LogMath.power10(LogMath.logProduct(log2, -logSum));
			double sum = prob0 + prob1 + prob2;
			prob0/=sum;
			prob1/=sum;
			prob2/=sum;
			genotypePosteriors[i][0] = prob0;
			genotypePosteriors[i][1] = prob1;
			genotypePosteriors[i][2] = prob2;
		}
	}
	
	private void imputeGenotypes(List<CalledSNV> genotypes, double[][] sums, int steps) {
		
		int imputed = 0;
		int inconsistent = 0;
		for(int i=0;i<sums.length;i++) {
			CalledSNV call = genotypes.get(i);
			byte g = call.getGenotype();
			byte maxG = 0;
			for(int j=0;j<sums[i].length;j++) {
				if(sums[i][maxG]<sums[i][j]) {
					maxG = (byte)j;
				}
			}
			if(call.isUndecided() ) {
				//Make genotype out of maximum allele
				call.setGenotype(maxG);
				double prob = sums[i][maxG]/steps;
				call.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(1-prob));
				imputed++;
			} else if(g!=maxG) {
				getLog().info("Genotype at "+call.getSequenceName()+":"+call.getFirst()+" inconsistent with prediction. Predicted: "+maxG+" given: "+g);
				inconsistent++;
			}
		}
		getLog().info("Imputed: "+imputed+" inconsistent: "+inconsistent);
	}
	
	
	
	private void initArrays(int k, int m) {
		if(forwardLogs.length!=m || forwardLogs[0].length!=k) {
			getLog().info("Creating array for forward probabilities of dimensions "+m+" x "+k);
			forwardLogs = new Double[m][k];
		}
		if(backwardLogs.length!=m || backwardLogs[0].length!=k) {
			getLog().info("Creating array for backward probabilities of dimensions "+m+" x "+k);
			backwardLogs = new Double[m][k];
		}
	}


	
}
