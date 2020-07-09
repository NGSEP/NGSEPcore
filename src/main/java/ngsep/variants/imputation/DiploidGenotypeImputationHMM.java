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
import java.util.List;
import java.util.Map;

import ngsep.hmm.RecombinationHMM;
import ngsep.math.LogMath;
import ngsep.math.NumberArrays;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledSNV;

public class DiploidGenotypeImputationHMM extends RecombinationHMM {
	
	private int startsBaumWelch = HaplotypeClustersHMM.DEF_STARTS_BAUM_WELCH;
	
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
				if(state1.getId()!=null && state2.getId()!=null) pairState.setId(state1.getId()+"/"+state2.getId());
				pairState.setLogStart(logStart);
				pairStates.add(pairState);
			}
		}
		return pairStates;
	}

	public void setAvgCMPerKbp(double avgCMPerKbp) {
		super.setAvgCMPerKbp(avgCMPerKbp);
		haploidBaseHMM.setAvgCMPerKbp(avgCMPerKbp);
	}

	public void setSkipTransitionsTraining(boolean skipTransitionsTraining) {
		super.setSkipTransitionsTraining(skipTransitionsTraining);
		haploidBaseHMM.setSkipTransitionsTraining(skipTransitionsTraining);
	}
	
	public void setTrainingData(List<List<? extends Object>> trainingData) {
		super.setTrainingData(trainingData);
		haploidBaseHMM.setTrainingData(trainingData);
	}
	
	public int getStartsBaumWelch() {
		return startsBaumWelch;
	}

	public void setStartsBaumWelch(int startsBaumWelch) {
		this.startsBaumWelch = startsBaumWelch;
	}

	public void imputeGenotypes (Map<String, List<CalledSNV>> genotypes, int [][][] outClusters) {
		List<String> sampleIds = new ArrayList<String>();
		sampleIds.addAll(genotypes.keySet());
		int n = sampleIds.size();
		int m = genotypes.values().iterator().next().size();
		int k = getNumStates();
		
		if(m!=getSteps()) throw new IllegalArgumentException("Number of variants: "+m+" in the set of genotypes does not coincide with steps of the HMM: "+getSteps());
		double [][][] sumGenotypeProbs = new double [n][m][3];
		double [][][] nextGenotypeProbs = new double [n][m][3];
		
		double [][] nextPosteriorsSample = new double [m][k];
		int [] nextViterbiPathSample = new int [m];
		for(int i=0;i<n;i++) {
			NumberArrays.initializeDoubleMatrix(sumGenotypeProbs[i]);
		}
		for(int h=0;h<startsBaumWelch;h++) {
			getLog().info("Training and sampling iteration: "+h);
			train();
			getLog().info("Model trained");
			for(int i=0;i<n;i++) {
				String sampleId = sampleIds.get(i);
				List<CalledSNV> genotypesSample = genotypes.get(sampleId);
				calculateGenotypePosteriors(genotypesSample, nextGenotypeProbs[i]);
				NumberArrays.accumulate(sumGenotypeProbs[i],nextGenotypeProbs[i]);
				
				//State posteriors for assignments
				calculatePosteriors(genotypesSample, nextPosteriorsSample);
				
				//Best viterbi path
				getViterbiPath(genotypesSample, nextViterbiPathSample);
				
				//Conciliate viterbi with posterior
				assignClusters (sampleId, genotypesSample, nextPosteriorsSample,nextViterbiPathSample,outClusters[h][i]);
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


	public void train() {
		haploidBaseHMM.train();
		getLog().info("Trained internal haploid model ");
		int n = getSteps();
		int kD = getNumStates();
		Double [][] logTransitionsStep = new Double [kD][kD];
		for(int step=0;step<n-1;step++) {
			for(int i = 0;i<kD; i++) {
				HaplotypePairHMMState statePair1 = (HaplotypePairHMMState)getState(i);
				for(int j = 0;j<kD; j++) {
					HaplotypePairHMMState statePair2 = (HaplotypePairHMMState)getState(j);
					Double t1 = haploidBaseHMM.getTransition(statePair1.getIndex1(), statePair2.getIndex1(), step);
					if(t1==null) {
						getLog().info("WARN: Zero transition between: "+statePair1.getIndex1()+" and "+statePair2.getIndex1()+" at step: "+step+" value: "+t1);
					}
					else if(statePair1.getIndex1()!=statePair2.getIndex1() && t1 > -1) {
						getLog().info("WARN: Abnormally high transition between: "+statePair1.getIndex1()+" and "+statePair2.getIndex1()+" at step: "+step+" value: "+t1);
					}
					Double t2 = haploidBaseHMM.getTransition(statePair1.getIndex2(), statePair2.getIndex2(), step);
					if(t2==null) {
						getLog().info("WARN: Zero transition between: "+statePair1.getIndex2()+" and "+statePair2.getIndex2()+" at step: "+step+" value: "+t2);
					}
					else if(statePair1.getIndex2()!=statePair2.getIndex2() && t2 > -1) {
						getLog().info("WARN: Abnormally high transition between: "+statePair1.getIndex2()+" and "+statePair2.getIndex2()+" at step: "+step+" value: "+t2);
					}
					logTransitionsStep[i][j] = LogMath.logProduct(t1, t2);
				}
			}
			//getLog().info("Setting transitions for step: "+step);
			setTransitions(logTransitionsStep, step);
		}
		getLog().info("Trained diploid model ");
	}

	public void calculateGenotypePosteriors(List<CalledSNV> genotypes, double[][] genotypePosteriors) {
		int m = genotypes.size();
		int k = getNumStates();
		Double [][] forwardLogs=calculateForward(genotypes);
		Double [][] backwardLogs=calculateBackward(genotypes);
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

	/**
	 * 
	 * @param nextPosteriors mxk posteriors per site per state
	 * @param viterbiPath m size pat over the states
	 * @param outClusters That will correspond to the states
	 */
	private void assignClusters(String sampleId, List<CalledSNV> genotypesSample, double[][] nextPosteriors, int[] viterbiPath, int[] outClusters) {
		if(viterbiPath[0]==-1) {
			getLog().warning("Sample "+sampleId+" Viterbi state not found.");
			return;
		}
		Byte a0 = 0;
		Byte a1 = 1;
		
		for(int i=0;i<nextPosteriors.length;i++) {
			outClusters[i] = viterbiPath[i];
			CalledSNV call = genotypesSample.get(i);
			//TODO: Phase sites with CN!=2
			if(call.getCopyNumber()!=2) continue;
			HaplotypePairHMMState bestState = (HaplotypePairHMMState) getState(outClusters[i]);
			HaplotypeClusterHMMState hapState1 = bestState.getState1();
			HaplotypeClusterHMMState hapState2 = bestState.getState2();
			//if(i>0 && outClusters[i]!=outClusters[i-1]) {
				//getLog().info("Crossover predicted for sample "+sampleId+" site "+call.getSequenceName()+":"+call.getFirst()+" Previous state: "+outClusters[i-1]+":"+getState(outClusters[i-1]).getId()+" new state "+outClusters[i]+":"+bestState.getId());
			//}
			int idxMaxPos = NumberArrays.getIndexMaximum(nextPosteriors[i]);
			
			if(idxMaxPos==-1) {
				getLog().warning("Sample "+sampleId+" site "+call.getSequenceName()+":"+call.getFirst()+" Viterbi state "+viterbiPath[i]+" best posterior state not found");
			} else if(idxMaxPos==viterbiPath[i]) {
				//Phase het variants
				if(call.isHeterozygous() && !call.isPhased()) {
					Double lp1 = LogMath.logProduct(hapState1.getEmission(a0, i),hapState2.getEmission(a1, i));
					Double lp2 = LogMath.logProduct(hapState1.getEmission(a1, i),hapState2.getEmission(a0, i));
					double p1 = LogMath.power10(lp1);
					double p2 = LogMath.power10(lp2);
					//double pMax = Math.max(p1, p2);
					call.setPhasingCN2(p2>p1);
					/*if( pMax < 0.5) {
						getLog().warning("Sample "+sampleId+" site "+call.getSequenceName()+":"+call.getFirst()+" Phasing with low probability "+pMax+" allele probabilities chosen states: "+hapState1.getEmission(a0, i)+" "+hapState1.getEmission(a1, i)+" "+hapState2.getEmission(a0, i)+" "+hapState2.getEmission(a1, i));
					}*/	
				} else if (call.isHomozygous()) {
					call.setPhasingCN2(!call.isHomozygousReference());
				}
				
			} /*else {
				double pViterbiState = nextPosteriors[i][viterbiPath[i]];
				double pMaxPosterior = nextPosteriors[i][idxMaxPos];
				if(pMaxPosterior-pViterbiState>0.05) {
					//getLog().warning("Sample "+sampleId+" site "+call.getSequenceName()+":"+call.getFirst()+" Viterbi state "+viterbiPath[i]+" different than best posterior state "+idxMaxPos+" viterbiProb: "+pViterbiState+"posterior prob: "+pMaxPosterior);
				}
			}*/
			
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
				//getLog().info("Genotype at "+call.getSequenceName()+":"+call.getFirst()+" inconsistent with prediction. Predicted: "+maxG+" given: "+g);
				inconsistent++;
			}
		}
		getLog().info("Imputed: "+imputed+" inconsistent: "+inconsistent);
	}


	
}
