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

import ngsep.math.NumberArrays;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledSNV;


public class GenotypeImputationHMM extends HaplotypeClustersHMM {
	
	private int startsBaumWelch = HaplotypeClustersHMM.DEF_STARTS_BAUM_WELCH;
	
	public GenotypeImputationHMM(List<? extends HaplotypeClusterHMMState> states, int numMarkers) {
		super(states, numMarkers);
	}
	
	public GenotypeImputationHMM(List<? extends HaplotypeClusterHMMState> states, int numMarkers, List<Integer> positions) {
		super(states, numMarkers,positions);
	}
	
	public static GenotypeImputationHMM createHMM (Map<String, List<CalledSNV>> genotypes, List<String> parentIds, int k, boolean inbreds) {
		List<CalledSNV> snvs = genotypes.values().iterator().next();
		List<Integer> positions = new ArrayList<Integer>();
		for(CalledSNV snv:snvs) positions.add(snv.getFirst());
		int m = snvs.size();
		List<HaplotypeClusterHMMState> states = HaplotypeClusterHMMState.createEmptyStates(m, k);
		if(inbreds) HaplotypeClusterHMMState.fillHMMStatesInbred(states, genotypes, parentIds);
		else HaplotypeClusterHMMState.fillHMMStatesPhasedDiploid(states, genotypes, parentIds);
		return new GenotypeImputationHMM(states, m, positions);
	}
	
	
	
	public int getStartsBaumWelch() {
		return startsBaumWelch;
	}

	public void setStartsBaumWelch(int startsBaumWelch) {
		this.startsBaumWelch = startsBaumWelch;
	}

	/**
	 * Impute the given set of genotypes assuming they are inbreds
	 * @param genotypes Map with one entry per individual. The key is the sample id and the value is a list of genotype calls
	 */
	public void imputeGenotypes (Map<String, List<CalledSNV>> genotypes, int [][][] outClusters) {
		List<String> sampleIds = new ArrayList<String>();
		sampleIds.addAll(genotypes.keySet());
		int n = sampleIds.size();
		
		int m = genotypes.values().iterator().next().size();
		int k = getNumStates();
		if(m!=getSteps()) throw new IllegalArgumentException("Number of variants: "+m+" in the set of genotypes does not coincide with steps of the HMM: "+getSteps());
		double [][][] sumAlleleProbs = new double [n][m][2];
		double [][][] nextAlleleProbs = new double [n][m][2];
		double [][] nextPosteriorsSample = new double [m][k];
		int [] nextViterbiPathSample = new int [m];
		for(int i=0;i<n;i++) {
			NumberArrays.initializeDoubleMatrix(sumAlleleProbs[i]);
		}
		for(int h=0;h<startsBaumWelch;h++) {
			NumberArrays.initializeIntMatrix(outClusters[h]);
			getLog().info("Training and sampling iteration: "+h);
			train();
			for(int i=0;i<n;i++) {
				String sampleId = sampleIds.get(i);
				List<CalledSNV> genotypesSample = genotypes.get(sampleId);
				List<Byte> haplotype = makeHaplotypeWithHomozygous(genotypesSample);
				
				//Allele posteriors for genotyping
				calculateAllelePosteriors(haplotype, nextAlleleProbs[i]);
				NumberArrays.accumulate(sumAlleleProbs[i],nextAlleleProbs[i]);
				
				//State posteriors for assignments
				calculatePosteriors(haplotype, nextPosteriorsSample);
				
				//Best viterbi path
				getViterbiPath(haplotype, nextViterbiPathSample);
				
				//Conciliate viterbi with posterior
				assignClusters (sampleId,genotypesSample, nextPosteriorsSample, nextViterbiPathSample, outClusters[h][i]);
			}
		}
		for(int i=0;i<n;i++) {
			String sampleId = sampleIds.get(i);
			getLog().info("Choosing best genotypes for sample: "+sampleId);
			imputeGenotypes(genotypes.get(sampleId),sumAlleleProbs[i],startsBaumWelch);
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
		for(int i=0;i<nextPosteriors.length;i++) {
			outClusters[i] = viterbiPath[i];
			int idxMaxPos = NumberArrays.getIndexMaximum(nextPosteriors[i]);
			CalledSNV call = genotypesSample.get(i);
			if(idxMaxPos==-1) {
				getLog().warning("Sample "+sampleId+" site "+call.getSequenceName()+":"+call.getFirst()+" Viterbi state "+viterbiPath[i]+" best posterior state not found");
			}/* else if(idxMaxPos==viterbiPath[i]) {
				if(i>0 && outClusters[i]!=outClusters[i-1]) {
					getLog().info("Crossover predicted for sample "+sampleId+" site "+call.getSequenceName()+":"+call.getFirst()+" Previous state: "+outClusters[i-1]+" new state: "+outClusters[i]);
				}
			} else {
				double pViterbiState = nextPosteriors[i][viterbiPath[i]];
				double pMaxPosterior = nextPosteriors[i][idxMaxPos];
				if(pMaxPosterior-pViterbiState>0.05) {
					getLog().warning("Sample "+sampleId+" site "+call.getSequenceName()+":"+call.getFirst()+" Viterbi state "+viterbiPath[i]+" different than best posterior state "+idxMaxPos+" viterbiProb: "+pViterbiState+"posterior prob: "+pMaxPosterior);
				}
			}*/
		}
	}

	private void imputeGenotypes(List<CalledSNV> genotypes, double[][] sums, int steps) {
		
		int imputed = 0;
		int heterozygous = 0;
		int inconsistent = 0;
		for(int i=0;i<sums.length;i++) {
			CalledSNV call = genotypes.get(i);
			byte g = call.getGenotype();
			int maxA = 0;
			for(int j=0;j<sums[i].length;j++) {
				if(sums[i][maxA]<sums[i][j]) {
					maxA = j;
				}
			}
			//Homozygous imputation
			byte maxG = (byte)(2*maxA);
			if(call.isUndecided() ) {
				//Make genotype out of maximum allele
				call.setGenotype(maxG);
				double prob = sums[i][maxA]/steps;
				call.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(1-prob));
				imputed++;
			} else if (call.isHeterozygous()) {
				heterozygous++;
			} else if(g!=maxG) {
				getLog().info("Genotype at "+call.getSequenceName()+":"+call.getFirst()+" inconsistent with prediction. Predicted: "+maxG+" given: "+g);
				inconsistent++;
			}
		}
		getLog().info("Imputed: "+imputed+" heterozygous: "+heterozygous+" inconsistent: "+inconsistent);
	}
}
