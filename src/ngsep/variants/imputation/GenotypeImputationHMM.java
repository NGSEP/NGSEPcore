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
import java.util.TreeMap;

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
	
	/**
	 * Impute the given set of genotypes assuming they are inbreds
	 * @param genotypes Map with one entry per individual. The key is the sample id and the value is a list of genotype calls
	 * @return Map<String,List<Integer>> For each individual, list of ids if the most likely HMM states for each site.
	 * 
	 */
	public Map<String,List<Integer>> imputeGenotypes (Map<String, List<CalledSNV>> genotypes) {
		List<String> sampleIds = new ArrayList<String>();
		sampleIds.addAll(genotypes.keySet());
		int n = sampleIds.size();
		
		int m = genotypes.values().iterator().next().size();
		if(m!=getSteps()) throw new IllegalArgumentException("Number of variants: "+m+" in the set of genotypes does not coincide with steps of the HMM: "+getSteps());
		double [][][] sumAlleleProbs = new double [n][m][2];
		double [][][] nextAlleleProbs = new double [n][m][2];
		double [][][] sumStateProbs = new double [n][m][getNumStates()];
		double [][][] nextStateProbs = new double [n][m][getNumStates()];
		for(int i=0;i<n;i++) {
			for(int j=0;j<m;j++) {
				Arrays.fill(sumAlleleProbs[i][j], 0.0);
				Arrays.fill(sumStateProbs[i][j], 0.0);
			}
		}
		if(getTrainingData() == null) setTrainingData(makeTrainingDataWithHomozygous(genotypes));
		for(int h=0;h<startsBaumWelch;h++) {
			getLog().info("Training and sampling iteration: "+h);
			train();
			for(int i=0;i<n;i++) {
				String sampleId = sampleIds.get(i);
				List<CalledSNV> genotypesSample = genotypes.get(sampleId);
				List<Byte> haplotype = makeHaplotypeWithHomozygous(genotypesSample);
				calculateAllelePosteriors(haplotype, nextAlleleProbs[i]);
				calculateStatePosteriors(haplotype, nextStateProbs[i]);
				NumberArrays.accumulate(sumAlleleProbs[i],nextAlleleProbs[i]);
				NumberArrays.accumulate(sumStateProbs[i],nextStateProbs[i]);
				
			}
		}
		Map<String,List<Integer>> assignments = new TreeMap<String, List<Integer>>();
		for(int i=0;i<n;i++) {
			String sampleId = sampleIds.get(i);
			getLog().info("Choosing best genotypes for sample: "+sampleId);
			imputeGenotypes(genotypes.get(sampleId),sumAlleleProbs[i],startsBaumWelch);
			assignments.put(sampleId, calculateFinalAssignments(sumStateProbs[i]));
		}
		
		return assignments;
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

	/**
	 * 
	 * @param votes Matrix with as many rows as markers and as many columns as HMM states
	 * @return List<Integer> For each site, index of the HMM state with the highest probability
	 */
	private List<Integer> calculateFinalAssignments(double[][] stateProbabilities) {
		List<Integer> assignments = new ArrayList<Integer>();
		for(int i=0;i<stateProbabilities.length;i++) {
			int maxJ = 0;
			for(int j=0;j<stateProbabilities[i].length;j++) {
				if(stateProbabilities[i][maxJ]<stateProbabilities[i][j]) {
					maxJ = j;
				}
			}
			assignments.add(maxJ);
		}
		return assignments;
	}
	
	private List<List<? extends Object>> makeTrainingDataWithHomozygous(Map<String, List<CalledSNV>> genotypes) {
		List<List<? extends Object>> haplotypes = new ArrayList<>();
		for (List<CalledSNV> genotypesSample:genotypes.values()) {
			List<? extends Object> hapSample = makeHaplotypeWithHomozygous(genotypesSample);
			haplotypes.add(hapSample);
		}
		return haplotypes;
	}

	private List<Byte> makeHaplotypeWithHomozygous(List<CalledSNV> genotypesSample) {
		List<Byte> hapSample = new ArrayList<>();
		for(CalledSNV call:genotypesSample) {
			byte allele = -1;
			if(call.isHomozygousReference()) {
				allele = 0;
			} else if (call.isHomozygous()) {
				allele = 1;
			}
			hapSample.add(allele);
		}
		return hapSample;
	}
	
	


}
