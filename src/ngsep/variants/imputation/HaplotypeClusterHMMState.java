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
import java.util.Random;

import ngsep.hmm.HMMState;
import ngsep.math.LogMath;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledSNV;


public class HaplotypeClusterHMMState implements HMMState {
	private String id = null;
	private byte [] haplotype = new byte [0]; //-1 for undecided, 0 for allele zero, 1 for allele 1
	private Double [] allele0Logs = new Double [0];
	private Double [] allele1Logs = new Double [0];
	private Double logStart=null;
	
	public static final Double LOGPROB_UNEXPECTED = Math.log10(0.01);
	public static final Double LOGPROB_EXPECTED = Math.log10(0.99);
	
	public HaplotypeClusterHMMState(int haplotypeLength) {
		initArrays(haplotypeLength);
	}
	
	private void initArrays(int m) {
		if(haplotype.length!=m) {
			haplotype = new byte [m];
			allele0Logs = new Double [m];
			allele1Logs = new Double [m];
		}
		Arrays.fill(haplotype, (byte)-1);
		setRandomEmissions(true);
	}
	public void fillEmissionsWithGenotypes(String id, List<CalledSNV> calls) {
		int m = calls.size();
		this.id = id;
		for(int i=0;i<m;i++) {
			CalledSNV call = calls.get(i);
			double successProb = 0.99;
			if(call.getGenotypeQuality()>0) successProb = 1.0 - PhredScoreHelper.calculateProbability(call.getGenotypeQuality());
			//TODO: Improve handling
			if(successProb > 0.999) successProb = 0.999;
			byte g = call.getGenotype();
			Double logError = LogMath.log10(1.0-successProb);
			Double logNoError = LogMath.log10(successProb);
			if(g==CalledSNV.GENOTYPE_HOMOREF) {
				haplotype[i] = 0;
				allele0Logs [i] = logNoError;
				allele1Logs [i] = logError;
			}
			if(g==CalledSNV.GENOTYPE_HOMOALT) {
				haplotype[i] = 1;
				allele0Logs [i] = logError;
				allele1Logs [i] = logNoError;
			}
		}
	}
	
	public void fillEmissionsWithHaplotype(String id, List<CalledSNV> calls, int hapId) {
		int m = calls.size();
		this.id = id;
		for(int i=0;i<m;i++) {
			CalledSNV call = calls.get(i);
			if(!call.isPhased()) continue;
			if(call.getCopyNumber()<=hapId) continue;
			//TODO: Read haplotype qualities
			double successProb = 0.99;
			if(call.getGenotypeQuality()>0) successProb = 1.0 - PhredScoreHelper.calculateProbability(call.getGenotypeQuality());
			//TODO: Improve handling
			if(successProb > 0.999) successProb = 0.999;
			byte [] idsPhasedAlleles = call.getIndexesPhasedAlleles();
			byte phasedAllele = idsPhasedAlleles[hapId];
			Double logError = LogMath.log10(1.0-successProb);
			Double logNoError = LogMath.log10(successProb);
			if(phasedAllele == 0) {
				haplotype[i] = 0;
				allele0Logs [i] = logNoError;
				allele1Logs [i] = logError;
			}
			else {
				haplotype[i] = 1;
				allele0Logs [i] = logError;
				allele1Logs [i] = logNoError;
			}
		}
	}

	

	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public void setRandomEmissions(boolean updateKnownSites) {
		Random r = new Random();
		for(int i=0;i<haplotype.length;i++) {
			if(updateKnownSites || haplotype[i]==-1) {
				double d = r.nextDouble()*0.8 + 0.1;
				allele0Logs [i] = LogMath.log10(1.0-d);
				allele1Logs [i] = LogMath.log10(d);
			}
			
		}
	}
	@Override
	public Double getLogStart() {
		return logStart;
	}
	
	@Override
	public void setLogStart(Double logStart) {
		this.logStart = logStart;
	}
	/**
	 * Changes the allele probabilities. Useful method for HMM training 
	 * @param logProbs Matrix with as many rows as sites and with two columns, one for allele zero and another for allele 1
	 * @param updateKnownSites True if probabilities should be updated for sites in which the
	 * haplotype was provided as an input
	 */
	public void setEmissionLogProbs(Double [][] logProbs, boolean updateKnownSites) {
		for(int i=0;i<logProbs.length;i++) {
			if(updateKnownSites || haplotype[i]==-1) {
				Double sum = LogMath.logSum(logProbs[i][0], logProbs[i][1]);
				if(sum!=null) {
					//sum == null implies that the expected counts for both allele are equal to zero, 
					//so the genotype information can not be used to reestimate emissions
					allele0Logs[i] = LogMath.logProduct(logProbs[i][0], -sum);
					allele1Logs[i] = LogMath.logProduct(logProbs[i][1], -sum);
				}
			}
		}
	}

	@Override
	public Double getEmission(Object value, int step) {
		if(value == null) return null;
		byte b = (byte) value;
		//TODO: take into account genotype quality
		//System.out.println("Allele: "+b+" allele 0 log: "+allele0Logs[step]+" allele 1 log: "+allele1Logs[step]);
		Double answer = null;
		if(b==0) answer = allele0Logs[step];
		else if(b==1) answer = allele1Logs[step];
		if(answer == null) {
			answer = HaplotypeClusterHMMState.LOGPROB_UNEXPECTED;
		} else {
			answer = LogMath.logProduct(answer, HaplotypeClusterHMMState.LOGPROB_EXPECTED);
		}
		return answer;
	}

	public static List<HaplotypeClusterHMMState> createEmptyStates(int m, int k) {
		List<HaplotypeClusterHMMState> states = new ArrayList<HaplotypeClusterHMMState>();
		for(int i=0;i<k;i++) {
			HaplotypeClusterHMMState state = new HaplotypeClusterHMMState(m);
			states.add(state);
		}
		return states;
	}
	
	public static void fillHMMStatesInbred(List<HaplotypeClusterHMMState> states, Map<String, List<CalledSNV>> genotypes, List<String> parentIds) {
		int i=0;
		for(String parentId:parentIds) {
			HaplotypeClusterHMMState state = states.get(i);
			state.fillEmissionsWithGenotypes(parentId, genotypes.get(parentId));
			i++;
		}
	}

	
	
	public static void fillHMMStatesPhasedDiploid(List<HaplotypeClusterHMMState> states, Map<String, List<CalledSNV>> genotypes, List<String> parentIds) {
		int i=0;
		for(String parentId:parentIds) {
			List<CalledSNV> genotypesParent = genotypes.get(parentId);
			HaplotypeClusterHMMState state0 = states.get(i);
			state0.fillEmissionsWithHaplotype(parentId+"_0", genotypesParent, 0);
			i++;
			HaplotypeClusterHMMState state1 = states.get(i);
			state1.fillEmissionsWithHaplotype(parentId+"_1", genotypesParent, 1);
			i++;
		}
	}


}
