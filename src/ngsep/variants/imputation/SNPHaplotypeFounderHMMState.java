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

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import ngsep.hmm.HMMState;
import ngsep.math.LogMath;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledSNV;


public class SNPHaplotypeFounderHMMState implements HMMState {
	private String id = null;
	private byte [] haplotype = new byte [0]; //-1 for undecided, 0 for allele zero, 1 for allele 1
	private Double [] allele0Logs = new Double [0];
	private Double [] allele1Logs = new Double [0];
	private Double logStart=null;
	
	public SNPHaplotypeFounderHMMState(int haplotypeLength) {
		initArrays(haplotypeLength);
	}
	public SNPHaplotypeFounderHMMState(List<CalledSNV> calls) {
		int m = calls.size();
		initArrays(m);
		for(int i=0;i<m;i++) {
			CalledSNV call = calls.get(i);
			double successProb = 1.0 - PhredScoreHelper.calculateProbability(call.getGenotypeQuality());
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

	private void initArrays(int m) {
		if(haplotype.length!=m) {
			haplotype = new byte [m];
			allele0Logs = new Double [m];
			allele1Logs = new Double [m];
		}
		Arrays.fill(haplotype, (byte)-1);
		setRandomEmissions(true);
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
		byte b = getGenotype(value);
		//TODO: take into account genotype quality
		//System.out.println("Genotype: "+b+" allele 0 log: "+allele0Logs[step]+" allele 1 log: "+allele1Logs[step]);
		if(b==CalledSNV.GENOTYPE_HOMOREF) return allele0Logs[step];
		if(b==CalledSNV.GENOTYPE_HOMOALT) return allele1Logs[step];
		return null;
	}

	public static byte getGenotype(Object value) {
		byte b = -1;
		if (value instanceof Byte) {
			b = (Byte)value;
		} else if (value instanceof CalledSNV) {
			b = ((CalledSNV)value).getGenotype();
		}
		return b;
	}


}
