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
package ngsep.discovery.rd;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.hmm.ConstantTransitionHMM;
import ngsep.hmm.HMM;
import ngsep.hmm.HMMState;
import ngsep.math.LogMath;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledCNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;


public abstract class AbstractHMMReadDepthAlgorithm implements SingleSampleReadDepthAlgorithm {

	private Logger log = Logger.getLogger(AbstractHMMReadDepthAlgorithm.class.getName());
	private ReadDepthDistribution readDepthDistribution;
	
	private byte normalPloidy = 2;
	private double changeProbability = 0.01;
	
	
	@Override
	public void setLog(Logger log) {
		this.log = log;
	}

	@Override
	public void setGenomeSize(long genomeSize) {
		//Parameter not needed in this implementation

	}

	@Override
	public void setNormalPloidy(byte normalPloidy) {
		this.normalPloidy = normalPloidy;

	}

	@Override
	public void setReadDepthDistribution(ReadDepthDistribution distribution) {
		this.readDepthDistribution = distribution;
	}

	@Override
	public List<CalledCNV> callCNVs() {
		List<CalledCNV> answer = new ArrayList<CalledCNV>();
		List<String> seqNames = readDepthDistribution.getSequences().getNamesStringList();
		log.info("Building HMM");
		HMM hmm = buildHMM ();
		for(String seqName:seqNames) {
			log.info("Calling CNVs for sequence "+seqName);
			List<ReadDepthBin> seqBins = readDepthDistribution.getBins(seqName);
			List<CalledCNV> cnvsSeq = callCNVsSequence(seqName,seqBins,hmm); 
			log.info("Called "+cnvsSeq.size()+" CNVs for sequence "+seqName);
			answer.addAll(cnvsSeq);
		}
		return answer;
	}
	
	
	public ReadDepthDistribution getReadDepthDistribution() {
		return readDepthDistribution;
	}

	public byte getNormalPloidy() {
		return normalPloidy;
	}
	
	

	public double getChangeProbability() {
		return changeProbability;
	}

	public void setChangeProbability(double changeProbability) {
		this.changeProbability = changeProbability;
	}

	protected abstract String getSource();
	//TODO: Design better
	protected abstract HMMState createHMMState(int copies, Double logStart);

	public Logger getLog() {
		return log;
	}

	private HMM buildHMM() {
		int nStates = 4*normalPloidy+1;
		
		List<HMMState> states = new ArrayList<HMMState>(nStates);
		double randomLogStart = LogMath.log10(1.0/nStates);
		for(int i=0;i<nStates;i++) {
			HMMState state = createHMMState(i, randomLogStart);
			states.add(state);
			
		}
		ConstantTransitionHMM answer = new ConstantTransitionHMM(states);
		answer.calculateUniformChangeTransitions(changeProbability);
		return answer;
	}

	
	

	private List<CalledCNV> callCNVsSequence(String seqName, List<ReadDepthBin> seqBins, HMM hmm) {
		int m = seqBins.size();
		int n = hmm.getNumStates();
		Double [] [] posteriorLogs = new Double [m][n];
		List<Double> observations = buildObservations(seqBins);
		hmm.calculatePosteriorLogs(observations, posteriorLogs);
		//hmm.calculateBackward(observations, posteriorLogs);
		if("chrI".equals(seqName)) printLogProbs(posteriorLogs);
		List<CalledCNV> answer = new ArrayList<CalledCNV>();
		int nextStartBin = -1;
		int copies = normalPloidy;
		for(int i=0;i<m;i++) {
			int state = chooseState(posteriorLogs[i]);
			if(state!=copies) {
				if(copies!=normalPloidy) {
					answer.add(createCNV(seqName,seqBins,posteriorLogs,nextStartBin,i-1,copies));
				}
				copies = state;
				nextStartBin=i;
			}
		}
		return answer;
	}

	private void printLogProbs(Double[][] logProbs) {
		for(int i=0;i<logProbs.length;i++) {
			for(int j=0;j<logProbs[i].length;j++) {
				System.out.print(" "+logProbs[i][j]);
			}
			System.out.println();
		}
	}

	private List<Double> buildObservations(List<ReadDepthBin> seqBins) {
		List<Double> observations = new ArrayList<Double>();
		for(ReadDepthBin bin:seqBins) {
			observations.add(bin.getCorrectedReadDepth());
		}
		return observations;
	}
	
	private int chooseState(Double[] logProbs) {
		int maxI = normalPloidy;
		double maxVal = Double.MIN_VALUE;
		if(logProbs[normalPloidy]!=null) maxVal = logProbs[normalPloidy];
		for(int i=0;i<logProbs.length;i++) {
			if(logProbs[i]!=null && maxVal<logProbs[i]) {
					maxI = i;
					maxVal = logProbs[i];
			}
		}
		return maxI;
	}

	private CalledCNV createCNV(String seqName, List<ReadDepthBin> seqBins, Double [][] posteriorLogs, int firstI, int lastI, int copies) {
		ReadDepthBin firstBin = seqBins.get(firstI);
		ReadDepthBin lastBin = seqBins.get(lastI);
		int bins = 0;
		int fragments = 0;
		double avgProb = 0;
		for(int i=firstI;i<=lastI;i++) {
			fragments+=seqBins.get(i).getRawReadDepth();
			avgProb += LogMath.power10(posteriorLogs[i][copies]); 
			bins++;
		}
		avgProb/=bins;
		GenomicVariantImpl cnv = new GenomicVariantImpl(seqName, firstBin.getFirst(), lastBin.getLast(),GenomicVariant.TYPE_CNV);
		CalledCNV call = new CalledCNV(cnv,copies);
		call.setTotalReadDepth(fragments);
		call.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(1-avgProb));
		call.setSource(getSource());
		return call;
	}
}