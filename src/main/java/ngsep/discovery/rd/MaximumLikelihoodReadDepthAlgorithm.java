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
import JSci.maths.statistics.NormalDistribution;
import ngsep.hmm.ConstantTransitionHMM;
import ngsep.hmm.HMM;
import ngsep.hmm.HMMState;
import ngsep.math.LogMath;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledCNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;

/**
 * 
 * @author Laura Castro
 *
 */
public class MaximumLikelihoodReadDepthAlgorithm implements SingleSampleReadDepthAlgorithm{

	public static final String SOURCE_MAXIMUMLIKELIHOOD = "MAXIMUMLIKELIHOOD";

	private Logger log = Logger.getLogger(SingleSampleReadDepthAlgorithm.class.getName());

	private ReadDepthDistribution readDepthDistribution;

	private byte normalPloidy = 2;
	private double changeProbability = 0.01;
	private HMM hmm;

	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	protected String getSource() {
		return SOURCE_MAXIMUMLIKELIHOOD;
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

	public void setNormalPloidy(byte normalPloidy) {
		this.normalPloidy = normalPloidy;

	}

	public void setReadDepthDistribution(ReadDepthDistribution distribution) {
		this.readDepthDistribution = distribution;
	}


	public List<CalledCNV> callCNVs() {
		List<CalledCNV> answer = new ArrayList<CalledCNV>();
		List<String> seqNames = readDepthDistribution.getSequences().getNamesStringList();
		log.info("Building HMM");
		buildHMM();
		for(String seqName:seqNames) {
			log.info("Calling CNVs for sequence "+seqName);
			List<ReadDepthBin> seqBins = readDepthDistribution.getBins(seqName);
			List<CalledCNV> cnvsSeq = callCNVsSequence(seqName,seqBins); 
			log.info("Called "+cnvsSeq.size()+" CNVs for sequence "+seqName);
			answer.addAll(cnvsSeq);
		}
		return answer;
	}

	private void buildHMM() {
		int nStates = 4*normalPloidy+1;

		List<HMMState> states = new ArrayList<HMMState>(nStates);
		double randomLogStart = LogMath.log10(1.0/nStates);
		for(int i=0;i<nStates;i++) {
			HMMState state = createHMMState(i, randomLogStart);
			states.add(state);

		}
		hmm = new ConstantTransitionHMM(states);
		((ConstantTransitionHMM) hmm).calculateUniformChangeTransitions(changeProbability);
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

	private CalledCNV createCNV(String seqName, List<ReadDepthBin> seqBins, Double [][] likelihoods, int firstI, int lastI, int copies) {
		ReadDepthBin firstBin = seqBins.get(firstI);
		ReadDepthBin lastBin = seqBins.get(lastI);
		int fragments = 0;
		double maxProb = 0;
		for(int i=firstI;i<=lastI;i++) {
			fragments+=seqBins.get(i).getRawReadDepth();
			Double logLike = likelihoods[i][copies];
			Double logNormalPloidy = likelihoods[i][normalPloidy];
			Double sum = LogMath.logSum(logLike, logNormalPloidy);
			double nextProb = LogMath.power10(logLike-sum);
			//System.out.println("--nextProb---" + nextProb + " ----binI---- " + binI + " -----binInormalPloidy---- " + binInormalPloidy);
			if(nextProb > maxProb) maxProb = nextProb;
		}
		GenomicVariantImpl cnv = new GenomicVariantImpl(seqName, firstBin.getFirst(), lastBin.getLast(),GenomicVariant.TYPE_CNV);
		CalledCNV call = new CalledCNV(cnv,copies);
		call.setTotalReadDepth(fragments);
		//System.out.println("-----Phred-----" + PhredScoreHelper.calculatePhredScore(1-maxProb) + " --maxProb---" + maxProb);
		call.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(1-maxProb));
		//System.out.println("---getGenoTypeQuality---" + call.getGenotypeQuality());
		call.setSource(getSource());
		return call;
	}

	protected HMMState createHMMState(int copies, Double logStart) {
		double avgNormalDepth = this.getReadDepthDistribution().getMeanReadDepth();
		double avgDepthState = avgNormalDepth*copies/getNormalPloidy();
		double varianza = Math.pow(this.getReadDepthDistribution().getSigmaReadDepth(),2);
		if(copies==0) avgDepthState = 1;
		HMMState state = new MaximumLikelihoodState(copies, avgDepthState, varianza, logStart);
		//System.out.println("Created state "+state.getId()+" with average depth "+avgDepthState+" log start "+logStart);
		return state; 
	}
	
	private void calculateLikelihood(List<Double> observations, Double[][] likelihoods){
		int m = observations.size();
		int k = hmm.getNumStates();
		if(likelihoods.length!=m) throw new IllegalArgumentException("Invalid rows of posterior logs. Expected: "+m+" Given: "+likelihoods.length);
		if(m>0 && likelihoods[0].length!=k) throw new IllegalArgumentException("Invalid columns of posterior logs. Expected: "+k+" Given: "+likelihoods[0].length);
		for(int i=0;i<likelihoods.length;i++) {
			for(int j=0;j<likelihoods[0].length;j++) {
				Double e = getEmission(j, observations.get(i));
				likelihoods[i][j] = e;
			}
		}
	}

	private List<CalledCNV> callCNVsSequence(String seqName, List<ReadDepthBin> seqBins){
		int m = seqBins.size();
		int n = hmm.getNumStates();
		Double [] [] likelihoods = new Double [m][n];
		List<Double> observations = buildObservations(seqBins);
		calculateLikelihood(observations, likelihoods);
		List<CalledCNV> answer = new ArrayList<CalledCNV>();
		int nextStartBin = -1;
		int copies = normalPloidy;
		for(int i=0;i<m;i++) {
			int state = chooseState(likelihoods[i]);
			if(state!=copies) {
				if(copies!=normalPloidy) {
					answer.add(createCNV(seqName,seqBins,likelihoods,nextStartBin,i-1,copies));
				}
				copies = state;
				nextStartBin=i;
			}
		}
		if(copies!=normalPloidy) {
			answer.add(createCNV(seqName,seqBins,likelihoods,nextStartBin,m-1,copies));
		}
		return answer;

	}

	public Double getEmission(int state, Object value) {
		return ((MaximumLikelihoodState) hmm.getState(state)).getEmission2(value);
	}

	@Override
	public void setGenomeSize(long genomeSize) {
		// TODO Auto-generated method stub
		
	}

}
class MaximumLikelihoodState implements HMMState{

	private int copies;
	private double averageDepth;
	private double variance;


	/**
	 * @param copies
	 * @param averageDepth
	 * @param logStart
	 */
	public MaximumLikelihoodState(int copies, double averageDepth, double variance, Double logStart) {
		super();
		this.copies = copies;
		this.averageDepth = averageDepth;
		this.variance = variance;
	}
	
	public Double getEmission2(Object value) {
		if(value == null || !(value instanceof Double)) return null;
		double depth = (Double)value;
		if(depth<1) depth = 1;
		NormalDistribution dist = new NormalDistribution(averageDepth,variance);
		double p = dist.cumulative(depth+0.5)-dist.cumulative(depth-0.5);
		//if(copies==0 && p<0.00001) System.out.println("Emission prob "+p+" cumulative 1: "+dist.cumulative(depth-0.5)+"cumulative 2 "+dist.cumulative(depth+0.5)+" depth "+depth);
		return LogMath.log10(p);
	}
	
	public double getAverageDepth(){
		return averageDepth;
	}

	public String getId() {
		return ""+copies;
	}

	@Override
	public Double getEmission(Object value, int step) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Double getLogStart() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setLogStart(Double logStart) {
		// TODO Auto-generated method stub
		
	}

}