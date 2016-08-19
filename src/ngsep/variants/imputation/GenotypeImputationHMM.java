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
import java.util.logging.Logger;

import ngsep.hmm.HMMState;
import ngsep.hmm.RecombinationHMM;
import ngsep.math.LogMath;
import ngsep.math.PhredScoreHelper;
import ngsep.variants.CalledSNV;


public class GenotypeImputationHMM extends RecombinationHMM {
	public static final double DEF_PCT_MISSING = 10;
	public static final double DEF_PCT_HETERO = 1;
	
	public static final int DEF_ITER_BAUM_WELCH = 10;
	public static final int DEF_STARTS_BAUM_WELCH = 5;
	
	private Logger log = Logger.getLogger(GenotypeImputationHMM.class.getName());
	
	private double pctMissingData=DEF_PCT_MISSING;
	private double pctHeterozygous=DEF_PCT_HETERO;
	private List<Integer> positions = null;
	private Double avgCMPerKbp = null;
	private boolean updateEmissionKnownSites = false;
	private boolean fixedTransitions = false;
	private int iterationsBaumWelch = DEF_ITER_BAUM_WELCH;
	private int startsBaumWelch = DEF_STARTS_BAUM_WELCH;
	private Map<String,List<CalledSNV>> referenceHaplotypes = null;
	
	private Double logProbMissing = Math.log10(DEF_PCT_MISSING/100.0);
	private Double logProbHetero = Math.log10(DEF_PCT_HETERO/100.0);
	private Double logProbHomo = Math.log10((100.0-DEF_PCT_HETERO-DEF_PCT_MISSING)/100.0);
	
	//Local arrays to save reallocation over many runs
	private Double [][] forwardLogs=new Double[0][0];
	private Double [][] backwardLogs=new Double[0][0];
	private Double [][][] logTransitions = new Double [0][0][0];
	private Double [][][] logEmissions = new Double [0][0][0];
	private Double [] logStarts = new Double [0];
	
	
	public GenotypeImputationHMM(List<? extends HMMState> states, int numMarkers) {
		super(states, numMarkers);
	}
	
	public GenotypeImputationHMM(List<? extends HMMState> states, int numMarkers, List<Integer> positions) {
		this(states, numMarkers);
		this.positions = positions;
	}
	
	public static GenotypeImputationHMM createHMM (Map<String, List<CalledSNV>> genotypes, List<String> parentIds, int k) {
		List<CalledSNV> snvs = genotypes.values().iterator().next();
		List<Integer> positions = new ArrayList<Integer>();
		for(CalledSNV snv:snvs) positions.add(snv.getFirst());
		int m = snvs.size();
		List<SNPHaplotypeFounderHMMState> states = createHMMStates(genotypes,m, k,parentIds);
		//TODO: Calculate percentage of missing and heterozygous data
		return new GenotypeImputationHMM(states, m, positions);
	}
	private static  List<SNPHaplotypeFounderHMMState> createHMMStates(Map<String, List<CalledSNV>> genotypes, int m, int k,List<String> parentIds) {
		List<SNPHaplotypeFounderHMMState> states = new ArrayList<SNPHaplotypeFounderHMMState>();
		for(String parentId:parentIds) {
			SNPHaplotypeFounderHMMState state = new SNPHaplotypeFounderHMMState(genotypes.get(parentId));
			state.setId(parentId);
			states.add(state);
		}
		while(states.size()<k) {
			SNPHaplotypeFounderHMMState state = new SNPHaplotypeFounderHMMState(m);
			states.add(state);
		}
		return states;
	}
	public Logger getLog() {
		return log;
	}
	
	public void setLog(Logger log) {
		this.log = log;
	}


	public double getPctMissingData() {
		return pctMissingData;
	}

	public void setPctMissingData(double pctMissingData) {
		this.pctMissingData = pctMissingData;
		logProbMissing = LogMath.log10(pctMissingData/100.0);
		updateLogProbHomo();
		
	}

	public double getPctHeterozygous() {
		return pctHeterozygous;
	}
	
	public void setPctHeterozygous(double pctHeterozygous) {
		this.pctHeterozygous = pctHeterozygous;
		logProbHetero = LogMath.log10(pctHeterozygous/100.0);
		updateLogProbHomo();
	}
	
	public double getAvgCMPerKbp() {
		return avgCMPerKbp;
	}

	public void setAvgCMPerKbp(double avgCMPerKbp) {
		this.avgCMPerKbp = avgCMPerKbp;
	}

	public boolean isUpdateEmissionKnownSites() {
		return updateEmissionKnownSites;
	}

	public void setUpdateEmissionKnownSites(boolean updateEmissionKnownSites) {
		this.updateEmissionKnownSites = updateEmissionKnownSites;
	}
	
	public boolean isFixedTransitions() {
		return fixedTransitions;
	}

	public void setFixedTransitions(boolean fixedTransitions) {
		this.fixedTransitions = fixedTransitions;
	}

	private void updateLogProbHomo() {
		double probHomo = (100.0 - pctMissingData - pctHeterozygous)/100;
		logProbHomo = LogMath.log10(probHomo);
	}



	@Override
	public Double getEmission(int state, Object value, int step) {
		Double logCondHomozygous = super.getEmission(state, value, step);
		byte genotype = SNPHaplotypeFounderHMMState.getGenotype(value);
		if(genotype == CalledSNV.GENOTYPE_UNDECIDED) {
			return logProbMissing;
		} else if (genotype == CalledSNV.GENOTYPE_HETERO) {
			return logProbHetero;
		}
		return LogMath.logProduct(logCondHomozygous, logProbHomo);	
	}
	/**
	 * Impute the given set of genotypes
	 * @param genotypes Map with one entry per individual. The key is the sample id and the value is a list of genotype calls
	 * @return Map<String,List<Integer>> If the samples are haploid, it returns the list of ids if the most likely HMM states for
	 * each site. If the samples are diploid, it just returns an empty map.
	 * TODO: Improve return type for diploids
	 */
	public Map<String,List<Integer>> imputeGenotypes (Map<String, List<CalledSNV>> genotypes) {
		List<String> sampleIds = new ArrayList<String>();
		sampleIds.addAll(genotypes.keySet());
		int n = sampleIds.size();
		
		int m = genotypes.values().iterator().next().size();
		if(m!=getSteps()) throw new IllegalArgumentException("Number of variants: "+m+" in the set of genotypes does not coincide with steps of the HMM: "+getSteps());
		double [][][] sumAlleleProbs = new double [n][m][3];
		double [][][] nextAlleleProbs = new double [n][m][3];
		//TODO: Check memory and use for diploid imputation
		double [][][] sumStateProbs = new double [n][m][getNumStates()];
		double [][][] nextStateProbs = new double [n][m][getNumStates()];
		for(int i=0;i<n;i++) {
			for(int j=0;j<m;j++) {
				Arrays.fill(sumAlleleProbs[i][j], 0.0);
				Arrays.fill(sumStateProbs[i][j], 0.0);
			}
		}
		for(int h=0;h<startsBaumWelch;h++) {
			log.info("Training and sampling iteration: "+h);
			if(referenceHaplotypes!=null) train(referenceHaplotypes,true);
			else train(genotypes,false);
			for(int i=0;i<n;i++) {
				String sampleId = sampleIds.get(i);
				List<CalledSNV> genotypesSample = genotypes.get(sampleId);
				calculatePosteriors(genotypesSample, nextAlleleProbs[i]);
				calculateAssignments(genotypesSample, nextStateProbs[i]);
				accumulate(sumAlleleProbs[i],nextAlleleProbs[i]);
				accumulate(sumStateProbs[i],nextStateProbs[i]);
				
			}
		}
		Map<String,List<Integer>> assignments = new TreeMap<String, List<Integer>>();
		for(int i=0;i<n;i++) {
			String sampleId = sampleIds.get(i);
			log.info("Choosing best genotypes for sample: "+sampleId);
			imputeGenotypes(genotypes.get(sampleId),sumAlleleProbs[i],startsBaumWelch);
			assignments.put(sampleId, calculateFinalAssignments(sumStateProbs[i]));
		}
		
		return assignments;
	}
	
	public void train(Map<String, List<CalledSNV>> genotypes, boolean forceHaploid) {
		log.info("Training model with "+genotypes.size()+" sequences");
		int n = getNumStates();
		setRandomTransitions();
		//printTransitions(0);
		//printTransitions(2);
		double logUniformStart = Math.log10(1.0/(double)n);
		for(int j=0;j<n;j++) {
			SNPHaplotypeFounderHMMState state = (SNPHaplotypeFounderHMMState) getState(j);
			state.setLogStart(logUniformStart);
			state.setRandomEmissions(updateEmissionKnownSites);
		}
		for(int h = 0; h < iterationsBaumWelch; h++) {
			log.info("Running "+h+" Baum-Welch iteration");
			runBaumWelchStep(genotypes,forceHaploid);
		}
		//printTransitions(0);
		//printTransitions(2);
	}
	
	public void printTransitions(int step) {
		int n = getNumStates();
		System.out.println("Transitions step: "+step);
		for(int j=0;j<n;j++) {
			System.out.print(getTransition(j, 0, step));
			for (int k=1;k<n;k++) {
				System.out.print("\t"+getTransition(j, k, step));
			}
			System.out.println();
		}
		
	}

	public void setRandomTransitions() {
		if(positions!=null && avgCMPerKbp!=null) {
			calculateTransitions(positions, avgCMPerKbp,!fixedTransitions);
		} else {
			super.setRandomTransitions();
		}
	}
	protected void runBaumWelchStep(Map<String, List<CalledSNV>> haplotypes, boolean forceHaploid) {
		initArrays(getNumStates(), getSteps());
		Arrays.fill(logStarts, 0.0);
		for(int i=0;i<logTransitions.length;i++) {
			for(int j=0;j<logTransitions[i].length;j++) {
				Arrays.fill(logTransitions[i][j], null);
			}
		}
		for(int j=0;j<logEmissions.length;j++) {
			for(int i=0;i<logEmissions[j].length;i++) {
				logEmissions[j][i][0] = logEmissions[j][i][1] = null;
			}
		}
		for (List<CalledSNV> haplotypesSample:haplotypes.values()) {
			Double logProb = calculateForward(haplotypesSample, forwardLogs);
			calculateBackward(haplotypesSample, backwardLogs);
			//Calculate new starts
			for(int j=0;j<logStarts.length;j++) {
				Object o = haplotypesSample.get(0);
				Double seqProduct = LogMath.logProduct(forwardLogs[0][j], backwardLogs[0][j]);
				seqProduct = LogMath.logProduct(seqProduct, getEmission(j, o, 0));
				seqProduct = LogMath.logProduct(seqProduct, -logProb);
				logStarts[j] = LogMath.logSum(logStarts[j], seqProduct);
			}
			//Calculate new transitions
			if(!fixedTransitions) {
				for(int i=0;i<logTransitions.length;i++) {
					for(int j=0;j<logTransitions[i].length;j++) {
						
						for(int k=0;k<logTransitions[i][j].length;k++) {
							Object o1 = haplotypesSample.get(i);
							Object o2 = haplotypesSample.get(i+1);
							
							Double seqProduct = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i+1][k]);
							seqProduct = LogMath.logProduct(seqProduct, getEmission(j, o1, i));
							seqProduct = LogMath.logProduct(seqProduct, getEmission(k, o2, i+1));
							seqProduct = LogMath.logProduct(seqProduct, getTransition(j, k, i));
							seqProduct = LogMath.logProduct(seqProduct, -logProb);
							logTransitions[i][j][k] = LogMath.logSum(logTransitions[i][j][k], seqProduct);
						}
					}
				}
			}
			
			//Calculate new emissions
			for(int j=0;j<logEmissions.length;j++) {
				for(int i=0;i<logEmissions[j].length;i++) {
					Object o = haplotypesSample.get(i);
					byte genotype = SNPHaplotypeFounderHMMState.getGenotype(o);
					if(genotype == CalledSNV.GENOTYPE_HOMOREF || genotype == CalledSNV.GENOTYPE_HOMOALT) {
						Double seqProduct = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i][j]);
						seqProduct = LogMath.logProduct(seqProduct, getEmission(j, o, i));
						seqProduct = LogMath.logProduct(seqProduct, -logProb);
						if(genotype == CalledSNV.GENOTYPE_HOMOREF) {
							logEmissions[j][i][0] = LogMath.logSum(logEmissions[j][i][0], seqProduct);
						} else {
							logEmissions[j][i][1] = LogMath.logSum(logEmissions[j][i][1], seqProduct);
						}
						
					}
				}
			}
			
		}
		//Normalize and update starts
		Double total = null;
		for(int j=0;j<logStarts.length;j++) total = LogMath.logSum(total, logStarts[j]);
		for(int j=0;j<logStarts.length;j++) getState(j).setLogStart(LogMath.logProduct(logStarts[j],-total));
		//Normalize and update transitions
		if(!fixedTransitions) {
			for(int i=0;i<logTransitions.length;i++) {
				setTransitions(logTransitions[i], i);
			}
		}
		
		//Normalize and update emissions
		for(int j=0;j<logEmissions.length;j++) {
			//TODO: Replace direct class casting with an interface
			SNPHaplotypeFounderHMMState state = (SNPHaplotypeFounderHMMState)getState(j);
			state.setEmissionLogProbs(logEmissions[j], updateEmissionKnownSites);
		}
	}
	
	private void initArrays(int k, int m) {
		if(forwardLogs.length!=m || forwardLogs[0].length!=k) forwardLogs = new Double[m][k];
		if(backwardLogs.length!=m || backwardLogs[0].length!=k) backwardLogs = new Double[m][k];
		if(logTransitions.length!=m-1 || logTransitions[0].length!=k) logTransitions = new Double [m-1][k][k];
		if(logEmissions.length!=k || logEmissions[0].length!=m) logEmissions = new Double [k][m][2];
		if(logStarts.length!=k) logStarts = new Double [k];
	}

	protected void calculatePosteriors(List<CalledSNV> genotypes, double[][] genotypePosteriors) {
		Byte b0 = CalledSNV.GENOTYPE_HOMOREF;
		Byte b1 = CalledSNV.GENOTYPE_HOMOALT;
		
		int m = genotypes.size();
		int n = getNumStates();
		calculateForward(genotypes, forwardLogs);
		calculateBackward(genotypes, backwardLogs);
		Double [] stateLogPosteriors = new Double [n];
		for(int i=0;i<m;i++) {
			CalledSNV csnv = genotypes.get(i);
			byte g = csnv.getGenotype();
			Double log0 = null;
			Double log1 = null;
			for(int j=0;j<n;j++) {
				Double fTimesB = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i][j]);
				log0 = LogMath.logSum(log0, LogMath.logProduct(fTimesB, getEmission(j, b0, i)));
				log1 = LogMath.logSum(log1, LogMath.logProduct(fTimesB, getEmission(j, b1, i)));
				stateLogPosteriors[j] = LogMath.logProduct(fTimesB, getEmission(j, g, i));
			}
			//Normalize and raise to calculate final probabilities of genotypes
			Double logSum = LogMath.logSum(log0, log1);
			double prob0 = LogMath.power10(LogMath.logProduct(log0, -logSum));
			double prob1 = LogMath.power10(LogMath.logProduct(log1, -logSum));
			double sum = prob0 + prob1;
			prob0/=sum;
			prob1/=sum;
			genotypePosteriors[i][CalledSNV.GENOTYPE_HOMOREF] = prob0;
			genotypePosteriors[i][CalledSNV.GENOTYPE_HETERO] = 0;
			genotypePosteriors[i][CalledSNV.GENOTYPE_HOMOALT] = prob1;
		}
	}
	
	private void calculateAssignments(List<CalledSNV> genotypes, double[][] statePosteriors) {
		int m = genotypes.size();
		int n = getNumStates();
		calculateForward(genotypes, forwardLogs);
		calculateBackward(genotypes, backwardLogs);
		Double [] stateLogPosteriors = new Double [n];
		for(int i=0;i<m;i++) {
			CalledSNV csnv = genotypes.get(i);
			byte g = csnv.getGenotype();
			for(int j=0;j<n;j++) {
				Double fTimesB = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i][j]);
				stateLogPosteriors[j] = LogMath.logProduct(fTimesB, getEmission(j, g, i));
			}
			LogMath.normalizeLogs(stateLogPosteriors);
			for(int j=0;j<n;j++) {
				statePosteriors[i][j] = LogMath.power10(stateLogPosteriors[j]);
			}
		}
	}

	private void accumulate(double[][] cumulative, double[][] next) {
		for(int i=0;i<cumulative.length;i++) {
			for(int j=0;j<cumulative[i].length;j++) {
				cumulative[i][j]+=next[i][j];
			}
		}
		
	}

	private void imputeGenotypes(List<CalledSNV> genotypes, double[][] sums, int steps) {
		
		int imputed = 0;
		int heterozygous = 0;
		int inconsistent = 0;
		for(int i=0;i<sums.length;i++) {
			CalledSNV call = genotypes.get(i);
			byte g = call.getGenotype();
			int maxG = 0;
			for(int j=0;j<sums[i].length;j++) {
				if(sums[i][maxG]<sums[i][j]) {
					maxG = j;
				}
			}
			
			if(call.isUndecided() ) {
				call.setGenotype((byte)maxG);
				double prob = sums[i][maxG]/steps;
				call.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(1-prob));
				imputed++;
			} else if (call.isHeterozygous()) {
				heterozygous++;
			} else if(g!=maxG) {
				log.info("Genotype at "+call.getSequenceName()+":"+call.getFirst()+" inconsistent with prediction. Predicted: "+maxG+" given: "+g);
				inconsistent++;
			}
		}
		log.info("Imputed: "+imputed+" heterozygous: "+heterozygous+" inconsistent: "+inconsistent);
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
}
