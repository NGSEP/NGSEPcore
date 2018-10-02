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
package ngsep.hmm;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import ngsep.math.LogMath;


public abstract class AbstractHMM implements HMM {
	
	public static final int DEF_STARTS_BAUM_WELCH = 5;
	public static final int DEF_ITER_BAUM_WELCH = 20;
	
	private Logger log = Logger.getLogger(AbstractHMM.class.getName());
	private Double [][] forwardLogs=new Double[0][0];
	private Double [][] backwardLogs=new Double[0][0];
	private Double [][] posteriorLogs=new Double[0][0];
	private Double [][] viterbiLogs = new Double [0][0];
	private int [][] viterbiBacktrace = new int [0][0];
	
	public Logger getLog() {
		return log;
	}
	
	public void setLog(Logger log) {
		this.log = log;
	}
	
	@Override
	public Double getEmission(int state, Object value, int step) {
		return getState(state).getEmission(value,step);
	}
	
	@Override
	public Double getStart(int state) {
		return getState(state).getLogStart();
	}
	@Override
	public Double calculatePosteriorLogs(List<? extends Object> observations,Double[][] posteriorLogs) {
		int m = observations.size();
		int k = getNumStates();
		initArrays(m,k,false);
		if(posteriorLogs.length!=m) throw new IllegalArgumentException("Invalid rows of posterior logs. Expected: "+m+" Given: "+posteriorLogs.length);
		if(m>0 && posteriorLogs[0].length!=k) throw new IllegalArgumentException("Invalid columns of posterior logs. Expected: "+k+" Given: "+posteriorLogs[0].length);
		Double logProb = calculateForward(observations,forwardLogs);
		calculateBackward(observations,backwardLogs);
		for(int i=0;i<posteriorLogs.length;i++) {
			for(int j=0;j<posteriorLogs[0].length;j++) {
				Double f = forwardLogs[i][j];
				Double b = backwardLogs[i][j];
				Double e = getEmission(j, observations.get(i), i);
				Double fTimesE = LogMath.logProduct(f,e);
				posteriorLogs[i][j] = LogMath.logProduct(b, fTimesE);
				if(i==posteriorLogs.length-1) logProb=LogMath.logSum(logProb, fTimesE);
			}
		}
		return logProb;
	}
	protected Double [][] calculatePosteriorLogs (List<? extends Object> observations) {
		//Init arrays if not previously created to avoid reassignment of posteriorLogs cache attribute
		int m = observations.size();
		int k = getNumStates();
		initArrays(m,k,true);
		calculatePosteriorLogs(observations,posteriorLogs);
		return posteriorLogs;
	}
	

	@Override
	public void calculatePosteriors(List<? extends Object> observations, double[][] posteriors) {
		int m = observations.size();
		int k = getNumStates();
		if(posteriors.length!=m) throw new IllegalArgumentException("Invalid rows of posterior logs. Expected: "+m+" Given: "+posteriors.length);
		if(m>0 && posteriors[0].length!=k) throw new IllegalArgumentException("Invalid columns of posteriors. Expected: "+k+" Given: "+posteriors[0].length);
		Double [][] posteriorLogs = calculatePosteriorLogs(observations);
		for(int i=0;i<posteriorLogs.length;i++) {
			//System.out.println("i: "+i+" observation: "+observations.get(i)+" posteriorLogs: "+posteriorLogs[i][0]+" "+posteriorLogs[i][1]);
			LogMath.normalizeLogs(posteriorLogs[i]);
			for(int j=0;j<posteriorLogs[i].length;j++) {
				posteriors[i][j] = LogMath.power10(posteriorLogs[i][j]);
			}
		}
	}

	@Override
	public Double calculateForward(List<? extends Object> observations, Double [][] forwardLogs) {
		int m = observations.size();
		int n = getNumStates();
		if(forwardLogs.length!=m) throw new IllegalArgumentException("Invalid rows of forward logs. Expected: "+m+" Given: "+forwardLogs.length);
		if(m>0 && forwardLogs[0].length!=n) throw new IllegalArgumentException("Invalid columns of forwardLogs. Expected: "+n+" Given: "+forwardLogs.length);
		
		//Array to precalculate forward times emission
		Double [] fTimesE = new Double[n]; 
		for(int i=0;i<m;i++) {
			Object lastO = null;
			if(i>0) {
				lastO = observations.get(i-1);
				Arrays.fill(fTimesE, null);
				for(int k=0;k<n;k++) {
					Double e = getEmission(k, lastO, i-1);
					fTimesE[k] = LogMath.logProduct(forwardLogs[i-1][k],e);
				}
			}
			
			for(int j=0;j<n;j++) {
				if(lastO==null) forwardLogs[i][j] = getStart(j);
				else {
					//The sum of probabilities starts with zero which in logarithm is represented as null
					forwardLogs[i][j] = null;
					for(int k=0;k<n;k++) {
						Double t = getTransition(k, j, i-1);
						forwardLogs[i][j]=LogMath.logSum(forwardLogs[i][j],LogMath.logProduct(fTimesE[k],t));
					}
				}
			}	
		}
		//Calculate final probability
		return getSequenceLogProb(observations, forwardLogs);
	}
	
	/**
	 * Method that recalculates forward probabilities and returns the internal array where they get saved
	 * @param observations
	 * @return Double [][] array with forward probabilities with as many rows as observations and as many columns as states
	 */
	protected Double [][] calculateForward (List<? extends Object> observations) {
		int m = observations.size();
		int k = getNumStates();
		initArrays(m, k, true);
		calculateForward(observations,forwardLogs);
		return forwardLogs;
	}
	
	/**
	 * Calculates the total probability of the given observations
	 * @param observations list of m observations
	 * @param forwardLogs precalculated forward logs
	 * @return Double log of the total probability of the sequence of observations
	 */
	protected Double getSequenceLogProb(List<? extends Object> observations, Double [][] forwardLogs) {
		Double logProb = null;
		int m = observations.size();
		for(int j=0;j<forwardLogs[m-1].length;j++) {
			Double f = forwardLogs[m-1][j];
			Double e = getEmission(j, observations.get(m-1), m-1);
			logProb=LogMath.logSum(logProb, LogMath.logProduct(f,e));
		}
		return logProb;
	}

	@Override
	public void calculateBackward(List<? extends Object> observations, Double [][] backwardLogs) {
		int m = observations.size();
		int n = getNumStates();
		
		if(backwardLogs.length!=m) throw new IllegalArgumentException("Invalid rows of backwardLogs. Expected: "+m+" Given: "+backwardLogs.length);
		if(m>0 && backwardLogs[0].length!=n) throw new IllegalArgumentException("Invalid columns of backwardLogs. Expected: "+n+" Given: "+backwardLogs.length);
		
		Double [] bTimesE = new Double[n]; 
		for(int i=m-1;i>=0;i--) {
			Object lastO = null;
			if(i<m-1) {
				lastO = observations.get(i+1);
				for(int k=0;k<n;k++) {
					Double e = getEmission(k, lastO, i+1);
					bTimesE[k] = LogMath.logProduct(backwardLogs[i+1][k],e);
				}
			}
			for(int j=0;j<n;j++) {
				if(lastO==null) backwardLogs[i][j] = (double) 0;
				else {
					//The sum of probabilities starts with zero which in logarithm is represented as null
					backwardLogs[i][j]=null;
					for(int k=0;k<n;k++) {
						Double t = getTransition(j, k, i);
						backwardLogs[i][j]=LogMath.logSum(backwardLogs[i][j],LogMath.logProduct(bTimesE[k],t));
					}
				}
			}	
		}
	}
	
	/**
	 * Method that recalculates backward probabilities and returns the internal array where they get saved
	 * @param observations
	 * @return Double [][] array with backward probabilities with as many rows as observations and as many columns as states
	 */
	protected Double [][] calculateBackward (List<? extends Object> observations) {
		int m = observations.size();
		int k = getNumStates();
		initArrays(m, k, true);
		calculateBackward(observations, backwardLogs);
		return backwardLogs;
	}

	@Override
	public Double getViterbiPath(List<? extends Object> observations, int [] path) {
		int m = observations.size();
		int n = getNumStates();
		initArraysViterbi(m, n);
		//Array to precalculate viterbi times emission
		Double [] vTimesE = new Double[n];
		for(int i=0;i<m;i++) {
			Object lastO = null;
			if(i>0) {
				lastO = observations.get(i-1);
				Arrays.fill(vTimesE, null);
				for(int k=0;k<n;k++) {
					Double e = getEmission(k, lastO, i-1);
					vTimesE[k] = LogMath.logProduct(viterbiLogs[i-1][k],e);
				}
			}
			
			for(int j=0;j<n;j++) {
				if(lastO==null) {
					viterbiLogs[i][j] = getStart(j);
					viterbiBacktrace[i][j] = -1;
				}
				else {
					//The max probabilities starts with zero which in logarithm is represented as null
					viterbiLogs[i][j] = null;
					viterbiBacktrace[i][j] = -1;
					for(int k=0;k<n;k++) {
						Double t = getTransition(k, j, i-1);
						Double prob = LogMath.logProduct(vTimesE[k],t);
						if(prob != null && (viterbiLogs[i][j] == null || viterbiLogs[i][j] < prob) ) {
							viterbiLogs[i][j] = prob;
							viterbiBacktrace[i][j] = k;
						}
					}
				}
			}	
		}
		Double bestP = null;
		int bestState = -1;
		for(int j=0;j<n;j++) {
			Double p = viterbiLogs[m-1][j];
			p = LogMath.logProduct(p,getEmission(j, observations.get(m-1), m-1));
			if(p!=null && (bestP == null || bestP < p)) {
				bestState = j;
				bestP = p;
			}
		}
		if(bestP == null) {
			return bestP;
		}
		//Backtrace best path
		for(int i=m-1;i>=0;i--) {
			path[i] = bestState;
			bestState = viterbiBacktrace[i][bestState];
		}
		return bestP;
	}
	

	
	public static void calculateUniformChangeTransitions(double changeProbability, Double [][]transitions) {
		int n = transitions.length;
		
		double noChangeP = 1.0-changeProbability;
		Double logNoChange = LogMath.log10(noChangeP);
		//The probability of recombination is split uniformly across the parents
		Double logChange1 = LogMath.log10(changeProbability/(n-1));
		for(int j=0;j<n;j++) {
			for(int k=0;k<n;k++) {
				if(j==k) transitions[j][k] = logNoChange;
				else transitions[j][k] = logChange1;
			}
		}
	}
	private void initArrays(int m, int k, boolean initPosteriors) {
		if(forwardLogs.length!=m || forwardLogs[0].length!=k) {
			getLog().info("Creating array for forward probabilities of dimensions "+m+" x "+k);
			forwardLogs = new Double[m][k];
		}
		if(backwardLogs.length!=m || backwardLogs[0].length!=k) {
			getLog().info("Creating array for backward probabilities of dimensions "+m+" x "+k);
			backwardLogs = new Double[m][k];
		}
		if(initPosteriors && (posteriorLogs.length!=m || posteriorLogs[0].length!=k)) {
			getLog().info("Creating array for posterior probabilities of dimensions "+m+" x "+k);
			posteriorLogs = new Double[m][k];
		}
	}
	private void initArraysViterbi(int m, int k) {
		if(viterbiLogs.length!=m || viterbiLogs[0].length!=k) {
			getLog().info("Creating array for viterbi probabilities of dimensions "+m+" x "+k);
			viterbiLogs = new Double[m][k];
		}
		if(viterbiBacktrace.length!=m || viterbiBacktrace[0].length!=k) {
			getLog().info("Creating array for viterbi backtrack of dimensions "+m+" x "+k);
			viterbiBacktrace = new int[m][k];
		}
	}
	
}
