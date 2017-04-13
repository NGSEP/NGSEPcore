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
	
	private Logger log = Logger.getLogger(AbstractHMM.class.getName());
	private Double [][] forwardLogs=new Double[0][0];
	private Double [][] backwardLogs=new Double[0][0];
	
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
	public Double calculatePosteriors(List<? extends Object> observations,Double[][] posteriorLogs) {
		int m = observations.size();
		int n = getNumStates();
		if(forwardLogs.length!=m || forwardLogs[0].length!=n) {
			log.info("Creating matrix for forward probabilities of dimensions "+m+" x "+n);
			forwardLogs = new Double[m][n];
		}
		if(backwardLogs.length!=m || backwardLogs[0].length!=n) {
			log.info("Creating matrix for backward probabilities of dimensions "+m+" x "+n);
			backwardLogs = new Double[m][n];
		}
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

	@Override
	public Double calculateForward(List<? extends Object> observations, Double [][] forwardLogs) {
		int m = observations.size();
		int n = getNumStates();
		
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
		Double logProb = null;
		for(int j=0;j<n;j++) {
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

	@Override
	public double getViterbiPath(List<Object> observations, int[] path) {
		// TODO Auto-generated method stub
		return 0;
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
	
}
