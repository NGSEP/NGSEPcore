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

import java.util.List;
import java.util.Random;

import ngsep.math.LogMath;

public class VariableTransitionHMM extends AbstractHMM {

	private List<? extends HMMState> states;
	private int steps;
	private int numStates;
	private Double [][][] transitions;
	
	
	/**
	 * Creates a new VariableTransitionHMM with the given states and the given steps
	 * @param states States of the HMM
	 * @param steps Number of steps at which the HMM will run with variable transitions among the states
	 */
	public VariableTransitionHMM(List<? extends HMMState> states, int steps) {
		this.states = states;
		numStates = states.size();
		this.steps = steps;
		transitions = new Double [numStates][numStates][steps-1];
	}
	
	public void setTransitions(Double[][] logTransitions, int step) {
		if(logTransitions.length!=numStates) throw new IllegalArgumentException("Transitions matrix should have the same number of rows as states of the HMM. States: "+numStates+" rows: "+logTransitions.length);
		for(int i=0;i<numStates;i++) {
			if(logTransitions[i].length!=numStates) throw new IllegalArgumentException("Transitions matrix should have the same number of columns as states of the HMM. States: "+numStates+" columns: "+logTransitions[i].length);
			//Normalize before update
			LogMath.normalizeLogs(logTransitions[i]);
			for(int j=0;j<numStates;j++)  this.transitions[i][j][step] = logTransitions[i][j];
		}
	}
	
	public void setRandomTransitions() {
		Double[][] logRandom = new Double[numStates][numStates];
		Random r = new Random();
		for(int i=0;i<steps-1;i++) {
			for(int j=0;j<numStates;j++) {
				for(int k=0;k<numStates;k++) {
					//TODO: Improve sampling from Direlecht
					logRandom[j][k] = Math.log10(r.nextDouble()*0.6+0.2);
				}
			}
			setTransitions(logRandom, i);
		}
	}
	
	public void calculateUniformChangeTransitions(double[] changeProbabilities) {
		int n = this.getNumStates();
		int m = getSteps();
		if(m-1!=changeProbabilities.length) throw new IllegalArgumentException("Length of changes vector "+changeProbabilities.length+" is not consistent with the number of steps "+m+". It should be "+(m-1));
		Double [][] transitions = new Double [n][n];
		for(int i=0;i<m-1;i++) {
			double p = changeProbabilities[i];
			AbstractHMM.calculateUniformChangeTransitions(p, transitions);
			setTransitions(transitions, i);
		}	
	}

	@Override
	public Double getTransition(int source, int dest, int step) {
		return transitions[source][dest][step];
	}
	
	@Override
	public HMMState getState(int state) {
		return states.get(state);
	}

	@Override
	public int getNumStates() {
		return numStates;
	}

	public int getSteps() {
		return steps;
	}
	
}
