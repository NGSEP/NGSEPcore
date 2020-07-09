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

public class ConstantTransitionHMM extends AbstractHMM {
	private List<? extends HMMState> states;
	private int n;
	private Double [][] transitions;
	
	/**
	 * @param states
	 */
	public ConstantTransitionHMM(List<? extends HMMState> states) {
		super();
		this.states = states;
		n = states.size();
		transitions = new Double[n][n];
	}

	public void setTransitions(Double[][] transitions) {
		if(transitions.length!=n) throw new IllegalArgumentException("Transitions matrix should have the same number of rows as states of the HMM. States: "+n+" rows: "+transitions.length);
		for(int i=0;i<n;i++) {
			if(transitions[i].length!=n) throw new IllegalArgumentException("Transitions matrix should have the same number of columns as states of the HMM. States: "+n+" columns: "+transitions[i].length);
			for(int j=0;j<n;j++) this.transitions[i][j] = transitions[i][j];
		}
	}


	@Override
	public Double getTransition(int source, int dest, int step) {
		return transitions[source][dest];
	}
	
	@Override
	public HMMState getState(int state) {
		return states.get(state);
	}

	@Override
	public int getNumStates() {
		return n;
	}
	public void calculateUniformChangeTransitions(double changeProbability) {
		AbstractHMM.calculateUniformChangeTransitions(changeProbability, this.transitions);
	}
	

}
