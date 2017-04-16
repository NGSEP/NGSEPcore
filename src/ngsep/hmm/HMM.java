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

public interface HMM {

	/**
	 * Returns the logarithm (base 10) of the probability of transition between the source and
	 * the dest states at the given step
	 * @param source First state
	 * @param dest Second state
	 * @param step Step at which the transition will happen
	 * @return Double log10 of the transition probability between source and dest at step.
	 * Null if the probability is zero 
	 */
	public Double getTransition(int source, int dest, int step );
	/**
	 * Returns the logarithm (base 10) of the emission probability of the given value by the given state
	 * at the given step 
	 * @param state From which the value is emitted
	 * @param value observed value
	 * @param step At which the value is emitted
	 * @return Double log10 of the emission probability of the given value by the given state
	 * at the given step
	 * Null if the probability is zero
	 */
	public Double getEmission(int state, Object value, int step);
	/**
	 * Returns the logarithm (base 10) of the initial probability of the given state
	 * @param state Potential initial state
	 * @return Double log10 of the probability of starting at the given state
	 * Null if the probability is zero
	 */
	public Double getStart(int state);
	/**
	 * Returns the state at the given position
	 * @param state Position of the state in the HMM
	 * @return HMMState Object representing the state
	 */
	public HMMState getState(int state);
	/**
	 * Returns the number of states
	 * @return int number of states
	 */
	public int getNumStates();
	/**
	 * Run the forward-backward algorithm to calculate posterior probabilities of states given a set of observations
	 * @param observations List of observed values
	 * @param posteriorLogs Output matrix with as many rows as observations and as many columns as states. It is
	 * designed as a parameter instead of a return value to avoid constant reallocation and to allow returning the
	 * probability of the data as a return value. Zero probabilities are represented as null objects
	 * @return Double log10 of the probability of the data given the HMM
	 * Null if the probability is zero   
	 */
	public Double calculatePosteriorLogs (List<? extends Object> observations, Double [][] posteriorLogs);
	
	/**
	 * Run the forward-backward algorithm to calculate posterior probabilities of states given a set of observations
	 * @param observations List of observed values
	 * @param posteriors Output matrix with as many rows as observations and as many columns as states. It is
	 * designed as a parameter instead of a return value to avoid constant reallocation and to allow returning the
	 * probability of the data as a return value. Unlike the method calculatePosteriorLogs, these are actual probabilities
	 * normalized row by row
	 */
	public void calculatePosteriors (List<? extends Object> observations, double [][] posteriors);
	
	/**
	 * Calculate forward log probabilities for each state at each step
	 * @param observations List of observed values
	 * @param forwardLogs Output matrix with as many rows as observations and as many columns as states.
	 * It is designed as a parameter instead of a return value to avoid constant reallocation
	 * To facilitate calculations, forward probabilities do not include the emission probability at
	 * each state.
	 * @return Double log10 of the probability of the data given the HMM
	 */
	public Double calculateForward(List<? extends Object> observations, Double [][] forwardLogs);
	/**
	 * Calculate backward log probabilities for each state at each step
	 * @param observations List of observed values
	 * @param backwardLogs Output matrix with as many rows as observations and as many columns as states.
	 * It is designed as a parameter instead of a return value to avoid constant reallocation
	 */
	public void calculateBackward(List<? extends Object> observations, Double [][] backwardLogs);
	/**
	 * 
	 * @param observations to calculate the path with the best probability
	 * @param path Output path
	 * @return Double Logarithm of th probability of the best path. Null if all paths have zero probability
	 */
	public Double getViterbiPath (List<? extends Object> observations, int [] path );
	
	
	
}
