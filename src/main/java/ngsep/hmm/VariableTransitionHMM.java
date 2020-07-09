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
import java.util.Random;

import ngsep.math.LogMath;

public class VariableTransitionHMM extends AbstractHMM {

	private int iterationsBaumWelch = DEF_ITER_BAUM_WELCH;
	
	private List<? extends HMMState> states;
	private int steps;
	private int numStates;
	private Double [][][] logTransitions;
	private boolean skipTransitionsTraining = false;
	//Local arrays to save reallocation over many runs
	private Double [] logStarts = new Double [0];
	private Double [][][] logTransitionsTrain = new Double [0][0][0];
	private List<List<? extends Object>> trainingData = null;
	
	/**
	 * Creates a new VariableTransitionHMM with the given states and the given steps
	 * @param states States of the HMM
	 * @param steps Number of steps at which the HMM will run with variable transitions among the states
	 */
	public VariableTransitionHMM(List<? extends HMMState> states, int steps) {
		this.states = states;
		numStates = states.size();
		this.steps = steps;
		getLog().info("Creating array for transitions of dimensions "+(steps-1)+" x "+numStates+" x "+numStates);
		logTransitions = new Double [steps-1][numStates][numStates];
	}
	
	public int getIterationsBaumWelch() {
		return iterationsBaumWelch;
	}


	public void setIterationsBaumWelch(int iterationsBaumWelch) {
		this.iterationsBaumWelch = iterationsBaumWelch;
	}


	public boolean isSkipTransitionsTraining() {
		return skipTransitionsTraining;
	}

	public void setSkipTransitionsTraining(boolean skipTransitionsTraining) {
		this.skipTransitionsTraining = skipTransitionsTraining;
	}
	
	public List<List<? extends Object>> getTrainingData() {
		return trainingData;
	}

	public void setTrainingData(List<List<? extends Object>> trainingData) {
		this.trainingData = trainingData;
	}

	public void setTransitions(Double[][] logTransitions, int step) {
		if(logTransitions.length!=numStates) throw new IllegalArgumentException("Transitions matrix should have the same number of rows as states of the HMM. States: "+numStates+" rows: "+logTransitions.length);
		for(int i=0;i<numStates;i++) {
			if(logTransitions[i].length!=numStates) throw new IllegalArgumentException("Transitions matrix should have the same number of columns as states of the HMM. States: "+numStates+" columns: "+logTransitions[i].length);
			//Normalize before update
			LogMath.normalizeLogs(logTransitions[i]);
			for(int j=0;j<numStates;j++)  this.logTransitions[step][i][j] = logTransitions[i][j];
		}
	}
	
	public void setRandomTransitions() {
		Double[][] logRandom = new Double[numStates][numStates];
		Random r = new Random();
		//System.out.println("Random transitions for "+steps+" steps");
		for(int i=0;i<steps-1;i++) {
			for(int j=0;j<numStates;j++) {
				for(int k=0;k<numStates;k++) {
					//TODO: Improve sampling from Direlecht
					logRandom[j][k] = Math.log10(r.nextDouble()*0.6+0.2);
				}
			}
			setTransitions(logRandom, i);
		}
		//printTransitions(0);
	}
	
	public void calculateUniformChangeTransitions(double[] changeProbabilities) {
		int n = this.getNumStates();
		int m = getSteps();
		if(m-1!=changeProbabilities.length) throw new IllegalArgumentException("Length of changes vector "+changeProbabilities.length+" is not consistent with the number of steps "+m+". It should be "+(m-1));
		if(changeProbabilities.length>0)getLog().info("Using change probabilities to infer transitions. First probability: "+changeProbabilities[0]);
		Double [][] transitions = new Double [n][n];
		for(int i=0;i<m-1;i++) {
			double p = changeProbabilities[i];
			AbstractHMM.calculateUniformChangeTransitions(p, transitions);
			setTransitions(transitions, i);
			//if (i==0) System.out.println("Transition between "+states.get(0).getId()+" and "+states.get(1).getId()+": "+this.transitions[0][1][0]);
		}	
	}

	@Override
	public Double getTransition(int source, int dest, int step) {
		return logTransitions[step][source][dest];
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
	
	public void train() {
		if (trainingData==null) throw new IllegalStateException("Training data must be provided");
		getLog().info("Training model with "+trainingData.size()+" sequences");
		int n = getNumStates();
		setRandomTransitions();
		//printTransitions(0);
		//printTransitions(2);
		double logUniformStart = Math.log10(1.0/(double)n);
		for(int j=0;j<n;j++) {
			HMMState state = getState(j);
			state.setLogStart(logUniformStart);
			randomizeEmissions(j);
		}
		for(int h = 0; h < iterationsBaumWelch; h++) {
			getLog().info("Running "+h+" Baum-Welch iteration");
			runBaumWelchStep();
		}
		//printTransitions(0);
		//printTransitions(2000);
	}
	
	
	/**
	 * Runs a step of baum-welch training with the attribute training data
	 */
	protected void runBaumWelchStep() {
		initArrays();
		Arrays.fill(logStarts, null);
		for(int i=0;i<logTransitionsTrain.length;i++) {
			for(int j=0;j<logTransitionsTrain[i].length;j++) {
				Arrays.fill(logTransitionsTrain[i][j], null);
			}
		}
		initEmissionsBaumWelch();
		//int datumIdx = 0;
		for (List<? extends Object> trainingDatum:trainingData) {
			Double [][] forwardLogs=calculateForward(trainingDatum);
			Double [][] backwardLogs=calculateBackward(trainingDatum);
			Double logProb = getSequenceLogProb(trainingDatum, forwardLogs);
			//Calculate new starts
			for(int j=0;j<logStarts.length;j++) {
				Object o = trainingDatum.get(0);
				Double seqProduct = LogMath.logProduct(forwardLogs[0][j], backwardLogs[0][j]);
				Double emission = getEmission(j, o, 0);
				seqProduct = LogMath.logProduct(seqProduct, emission);
				seqProduct = LogMath.logProduct(seqProduct, -logProb);
				//if(seqProduct>-0.5) System.out.println("Datum: "+datumIdx+". Next most likely start: "+j+" forward: "+forwardLogs[0][j]+" backward: "+backwardLogs[0][j]+" emission: "+emission+" logProb: "+logProb+"seq product: "+seqProduct);
				logStarts[j] = LogMath.logSum(logStarts[j], seqProduct);
			}
			//Calculate new transitions
			if(!skipTransitionsTraining) {
				for(int i=0;i<logTransitionsTrain.length;i++) {
					for(int j=0;j<logTransitionsTrain[i].length;j++) {
						
						for(int k=0;k<logTransitionsTrain[i][j].length;k++) {
							Object o1 = trainingDatum.get(i);
							Object o2 = trainingDatum.get(i+1);
							
							Double seqProduct = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i+1][k]);
							seqProduct = LogMath.logProduct(seqProduct, getEmission(j, o1, i));
							seqProduct = LogMath.logProduct(seqProduct, getEmission(k, o2, i+1));
							seqProduct = LogMath.logProduct(seqProduct, getTransition(j, k, i));
							seqProduct = LogMath.logProduct(seqProduct, -logProb);
							logTransitionsTrain[i][j][k] = LogMath.logSum(logTransitionsTrain[i][j][k], seqProduct);
						}
					}
				}
			}
			
			//Calculate new emissions
			for(int i=0;i<steps;i++) {
				Object o = trainingDatum.get(i);
				for(int j=0;j<numStates;j++) {
					Double seqProduct = LogMath.logProduct(forwardLogs[i][j], backwardLogs[i][j]);
					seqProduct = LogMath.logProduct(seqProduct, getEmission(j, o, i));
					seqProduct = LogMath.logProduct(seqProduct, -logProb);
					accumulateEmissionBaumWelch(i,j,o,seqProduct);
				}
			}
			//datumIdx++;
		}
		//Normalize and update starts
		Double total = null;
		for(int j=0;j<logStarts.length;j++) total = LogMath.logSum(total, logStarts[j]);
		for(int j=0;j<logStarts.length;j++) getState(j).setLogStart(LogMath.logProduct(logStarts[j],-total));
		//Normalize and update transitions
		if(!skipTransitionsTraining) {
			for(int i=0;i<logTransitionsTrain.length;i++) {
				setTransitions(logTransitionsTrain[i], i);
			}
		}
		
		//Normalize and update emissions
		for(int j=0;j<numStates;j++) {
			updateEmissionsBaumWelch(j);
		}
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
	/**
	 * Randomize emission parameters
	 * WARN: This method will throw a RuntimeException because it should be redefined to use automated training
	 * @param stateIndex Index of the state to randomize
	 */
	protected void randomizeEmissions(int stateIndex) {
		throw new RuntimeException("This method should be redefined in a subclass to use automated Baum-Welch training");
	}

	/**
	 * Initializes training emissions for a baum welch step
	 * WARN: This method will throw a RuntimeException because it should be redefined to use automated training
	 */
	protected void initEmissionsBaumWelch() {
		throw new RuntimeException("This method should be redefined in a subclass to use automated Baum-Welch training");	
	}

	/**
	 * Registers a posterior probability of a fixed emission during Baum-Welch training
	 * WARN: This method will throw a RuntimeException because it should be redefined to use automated training
	 * @param step at which the observation is registered
	 * @param stateIndex Index of the state where the posterior was calculated
	 * @param datum observed value
	 * @param logPosterior Logarithm of the posterior probability
	 */
	protected void accumulateEmissionBaumWelch(int step, int stateIndex, Object datum, Double logPosterior) {
		throw new RuntimeException("This method should be redefined in a subclass to use automated Baum-Welch training");
	}

	/**
	 * Updates the emission probabilities during Baum-Welch training
	 * WARN: This method will throw a RuntimeException because it should be redefined in a subclass to use automated training
	 * @param stateIndex Index of the state to be updated
	 */
	protected void updateEmissionsBaumWelch(int stateIndex) {
		throw new RuntimeException("This method should be redefined in a subclass to use automated Baum-Welch training");
	}

	private void initArrays() {
		if(logTransitionsTrain.length!=steps-1 || logTransitionsTrain[0].length!=numStates) {
			getLog().info("Creating array for transitions of dimensions "+(steps-1)+" x "+numStates+" x "+numStates);
			logTransitionsTrain = new Double [steps-1][numStates][numStates];
		}
		if(logStarts.length!=numStates) logStarts = new Double [numStates];
	}
	
}
