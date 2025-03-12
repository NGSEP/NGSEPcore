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

import ngsep.math.LogMath;

/**
 * @author Leonidas Villamil
 * @author Jorge Duitama
 */
public class ProfileAlignmentHMM extends AbstractHMM {
	private ProfileAlignmentNullModel nullmodel;
	private String id;
	private String name;
	private String domainCode;
	private List<? extends HMMState> states;
	private int numStates;
	private int steps;
	private Double [][][] transitionMatrix;
	private Double [][][] viterbiMatrix;
	private byte [][][] tracebackMatrix;
	private Double miu;
	private Double lambda;
	

	public ProfileAlignmentHMM(String id , int steps, List<? extends HMMState> states, ProfileAlignmentNullModel nullModel) {
		this.id=id;
		this.steps=steps;
		this.numStates=states.size();
		this.states=states;
		this.nullmodel = nullModel;
	}
	
	
	public String getId() {
		return id;
	}
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	
	public String getDomainCode() {
		return domainCode;
	}
	public void setDomainCode(String domainCode) {
		this.domainCode = domainCode;
	}


	public Double getMiu() {
		return miu;
	}

	public void setMiu(Double miu) {
		this.miu = miu;
	}

	public Double getLambda() {
		return lambda;
	}

	public void setLambda(Double lambda) {
		this.lambda = lambda;
	}

	public void setTransitionMatrix(Double[][][] transitionMatrix) {
		this.transitionMatrix = transitionMatrix;
	}

	public Double[][][] getTransitionMatrix() {
		return transitionMatrix;
	}
	
	public int getSteps() {
		return steps;
	}


	@Override
	public int getNumStates() {
		return numStates;
	}


	@Override
	public HMMState getState(int index) {
		return states.get(index);

	}

	@Override
	public Double getTransition(int source, int dest, int step) {
		// Answer with match probability
		return transitionMatrix[step][source][dest];
	}

	public ProfileAlignmentDomain findDomain(String query) {
		ProfileAlignmentDomain domain = null;
		calculateViterbiPathMatrix(query);
		StringBuilder profile = new StringBuilder();
		StringBuilder alignmentQuery = new StringBuilder();
		char[] queryArray = query.toCharArray();
		
		// Start in the last position
		int k = 0; // end in match
		int j = steps-1;
		int i = getMaxProbabilityQuery();
		if(i<0) {
			return domain;
			//i = query.length();
			//j = getMaxProbabilityHMM();
			//if(j<0) return domain;
		}
		int queryEnd = i-1;
		
		// Look for max value
		Double maxScore =viterbiMatrix[i][j][k];
		//System.out.println("Max raw score: "+maxScore);

		while (i > 0 && j >0) {
			// Depending on the state k, it could be "match", "insertion" o "deletion"
			if (k == 0) {  // Match
				profile.append('X');
				alignmentQuery.append(queryArray[i-1]);
				i--;
				j--;
			} else if (k == 1) {  // Insertion
				profile.append('I');
				alignmentQuery.append(queryArray[i-1]);
				i--;
			} else if (k == 2) {  // Deletion
				profile.append('X');
				alignmentQuery.append('-');
				j--;
			}
			byte prevState = tracebackMatrix[i][j][k];
			if (prevState != -1) {
				k = prevState;  // Previous state
			} else {
				break;
			}
		}
		//System.out.println("Start sequence: "+i);
		int queryStart=i;
		int m = queryEnd-queryStart;
		Double eValue = calculateEValue(maxScore, m);
		//System.out.println("Limits query in profile: " +queryStart+ "\t"+ queryEnd+" evalue: "+eValue);
		if (passFilter(eValue)) {
			domain = new ProfileAlignmentDomain(query, queryStart, queryEnd, id);
			domain.setDomainCode(domainCode);
			domain.setEvalue(eValue);
			for (;j>0;j--) profile.append('S');
			domain.setAlignment(alignmentQuery.reverse().toString(), profile.reverse().toString());
		}
		return domain;
	}
	
	private void calculateViterbiPathMatrix(String query){
	    int querySize = query.length();
	    char[] queryArray = query.toCharArray();
	    viterbiMatrix = new Double[querySize+1][steps][3];
	    tracebackMatrix = new byte[querySize+1][steps][3];
	    // i represtenta el tamaño del prefijo del query 
	    for (int i = 0; i < querySize+1; i++) {
	        for (int j = 0; j < steps; j++) {
	            for (int k = 0; k < 3; k++) {
	                byte prevState = -1; // Usar -1 para indicar valor no asignado
	                // base cases
	                if (i==0) {
	                	if(j==0 && (k==0)) {
	                		viterbiMatrix[i][j][k] = 0.0;
		                    tracebackMatrix[i][j][k] = 0; // Inicializar como -1 en lugar de null
	                	}
	                }
	                else if (i>0 && j==0 && k==0) {
	                	viterbiMatrix[i][j][k] = 0.0;
	                    tracebackMatrix[i][j][k] = 0;
	                }
	                else {

	                	// cases for the first states after the start viterbi of the start is 1 -> 0 (log)
	                	if (j==1 && k==0) {
	                		Double emissionProbability=states.get(k).getEmission(queryArray[i-1], j);
	                		emissionProbability = emissionProbability-nullmodel.calculateScore(query.substring(i-1, i));
	                		Double fromOrigin=transitionMatrix[0][0][0];
	                		Double viterbiInsertion=viterbiMatrix[i-1][j-1][1];
		                	Double transition= transitionMatrix[0][1][0];
	                		Double fromInsertion=LogMath.logProduct(transition , viterbiInsertion);
	                		Double max = LogMath.logMax(fromOrigin, fromInsertion);
	                		viterbiMatrix[i][j][k]=LogMath.logProduct(emissionProbability, max);
	                		if (max != null) {
		                        prevState = (byte) (max.equals(fromOrigin) ? 0 : 1);
		                        tracebackMatrix[i][j][k] = prevState;
		                    } else {
		                        tracebackMatrix[i][j][k] = -1;
		                    }
	                	}
	                	else if (j==0 && k==1) {
	                		Double emissionProbability=states.get(k).getEmission(queryArray[i-1], j);
	                		emissionProbability = emissionProbability-nullmodel.calculateScore(query.substring(i-1, i));
	                		Double fromOrigin=transitionMatrix[0][0][1];
	                		Double transition=transitionMatrix[0][1][1];
		                	Double viterbi=viterbiMatrix[i-1][j][1];
	                		Double fromInsertion=LogMath.logProduct( transition, viterbi);
	                		Double max = LogMath.logMax(fromOrigin, fromInsertion);
	                		viterbiMatrix[i][j][k]=LogMath.logProduct(emissionProbability, max);
	                		if (max != null) {
		                        prevState = (byte) (max.equals(fromOrigin) ? 0 : 1);
		                        tracebackMatrix[i][j][k] = prevState;
		                    } else {
		                        tracebackMatrix[i][j][k] = -1;
		                    }
	                	}
	                	else if(j==1 && k==2) {
	                		Double fromOrigin=transitionMatrix[0][0][2];
	                		Double viterbiInsertion=viterbiMatrix[i-1][j-1][1];
		                	Double transition= transitionMatrix[0][1][2];
	                		Double fromInsertion=LogMath.logProduct(transition , viterbiInsertion);
	                		Double max = LogMath.logMax(fromOrigin, fromInsertion);
	                		viterbiMatrix[i][j][k]=max;
	                		if (max != null) {
		                        prevState = (byte) (max.equals(fromOrigin) ? 0 : 1);
		                        tracebackMatrix[i][j][k] = prevState;
		                    } else {
		                        tracebackMatrix[i][j][k] = -1;
		                    }
	                	}
	                // match case (k = 0)
	                else if (j > 0 && k == 0) {
	                    Double emissionProbability = states.get(k).getEmission(queryArray[i - 1], j);
	                    emissionProbability = emissionProbability-nullmodel.calculateScore(query.substring(i-1, i));
	                    Double fromMatch = LogMath.logProduct(viterbiMatrix[i - 1][j - 1][0], transitionMatrix[j - 1][0][0]);
	                    Double fromInsertion = LogMath.logProduct(viterbiMatrix[i - 1][j - 1][1], transitionMatrix[j - 1][1][0]);
	                    Double fromDeletion = LogMath.logProduct(viterbiMatrix[i - 1][j - 1][2], transitionMatrix[j - 1][2][0]);
	                    
	                    Double max = LogMath.logMax(fromMatch, LogMath.logMax(fromInsertion, fromDeletion));
	                    viterbiMatrix[i][j][k] = LogMath.logProduct(emissionProbability, max);

	                    // Traceback
	                    if (max != null) {
	                        if (max.equals(fromMatch)) prevState = 0;
	                        else if (max.equals(fromInsertion)) prevState = 1;
	                        else prevState = 2;
	                        
	                        tracebackMatrix[i][j][k] = prevState;
	                    } else {
	                        tracebackMatrix[i][j][k] = -1;
	                    }
	                }
	                
	                // insertion case (k = 1)
	                else if (j > 0 && k == 1) {
	                    Double emissionProbability = states.get(k).getEmission(queryArray[i - 1], j);
	                    emissionProbability = emissionProbability-nullmodel.calculateScore(query.substring(i-1, i));
	                    Double fromMatch = LogMath.logProduct(viterbiMatrix[i - 1][j][0], transitionMatrix[j][0][1]);
	                    Double fromInsertion = LogMath.logProduct(viterbiMatrix[i - 1][j][1], transitionMatrix[j][1][1]);

	                    Double max = LogMath.logMax(fromInsertion, fromMatch);
	                    viterbiMatrix[i][j][k] = LogMath.logProduct(emissionProbability, max);

	                    if (max != null) {
	                        prevState = (byte) (max.equals(fromMatch) ? 0 : 1);
	                        tracebackMatrix[i][j][k] = prevState;
	                    } else {
	                        tracebackMatrix[i][j][k] = -1;
	                    }
	                }
	                
	                // deletion case (k = 2)
	                else if (j > 0 && k == 2) {
	                	Double penalty = 0.0; 
	                    Double fromMatch = LogMath.logProduct(viterbiMatrix[i][j - 1][0], transitionMatrix[j - 1][0][2]);
	                    Double fromDeletion = LogMath.logProduct(viterbiMatrix[i][j - 1][2], transitionMatrix[j - 1][2][2]);
	                    
	                    Double max = LogMath.logMax(fromDeletion, fromMatch);
	                    viterbiMatrix[i][j][k] = LogMath.logProduct(penalty, max);

	                    // Traceback
	                    if (max != null) {
	                        prevState = (byte) (max.equals(fromMatch) ? 0 : 2);
	                        tracebackMatrix[i][j][k] = prevState;
	                    } else {
	                        tracebackMatrix[i][j][k] = -1;
	                    }
	                    
	                }
	            }
	            }
	        }
	    }
	    
	    //Double fromMatch = LogMath.logProduct(viterbiMatrix[querySize-1][steps-1][0],transitionMatrix[steps-1][0][0]);
	    //Double fromInsertion = LogMath.logProduct(viterbiMatrix[querySize-1][steps-1][1],transitionMatrix[steps-1][1][0]);
	    //Double fromDeletion = LogMath.logProduct(viterbiMatrix[querySize-1][steps-1][2],transitionMatrix[steps-1][2][0]);
	    //finalNodeProbability=LogMath.logMax(fromDeletion, LogMath.logMax(fromInsertion, fromMatch));
	    /*if (finalNodeProbability != null) {
	    	finalNodeTraceback = (byte) (finalNodeProbability.equals(fromMatch) ? 0 : 1);
        } else {
        	finalNodeTraceback = -1;
        }*/
//	    printViterbiMatrix(0);
//	    printViterbiMatrix(2);
	}
	
	private int getMaxProbabilityQuery() {
		Double maxProbability= -5.0*steps;
		int position=-1;
		for (int i=0; i<viterbiMatrix.length;i++) {
			Double probability=viterbiMatrix[i][steps-1][0];
			//if("PF00077.25".equals(id) && probability!=null && probability>-5) System.out.println("probability "+i+": "+probability);
			if (probability!=null && probability>maxProbability) {
				maxProbability=probability;
				position=i;
			}
		}
		//if("PF00077.25".equals(id)) System.out.println("Maximum probability: "+maxProbability+" pos: "+position);
		return position;
	}
	private int getMaxProbabilityHMM() {
		Double maxProbability=-5.0*steps;
		int n = viterbiMatrix.length-1;
		int position=-1;
		for (int j=10; j<viterbiMatrix[0].length;j++) {
			Double probability=viterbiMatrix[n][j][0];
			//if("PF00077.25".equals(id)) System.out.println("probability "+j+": "+probability);
			if (probability!=null && probability>maxProbability) {
				maxProbability=probability;
				position=j;
				//System.out.println("Local maximum "+maxProbability+" pos: "+i);
			}
		}
		return position;
	}
	
	public void printViterbiMatrix(int k) {
		for (int i = 0; i < viterbiMatrix.length; i++) {
			for(int j=0 ;j<viterbiMatrix[i].length; j++) {
				Double log = viterbiMatrix[i][j][k];
				System.out.print(log+" ");
			}
			System.out.println();
		}
	}
	
	public void printTracebackMatrix(int k) {
		for (int i = 0; i < tracebackMatrix.length; i++) {
			for(int j=0 ;j<tracebackMatrix[i].length; j++) {
				byte value = tracebackMatrix[i][j][k];			
				System.out.print(value+" ");
			}
			System.out.println();
		}
	}
	
	public Double calculateEValue (Double maxScore,int m) {
		//Double density = 1- Math.exp(-Math.exp(-lambda*(maxScore-miu)));
		//Double eValue = 1- Math.exp(-density)*nModel;
		//System.out.println("New E value: "+ eValue);
		Double asd = (double) (m*(this.steps-1)*Math.pow(10,-lambda*(maxScore-miu)));
		//System.out.println("Old eValue: " + asd);
		return asd;
	}
	
	public boolean passFilter (Double evalue) {
		if (evalue<0.01) {return true;}
		else {return false;}
	}

}
