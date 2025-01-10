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

/**
 * @author Leonidas Villamil
 * @author Jorge Duitama
 */
public class ProfileAlignmentHMMState implements HMMState {
	private String id;
	private String alphabet;
	private int steps;
	private Double [][] emissionMatrix;
	
	public ProfileAlignmentHMMState (String id,String alphabet, int steps) {
		this.id=id;
		this.alphabet = alphabet;
		this.steps = steps;
		this.emissionMatrix=new Double[steps][alphabet.length()];
	}
	
	public String getAlphabet() {
		return alphabet;
	}

	public int getSteps() {
		return steps;
	}

	/**
	 * Returns the logarithm (base 10) of the probability of emission of the given value
	 * @param value that will be emitted
	 * @param step At which the value is emitted
	 * @return Double log10 of the probability of observing the given value
	 * Null if the probability is zero 
	 */
	@Override
    public Double getEmission(Object value, int step) {
        if (value instanceof Character) {
            char queryChar = (Character) value;
            int charIdx = alphabet.indexOf(queryChar);
            if(charIdx<0) throw new IllegalArgumentException("Unsupported character: "+queryChar+". Alphabet "+alphabet);
            return emissionMatrix[step][charIdx];
        } else {
            throw new IllegalArgumentException("Unsupported object: "+value+". class "+value.getClass().getName());
        }
    }

	@Override
	public String getId() {
		return id;
	}

	@Override
	public Double getLogStart() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setLogStart(Double arg0) {
		// TODO Auto-generated method stub
		
	}

	public Double[][] getEmissionMatrix() {
		return emissionMatrix;
	}

	public void setEmissionMatrix(Double[][] emissionMatrix) {
		this.emissionMatrix = emissionMatrix;
	}
	public void setStepEmissions(int step,Double [] stateEmissions) {
		emissionMatrix[step]=stateEmissions;
	}
}
