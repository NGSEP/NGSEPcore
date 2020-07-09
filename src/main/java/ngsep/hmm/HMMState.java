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

public interface HMMState {
	/**
	 * Returns the logarithm (base 10) of the probability of emission of the given value
	 * @param value that will be emitted
	 * @param step At which the value is emitted
	 * @return Double log10 of the probability of observing the given value
	 * Null if the probability is zero 
	 */
	public Double getEmission(Object value, int step);
	
	/**
	 * Returns the logarithm (base 10) of the probability of starting at this state
	 * @return double log10 of the probability of starting at this state
	 * Null if the probability is zero 
	 */
	public Double getLogStart();
	/**
	 * Changes the probability of starting at this state
	 * @param logStart log10 of the new probability of starting at this state
	 * Null if the probability is zero 
	 */
	public void setLogStart(Double logStart);
	/**
	 * Returns the id of the state
	 * @return String id assigned to the state
	 */
	public String getId();
}
