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
package ngsep.transposons;

import java.util.List;

import ngsep.genome.ReferenceGenome;

public interface DeNovoTransposableElementsFinder {

	/**
	 * Finds transposable elements in a refrence genome
	 * @param genome to analyze
	 * @return List<TransposableElement> List of sequences identified as transposable elements
	 */
	public List<TransposableElementAnnotation> findTransposons(ReferenceGenome genome);
}
