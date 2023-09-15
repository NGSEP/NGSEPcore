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
package ngsep.alignments;

/**
 * 
 * @author Jorge Duitama
 *
 */
public interface UngappedSearchHitsClusterAligner {
	public static final int ALIGNMENT_ALGORITHM_AFFINE_GAP = 1;
	public static final int ALIGNMENT_ALGORITHM_DYNAMIC_KMERS = 2;
	public static final int ALIGNMENT_ALGORITHM_SIMPLE_GAP = 3;
	public static final int ALIGNMENT_ALGORITHM_NAIVE = 4;
	public static final int ALIGNMENT_ALGORITHM_SHORT_READS = 5;
	public static final int ALIGNMENT_ALGORITHM_STATIC_BAND = 6;
	/**
	 * Performs the alignment process to build a read alignment from a cluster of ungapped hits
	 * @param query sequence to be aligned
	 * @param subject sequence for the alignment
	 * @param cluster of shared substrings between the query and the subject  
	 * @return ReadAlignment object describing the alignment between the query and the subject
	 */
	public ReadAlignment buildAlignment(CharSequence query, CharSequence subject, UngappedSearchHitsCluster hitsCluster);
}
