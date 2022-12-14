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
package ngsep.genome;

import java.util.List;

/**
 * 
 * @author Jorge Duitama
 *
 */
public interface PairwiseSyntenyBlocksFinder {
	public static final int DEF_MIN_BLOCK_LENGTH = 100000;
	public static final int DEF_MIN_HOMOLOGY_UNITS_BLOCK = 6;
	public static final int DEF_MAX_DISTANCE_BETWEEN_UNITS = 100000;
	/**
	 * Calculates synteny blocks based on orthogroups defined by homology clusters
	 * @param g1 First genome to compare
	 * @param g2 Second genome to compare
	 * @param clusters Homology clusters
	 * @return List<PairwiseSyntenyBlock> List of pairwise synteny homolog clusters between the two groups 
	 */
	public List<PairwiseSyntenyBlock> findSyntenyBlocks(AnnotatedReferenceGenome g1, AnnotatedReferenceGenome g2, List<HomologyCluster> clusters);
}
