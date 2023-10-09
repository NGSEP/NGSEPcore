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
package ngsep.clustering;
/**
 * @author Sebastian Lemus
 * @author Jorge Duitama
 *
 */
public interface DistanceMatrixClustering {
	/**
	 * Builds a dendrogram representing the clustering given the current distance matrix
	 * @param distances Matrix of distances
	 * @return Dendrogram representing the clusters
	 */
	public Dendrogram buildDendrogram(DistanceMatrix distances);
}
