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
 * @author Laura Gonzalez
 *
 */
public class SyntenyVertex {
	private LocalHomologyCluster localRegion1;
	private LocalHomologyCluster localRegion2;
	public SyntenyVertex(LocalHomologyCluster localRegion1, LocalHomologyCluster localRegion2) {
		super();
		this.localRegion1 = localRegion1;
		this.localRegion2 = localRegion2;
	}
	public HomologyCluster getHomologyCluster() {
		return localRegion1.getParent();
	}
	public LocalHomologyCluster getLocalRegion1() {
		return localRegion1;
	}
	public LocalHomologyCluster getLocalRegion2() {
		return localRegion2;
	}
	public double getMaximumEdgePCTSharedKmers() {
		List<HomologyUnit> units1 = localRegion1.getHomologyUnitsCluster();
		List<HomologyUnit> units2 = localRegion2.getHomologyUnitsCluster();
		double maxPCT = 0;
		for(HomologyUnit u1:units1) {
			for(HomologyUnit u2: units2) {
				HomologyEdge edge1 = u1.getHomologyEdge(u2);
				if(edge1!=null) maxPCT = Math.max(maxPCT, edge1.getScore());
				HomologyEdge edge2 = u2.getHomologyEdge(u1);
				if(edge2!=null) maxPCT = Math.max(maxPCT, edge2.getScore());
			}
			
		}
		return maxPCT;
	}
	
}
