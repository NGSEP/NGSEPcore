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
package ngsep.assembly;

import java.util.Comparator;

/**
 * @author Jorge Duitama
 */
public class LayoutBuilderGreedyMinCost extends LayoutBuilderGreedy {

	@Override
	public void findPaths(AssemblyGraph graph) {
		Comparator<AssemblyEdge> comparator = new Comparator<AssemblyEdge>() {
			@Override
			public int compare(AssemblyEdge edge1, AssemblyEdge edge2) {
				return edge1.getCost()-edge2.getCost();
			}
		};
		setComparator(comparator);
		super.findPaths(graph);
	}
}
