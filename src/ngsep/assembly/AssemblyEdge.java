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

import java.io.Serializable;

/**
 * @author Jorge Duitama
 */
public class AssemblyEdge implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -52781446837930525L;
	
	private AssemblyVertex vertex1;
	private AssemblyVertex vertex2;
	private int cost;
	private int overlap;

	public AssemblyEdge(AssemblyVertex vertex1, AssemblyVertex vertex2, int cost, int overlap) {
		this.vertex1 = vertex1;
		this.vertex2 = vertex2;
		this.cost = cost;
		this.overlap = overlap;
	}

	/**
	 * @return the vertex1
	 */
	public AssemblyVertex getVertex1() {
		return vertex1;
	}

	/**
	 * @return the vertex2
	 */
	public AssemblyVertex getVertex2() {
		return vertex2;
	}

	/**
	 * @return the cost
	 */
	public int getCost() {
		return cost;
	}

	/**
	 * @param cost the cost to set
	 */
	public void setCost(int cost) {
		this.cost = cost;
	}

	/**
	 * @return the overlap
	 */
	public int getOverlap() {
		return overlap;
	}

	/**
	 * @param overlap the overlap to set
	 */
	public void setOverlap(int overlap) {
		this.overlap = overlap;
	}

	public AssemblyVertex getConnectingVertex(AssemblyVertex vertex) {
		if(vertex1==vertex) return vertex2;
		if(vertex2==vertex) return vertex1;
		return null;
	}
	
	public AssemblyVertex getSharedVertex (AssemblyEdge edge2) {
		if(vertex1==edge2.vertex1 || vertex1 == edge2.vertex2) return vertex1;
		if(vertex2==edge2.vertex1 || vertex2 == edge2.vertex2) return vertex2;
		return null;
	}
	
	

	

}
