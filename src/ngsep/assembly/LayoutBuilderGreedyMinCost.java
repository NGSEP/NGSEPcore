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

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * @author Jorge Duitama
 */
public class LayoutBuilderGreedyMinCost implements LayourBuilder {

	@Override
	public void findPaths(AssemblyGraph graph) {
		Set<Integer> sequencesInPaths = new HashSet<>();
		List<AssemblyVertex> vertices = graph.getVertices();
		for(int i=0;i<vertices.size();i++) {
			AssemblyVertex vertex = findNextUncoveredVertex(vertices,sequencesInPaths);
			if(vertex==null) return;
			
			
			LinkedList<AssemblyEdge> currentPath = new LinkedList<AssemblyEdge>();
			AssemblyEdge nextEdgeLeft = graph.getSameSequenceEdge (vertex);
			AssemblyEdge nextEdgeRight = nextEdgeLeft;
			currentPath.add(nextEdgeRight);
			sequencesInPaths.add(vertex.getSequenceIndex());
			AssemblyVertex vertexLeft = nextEdgeRight.getVertex1();
			nextEdgeLeft = findBestUncoveredEdge(graph, vertexLeft, sequencesInPaths);
			AssemblyVertex vertexRight = nextEdgeRight.getVertex2();
			nextEdgeRight = findBestUncoveredEdge(graph, vertexRight, sequencesInPaths);
			while(nextEdgeLeft!=null || nextEdgeRight!=null) {
				if(nextEdgeRight==null || (nextEdgeLeft!=null && nextEdgeLeft.getCost()<nextEdgeRight.getCost())) {
					currentPath.add(0,nextEdgeLeft);
					vertexLeft = nextEdgeLeft.getConnectingVertex(vertexLeft);
					nextEdgeLeft = graph.getSameSequenceEdge (vertexLeft);
					currentPath.add(0,nextEdgeLeft);
					vertexLeft = nextEdgeLeft.getConnectingVertex(vertexLeft);
					sequencesInPaths.add(vertexLeft.getSequenceIndex());
				} else if(nextEdgeRight!=null) {
					currentPath.add(nextEdgeRight);
					vertexRight = nextEdgeRight.getConnectingVertex(vertexRight);
					nextEdgeRight = graph.getSameSequenceEdge (vertexRight);
					currentPath.add(nextEdgeRight);
					vertexRight = nextEdgeRight.getConnectingVertex(vertexRight);
					sequencesInPaths.add(vertexRight.getSequenceIndex());
				}
				nextEdgeLeft = findBestUncoveredEdge(graph, vertexLeft, sequencesInPaths);
				nextEdgeRight = findBestUncoveredEdge(graph, vertexRight, sequencesInPaths);
				//System.out.println("Vertex left "+vertexLeft.getIndex()+" vertex right: "+vertexRight.getIndex());
			}
			if(currentPath.size()>1) {
				System.out.println("Found path of size "+currentPath.size());
				//printPath(currentPath);
				graph.addPath(currentPath);
			}
		}
		
		
		
		
	}

	public void printPath(LinkedList<AssemblyEdge> path) {
		for(AssemblyEdge edge:path) {
			AssemblyVertex v1 = edge.getVertex1();
			AssemblyVertex v2 = edge.getVertex2();
			System.out.println("Edge between "+v1.getSequenceIndex()+"-"+v1.isStart()+" and "+v2.getSequenceIndex()+"-"+v2.isStart());
		}	
	}

	private AssemblyVertex findNextUncoveredVertex(List<AssemblyVertex> vertices, Set<Integer> sequencesInPaths) {
		for(AssemblyVertex vertex:vertices) {
			if(!sequencesInPaths.contains(vertex.getSequenceIndex())) return vertex;
		}
		return null;
	}
	private AssemblyEdge findBestUncoveredEdge(AssemblyGraph graph, AssemblyVertex vertex, Set<Integer> sequencesInPaths) {
		List<AssemblyEdge> edgesVertex = graph.getEdges(vertex);
		int minCost = 0;
		AssemblyEdge minEdge = null;
		for(AssemblyEdge edge:edgesVertex) {
			AssemblyVertex connectingVertex = edge.getConnectingVertex(vertex);
			if(!sequencesInPaths.contains(connectingVertex.getSequenceIndex())) {
				if(minEdge == null || minCost>edge.getCost()) {
					minCost = edge.getCost();
					minEdge = edge;
				}
			}
		}
		return minEdge;
	}

}
