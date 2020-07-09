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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import ngsep.math.Distribution;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class LayoutBuilderModifiedKruskal implements LayoutBuilder {

	@Override
	public void findPaths(AssemblyGraph graph) {
		List<List<AssemblyEdge>> edgesSTConnectedComponents = buildSpanningTree(graph.getEdges());
		System.out.println("Number of connected components after spanning tree: "+edgesSTConnectedComponents.size());
		for(List<AssemblyEdge> edges:edgesSTConnectedComponents) {
			if(edges.size()<=1) continue;
			buildPaths (graph, edges);
			//List<List<AssemblyEdge>> subClusters = removeEdgesHighDegreeVertices(edges);
			//for(List<AssemblyEdge> edges2:subClusters) {
				//if(edges.size()<=1) continue;
				//buildPaths (graph, edges2);
			//}		
		}	
	}

	private List<List<AssemblyEdge>> buildSpanningTree(List<AssemblyEdge> initialEdges) {
		List<AssemblyEdge> edgesDifferentSequences = new ArrayList<AssemblyEdge>();
		//Set of connected vertices indexed by each vertex  
		Map<Integer, Set<Integer>> vertexSets = new HashMap<Integer, Set<Integer>>();
		//Map of connecting edges, also indexed by each vertex
		Map<Integer, List<AssemblyEdge>> forest = new HashMap<Integer, List<AssemblyEdge>>();
		
		for(AssemblyEdge edge:initialEdges) {
			int i1 = edge.getVertex1().getUniqueNumber();
			int i2 = edge.getVertex2().getUniqueNumber();
			if(edge.isSameSequenceEdge()) {
				Set<Integer> singleEdgeSet = new TreeSet<Integer>();
				singleEdgeSet.add(i1);
				singleEdgeSet.add(i2);
				vertexSets.put(i1, singleEdgeSet);
				vertexSets.put(i2, singleEdgeSet);
				List<AssemblyEdge> singleEdgeList = new ArrayList<AssemblyEdge>();
				singleEdgeList.add(edge);
				forest.put(i1, singleEdgeList);
				forest.put(i2, singleEdgeList);
			} else {
				edgesDifferentSequences.add(edge);
			}
		}
		System.out.println("Edges different sequences: "+edgesDifferentSequences.size());
		Collections.sort(edgesDifferentSequences, (e1,e2)-> e1.getCost()-e2.getCost());
		for(AssemblyEdge edge: edgesDifferentSequences) {
			int i1 = edge.getVertex1().getUniqueNumber();
			int i2 = edge.getVertex2().getUniqueNumber();
			Set<Integer> connected1 = vertexSets.get(i1);
			Set<Integer> connected2 = vertexSets.get(i2);
			List<AssemblyEdge> edges1 = forest.get(i1);
			List<AssemblyEdge> edges2 = forest.get(i2);
			if(connected1==null) {
				//This should never happen because all vertices should at least be connected with its same sequence vertex
				System.err.println("WARN: Same sequence edge absent for vertex: "+edge.getVertex1().getSequenceIndex()+"-"+edge.getVertex1().isStart());
				continue;
			}
			if(connected2==null) {
				//This should never happen because all vertices should at least be connected with its same sequence vertex
				System.err.println("WARN: Same sequence edge absent for vertex: "+edge.getVertex2().getSequenceIndex()+"-"+edge.getVertex2().isStart());
				continue;
			}
			if(connected1!=connected2) {
				if(!Collections.disjoint(connected1, connected2)) {
					System.err.println("WARN: Intersection detected between different vertex sets");
					continue;
				}
				connected1.addAll(connected2);
				edges1.addAll(edges2);
				edges1.add(edge);
				for(int i:connected1) {
					vertexSets.put(i, connected1);
					forest.put(i, edges1);
				}
			}
		}
		List<List<AssemblyEdge>> finalList = new ArrayList<List<AssemblyEdge>>();
		Set<Integer> unionVertexSets = new HashSet<Integer>();
		for(int i:vertexSets.keySet()) {
			if(unionVertexSets.contains(i)) continue;
			Set<Integer> connected = vertexSets.get(i);
			unionVertexSets.addAll(connected);
			List<AssemblyEdge> nextList = forest.get(i);
			System.out.println("Adding next tree with: "+nextList.size()+" edges");
			finalList.add(nextList);
		}
		Collections.sort(finalList, (l1,l2)-> l2.size()-l1.size());
		return finalList;
	}

	private List<List<AssemblyEdge>> removeEdgesHighDegreeVertices(List<AssemblyEdge> edges) {
		Map<Integer, List<AssemblyEdge>> adjacencyList = buildAdjacencyList(edges);
		Distribution degreeDistribution = new Distribution(0, 10, 1);
		for(List<AssemblyEdge> connected:adjacencyList.values()) degreeDistribution.processDatapoint(connected.size());
		System.out.println("Processing list of "+edges.size()+" edges. Degree distribution");
		degreeDistribution.printDistributionInt(System.out);
		// Edges to remove indexed by vertex number
		Map<Integer, List<AssemblyEdge>>  toRemove = new HashMap<Integer, List<AssemblyEdge>>();
		for(int i:adjacencyList.keySet()) {
			List<AssemblyEdge> edgesI = adjacencyList.get(i);
			if(edgesI.size()<3) {
				continue;
			}
			List<AssemblyEdge> toRemoveI = chooseEdgesToRemove (edgesI, adjacencyList);
			for(AssemblyEdge edge: toRemoveI) {
				int j = edge.getVertex1().getUniqueNumber();
				List<AssemblyEdge> toRemoveJ = toRemove.computeIfAbsent(j, l-> new ArrayList<AssemblyEdge>());
				toRemoveJ.add(edge);
			}
		}
		List<AssemblyEdge> finalEdges = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge edge:edges) {
			int i = edge.getVertex1().getUniqueNumber();
			List<AssemblyEdge> toRemoveI = toRemove.get(i);
			if(toRemoveI==null || !toRemoveI.contains(edge)) finalEdges.add(edge);
		}
		System.out.println("Final number of edges: "+finalEdges.size());
		return buildSpanningTree(finalEdges);
	}

	private Map<Integer, List<AssemblyEdge>> buildAdjacencyList(List<AssemblyEdge> edges) {
		Map<Integer, List<AssemblyEdge>> answer = new HashMap<Integer, List<AssemblyEdge>>();
		for(AssemblyEdge edge:edges) {
			int i1 = edge.getVertex1().getUniqueNumber();
			int i2 = edge.getVertex2().getUniqueNumber();
			List<AssemblyEdge> l1 = answer.computeIfAbsent(i1, l->new ArrayList<AssemblyEdge>());
			l1.add(edge);
			List<AssemblyEdge> l2 = answer.computeIfAbsent(i2, l->new ArrayList<AssemblyEdge>());
			l2.add(edge);
		}
		return answer;
	}

	private List<AssemblyEdge> chooseEdgesToRemove(List<AssemblyEdge> edges, Map<Integer, List<AssemblyEdge>> adjacencyList) {
		List<AssemblyEdge> toRemove = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			int i1 = edge.getVertex1().getUniqueNumber();
			int i2 = edge.getVertex2().getUniqueNumber();
			//TODO: Improve this algorithm. By now it only would remove edges connecting two high degree nodes
			if(adjacencyList.get(i1).size()<3) continue;
			if(adjacencyList.get(i2).size()<3) continue;
			toRemove.add(edge);
		}
		return toRemove;
	}

	private void buildPaths(AssemblyGraph graph, List<AssemblyEdge> edges) {
		Map<Integer, List<AssemblyEdge>> adjacencyList = buildAdjacencyList(edges);
		System.out.println("Building paths from tree built from "+edges.size()+" edges. Connected vertices: "+adjacencyList.size());
		Set<Integer> addedVertices = new HashSet<Integer>();
		while(addedVertices.size() < adjacencyList.size()) {
			List<AssemblyEdge> pathMaxLength = null;
			int maxLength = 0;
			for(int i:adjacencyList.keySet()) {
				if(addedVertices.contains(i)) continue;
				List<AssemblyEdge> edgesI = adjacencyList.get(i);
				if(edgesI.size()==1 && edgesI.get(0).isSameSequenceEdge()) {
					List<AssemblyEdge> path = buildNewPath(adjacencyList, i, edgesI.get(0),addedVertices);
					if(pathMaxLength==null || maxLength<path.size()) {
						pathMaxLength = path;
						maxLength = path.size();
					}
					
				}
			}
			if(pathMaxLength==null || maxLength<2) break;
			for(AssemblyEdge edge:pathMaxLength) {
				addedVertices.add(edge.getVertex1().getUniqueNumber());
				addedVertices.add(edge.getVertex2().getUniqueNumber());
			}
			if(pathMaxLength.size()>1) {	
				System.out.println("Adding path of length: "+pathMaxLength.size());
				graph.addPath(pathMaxLength);
			}
			
		}
	}
	
	public List<AssemblyEdge> buildNewPath (Map<Integer, List<AssemblyEdge>> adjacencyList, int start, AssemblyEdge startEdge, Set<Integer> addedVertices) {
		List<AssemblyEdge> newPath = new ArrayList<AssemblyEdge>();
		Set<Integer> addedVertices2 = new HashSet<Integer>();
		addedVertices2.addAll(addedVertices);
		newPath.add(startEdge);
		addedVertices2.add(start);
		AssemblyVertex nextVertex = startEdge.getVertex1();
		if(nextVertex.getUniqueNumber()==start) nextVertex = startEdge.getVertex2();
		addedVertices2.add(nextVertex.getUniqueNumber());
		AssemblyEdge nextEdge = chooseNextEdge(nextVertex, adjacencyList.get(nextVertex.getUniqueNumber()),addedVertices2, false);
		while(nextEdge !=null) {
			newPath.add(nextEdge);
			nextVertex = nextEdge.getConnectingVertex(nextVertex);
			int i = nextVertex.getUniqueNumber();
			addedVertices2.add(i);
			//Add next same sequence edge
			nextEdge = chooseNextEdge(nextVertex, adjacencyList.get(i),addedVertices2, true);
			newPath.add(nextEdge);
			nextVertex = nextEdge.getConnectingVertex(nextVertex);
			int i2 = nextVertex.getUniqueNumber();
			addedVertices2.add(i2);
			nextEdge = chooseNextEdge(nextVertex, adjacencyList.get(i2),addedVertices2, false);
		}
		return newPath;
	}

	private AssemblyEdge chooseNextEdge(AssemblyVertex vertex, List<AssemblyEdge> edges, Set<Integer> addedVertices, boolean sameSequence) {
		AssemblyEdge answer = null;
		int minCost = 0;
		if(edges==null) return null;
		for(AssemblyEdge edge:edges) {
			if (sameSequence && edge.isSameSequenceEdge()) return edge;
			if (edge.isSameSequenceEdge()) continue;
			AssemblyVertex connecting = edge.getConnectingVertex(vertex);
			//System.out.println("Next Vertex:"+vertex.getUniqueNumber()+" Edge v1: "+edge.getVertex1().getUniqueNumber()+" v2: "+edge.getVertex2().getSequenceIndex()+" vertex: "+vertex.getSequenceIndex()+" connecting "+connecting);
			if (addedVertices.contains(connecting.getUniqueNumber())) continue;
			if(answer ==null || minCost> edge.getCost()) {
				answer = edge;
				minCost = edge.getCost();
			}
			
		}
		return answer;
	}

}
