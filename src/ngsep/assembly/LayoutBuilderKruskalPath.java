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
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

import JSci.maths.statistics.NormalDistribution;
import ngsep.math.Distribution;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class LayoutBuilderKruskalPath implements LayoutBuilder {
	
	private Logger log = Logger.getLogger(LayoutBuilderKruskalPath.class.getName());

	private int minPathLength = 6;
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	public int getMinPathLength() {
		return minPathLength;
	}
	public void setMinPathLength(int minPathLength) {
		this.minPathLength = minPathLength;
	}
	
	@Override
	public void findPaths(AssemblyGraph graph) {
		List<AssemblyEdge> pathEdges = graph.selectSafeEdges();
		log.info("Number of safe edges: "+pathEdges.size());
		
		List<LinkedList<AssemblyEdge>> safePaths = buildPaths(graph, pathEdges);
		//Map<Integer,Integer> vertexPathIds = calculateVertexPathIdsMap(safePaths);
		log.info("Number of paths safe edges: "+safePaths.size());
		
		//Algorithms to resolve conflicts between almost safe close edges
		//addEdges2(graph, safePaths, pathEdges);
		//log.info("Number of edges after second algorithm: "+pathEdges.size());
		//safePaths = buildPaths(graph, pathEdges);
		//log.info("Updated number of paths: "+safePaths.size());
		
		//addEdges3(graph, safePaths, pathEdges);
		//log.info("Number of edges after third algorithm: "+pathEdges.size());
		//safePaths = buildPaths(graph, pathEdges);
		//log.info("Updated number of paths: "+safePaths.size());
		
		addConnectingEdges(graph, safePaths, pathEdges);
		List<LinkedList<AssemblyEdge>> paths = buildPaths(graph, pathEdges);
		log.info("Paths costs algorithm: "+paths.size());
		List<List<AssemblyEdge>> finalPaths = mergeClosePaths(graph, paths, pathEdges);
		for(List<AssemblyEdge> path:finalPaths) {
			if(path.size()<minPathLength) continue;
			graph.addPath(path);
		}
		log.info("Final number of paths: "+graph.getPaths().size());
		System.out.println("Estimated N statistics");
		long [] stats = graph.estimateNStatisticsFromPaths();
		if(stats!=null) NStatisticsCalculator.printNStatistics(stats, System.out);
	}

	private List<LinkedList<AssemblyEdge>> buildPaths(AssemblyGraph graph, List<AssemblyEdge> edges) {
		Map<Integer,AssemblyEdge> edgesByVertex = new HashMap<Integer, AssemblyEdge>(); 
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			AssemblyVertex v1 = edge.getVertex1();
			if(edgesByVertex.containsKey(v1.getUniqueNumber())) System.err.println("WARN: two edges for vertex. E1 "+edgesByVertex.get(v1.getUniqueNumber())+" E2: "+edge);
			AssemblyVertex v2 = edge.getVertex2();
			if(edgesByVertex.containsKey(v2.getUniqueNumber())) System.err.println("WARN: two edges for vertex. E1 "+edgesByVertex.get(v2.getUniqueNumber())+" E2: "+edge);
			edgesByVertex.put(v1.getUniqueNumber(), edge);
			edgesByVertex.put(v2.getUniqueNumber(), edge);
		}
		int n = graph.getNumSequences();
		List<LinkedList<AssemblyEdge>> paths = new ArrayList<LinkedList<AssemblyEdge>>();
		Set<Integer> sequencesInPaths = new HashSet<>();
		for(int i=0;i<n;i++) {
			if(graph.getVertex(i, true)==null || graph.getVertex(i, false)==null) continue;
			if(graph.isEmbedded(i)) continue;
			
			if(sequencesInPaths.contains(i)) continue;
			AssemblyVertex nextVertex = graph.getVertex(i, true);
			if(graph.getEdges(nextVertex).size()==0) continue;
			LinkedList<AssemblyEdge> nextPath = new LinkedList<AssemblyEdge>();
			nextPath.add(graph.getSameSequenceEdge(nextVertex));
			//Expand v1
			while(nextVertex!=null) {
				AssemblyEdge nextEdge = edgesByVertex.get(nextVertex.getUniqueNumber());
				if(nextEdge == null) nextVertex=null;
				else {
					nextVertex = nextEdge.getConnectingVertex(nextVertex);
					if(sequencesInPaths.contains(nextVertex.getSequenceIndex())) {
						System.err.println("WARN: Cycle detected building paths. Next edge: "+nextEdge.getVertex1().getRead().getName()+" "+nextEdge.getVertex2().getRead().getName());
						break;
					}
					nextPath.add(0, nextEdge);
					AssemblyEdge seqEdge = graph.getSameSequenceEdge(nextVertex);
					nextPath.add(0, seqEdge);
					sequencesInPaths.add(nextVertex.getSequenceIndex());
					nextVertex = seqEdge.getConnectingVertex(nextVertex);
				}
			}
			//Expand v2
			nextVertex = graph.getVertex(i, false);
			while(nextVertex!=null) {
				AssemblyEdge nextEdge = edgesByVertex.get(nextVertex.getUniqueNumber());
				if(nextEdge == null) nextVertex=null;
				else {
					nextVertex = nextEdge.getConnectingVertex(nextVertex);
					if(sequencesInPaths.contains(nextVertex.getSequenceIndex())) {
						System.err.println("WARN: Cycle detected building paths. Next edge: "+nextEdge.getVertex1().getRead().getName()+" "+nextEdge.getVertex2().getRead().getName());
						break;
					}
					nextPath.add(nextEdge);
					AssemblyEdge seqEdge = graph.getSameSequenceEdge(nextVertex);
					nextPath.add(seqEdge);
					sequencesInPaths.add(nextVertex.getSequenceIndex());
					nextVertex = seqEdge.getConnectingVertex(nextVertex);
				}
			}
			paths.add(nextPath);
		} 
		return paths;
	}
 
	
	private AssemblyVertex [] extractEndVertices (List<LinkedList<AssemblyEdge>> paths) {
		int p = paths.size();
		AssemblyVertex [] vertices = new AssemblyVertex[2*p];
		int v=0;
		for(int i=0;i<p;i++) {
			List<AssemblyEdge> path = paths.get(i);
			if(path.size()==1) {
				vertices[v] = path.get(0).getVertex1();
				v++;
				vertices[v] = path.get(0).getVertex2();
				v++;
				continue;
			}
			AssemblyEdge edge0 = path.get(0);
			AssemblyEdge edge1 = path.get(1);
			AssemblyVertex v1 = edge0.getSharedVertex(edge1);
			vertices[v] = edge0.getConnectingVertex(v1);
			v++;
			AssemblyEdge edgef = path.get(path.size()-1);
			AssemblyEdge edgef2 = path.get(path.size()-2);
			AssemblyVertex v2 = edgef.getSharedVertex(edgef2);
			vertices[v] = edgef.getConnectingVertex(v2);
			v++;
		}
		return vertices;
	}
	private void addConnectingEdges(AssemblyGraph graph, List<LinkedList<AssemblyEdge>> paths, List<AssemblyEdge> pathEdges) {
		NormalDistribution distIKBP = graph.estimateDistributions(pathEdges, new HashSet<Integer>())[5];
		AssemblyVertex [] vertices = extractEndVertices(paths);
		log.info("KruskalPathAlgorithm. Extracted "+vertices.length+" end vertices");
		List<AssemblyEdge> candidateEdges = new ArrayList<AssemblyEdge>();
		for(int i=0;i<vertices.length;i++) {
			List<AssemblyEdge> edgesVI = graph.getEdges(vertices[i]);
			for(AssemblyEdge edge:edgesVI) {
				if(edge.isSameSequenceEdge()) continue;
				if(graph.isEmbedded(edge.getVertex1().getSequenceIndex())) continue;
				if(graph.isEmbedded(edge.getVertex2().getSequenceIndex())) continue;
				//Add each candidate vertex only one time
				candidateEdges.add(edge);
			}
		}
		log.info("KruskalPathAlgorithm. selected "+candidateEdges.size()+" candidate edges");
		Collections.sort(candidateEdges,(e1,e2)->e1.getCost()-e2.getCost());
		//Collections.sort(candidateEdges,(e1,e2)->e2.getScore()-e1.getScore());
		log.info("KruskalPathAlgorithm. Sorted "+candidateEdges.size()+" edges");
		List<AssemblyEdge> selectedEdges = selectEdgesToMergePaths(candidateEdges,vertices, distIKBP);
		log.info("KruskalPathAlgorithm. Selected "+selectedEdges.size()+" edges for paths");
		pathEdges.addAll(selectedEdges);
	}
	private List<AssemblyEdge> selectEdgesToMergePaths(List<AssemblyEdge> candidateEdges, AssemblyVertex [] vertices, NormalDistribution distIKBP) {
		double limitIKBP = distIKBP.getMean()+7*Math.sqrt(distIKBP.getVariance());
		log.info("Limit for IKBP: "+limitIKBP+" Average: "+distIKBP.getMean()+" variance: "+distIKBP.getVariance());
		int n = vertices.length;
		Map<Integer,Integer> verticesPos = new HashMap<Integer, Integer>();
		for(int i=0;i<n;i++) {
			verticesPos.put(vertices[i].getUniqueNumber(), i);
		}
		int [] clusters = new int[n];
		boolean [] used = new boolean[n];
		Arrays.fill(used, false);
		for(int i=0;i<n;i++) {
			clusters[i] = i/2;
		}
		List<AssemblyEdge> answer = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge nextEdge:candidateEdges) {
			//if(nextEdge.getEdgeAssemblyGraph().getVertex1().getUniqueNumber()==2854) System.out.println("score: "+calculateCost(nextEdge, edgesStats)+" edge: "+nextEdge.getEdgeAssemblyGraph());
			//if(answer.size()<20) System.out.println("score: "+calculateCost(nextEdge, edgesStats)+ " edge: "+nextEdge.getEdgeAssemblyGraph());
			AssemblyVertex v1 = nextEdge.getVertex1();
			AssemblyVertex v2 = nextEdge.getVertex2();
			if(v1==null || v2 == null) continue;
			Integer posV1 = verticesPos.get(v1.getUniqueNumber());
			Integer posV2 = verticesPos.get(v2.getUniqueNumber());
			if(posV1==null || posV2==null) continue;
			if(v1.getUniqueNumber()==-97856 || v2.getUniqueNumber()==69473) System.out.println("SelectingEdges. next edge: "+nextEdge+" used: "+used[posV1]+" "+used[posV2]+" clusters: "+clusters[posV1]+" "+clusters[posV2]);
			
			if(used[posV1] || used[posV2]) continue;
			if(nextEdge.getIndelsPerKbp()>limitIKBP) continue;
			int c1 = clusters[posV1];
			int c2 = clusters[posV2];
			
			if(c1!=c2) {
				answer.add(nextEdge);
				used[posV1] =used[posV2] = true;
				for(int i=0;i<n;i++) {
					if(clusters[i]==c2) clusters[i] = c1;
				}
			}
		}
		return answer;
	}
	private List<List<AssemblyEdge>> mergeClosePaths (AssemblyGraph graph, List<LinkedList<AssemblyEdge>> paths, List<AssemblyEdge> pathEdges) {
		Map<Integer,VertexPathLocation> vertexPositions = new HashMap<Integer,VertexPathLocation>();
		Map<Integer,Integer> vertexEnds = new HashMap<Integer, Integer>();
		Distribution costsDistribution = new Distribution(0, 1000000, 1000);
		for(int i=0;i<paths.size();i++) {
			List<AssemblyEdge> path = new ArrayList<AssemblyEdge>(paths.get(i));
			int pathId = i+1;
			int n = path.size();
			if(n<5) continue;
			AssemblyEdge edge0 = path.get(0);
			AssemblyEdge edge1 = path.get(1);
			AssemblyVertex lastVertex = edge0.getSharedVertex(edge1);
			AssemblyVertex firstVertexOut = edge0.getConnectingVertex(lastVertex);
			vertexEnds.put(firstVertexOut.getUniqueNumber(),pathId);
			vertexPositions.put(firstVertexOut.getUniqueNumber(), new VertexPathLocation(firstVertexOut.getUniqueNumber(), pathId, n, 0));
			vertexPositions.put(lastVertex.getUniqueNumber(), new VertexPathLocation(lastVertex.getUniqueNumber(), pathId, n, 1));
			int k = 2;
			for(int j=1;j<path.size();j++) {
				AssemblyEdge edge = path.get(j);
				if(!edge.isSameSequenceEdge()) costsDistribution.processDatapoint(edge.getCost());
				AssemblyVertex nextVertex = edge.getConnectingVertex(lastVertex);
				vertexPositions.put(nextVertex.getUniqueNumber(), new VertexPathLocation(nextVertex.getUniqueNumber(), pathId, n, k));
				if(k==2 || k==n-2) vertexEnds.put(nextVertex.getUniqueNumber(),pathId);
				k++;
				lastVertex = nextVertex;
			}
			vertexEnds.put(lastVertex.getUniqueNumber(),pathId);
		}
		log.info("Average cost path edges: "+costsDistribution.getAverage());
		Map<String,PathEndJunctionEdge> pathEndEdges = new HashMap<String, PathEndJunctionEdge>();
		
		for(int vertexId:vertexEnds.keySet()) {
			AssemblyVertex v1 = graph.getVertexByUniqueId(vertexId);
			VertexPathLocation v1Loc = vertexPositions.get(v1.getUniqueNumber());
			AssemblyEdge edge = getEdgeBestOverlap(graph,v1);
			//if(vertexId==-32163) System.out.println("FInding edges to merge paths. Edge: "+edge);
			if(edge == null) continue;
			if(edge.getCost()>2*costsDistribution.getAverage()) continue;
			AssemblyVertex v2 = edge.getConnectingVertex(v1);
			VertexPathLocation v2Loc = vertexPositions.get(v2.getUniqueNumber());
			if(v2Loc == null) continue;
			if (v1Loc.getPathId() == v2Loc.getPathId()) continue;
			int pathEndV1 = v1Loc.getPathId();
			if(v1Loc.getPathPosition()>v1Loc.getPathLength()-4) pathEndV1 = -pathEndV1;
			else if (v1Loc.getPathPosition()>4) continue;
			int pathEndV2 = v2Loc.getPathId();
			if(v2Loc.getPathPosition()>v2Loc.getPathLength()-4) pathEndV2 = -pathEndV2;
			else if (v2Loc.getPathPosition()>4) continue;
			int minId = Math.min(pathEndV1, pathEndV2);
			int maxId = Math.max(pathEndV1, pathEndV2);
			String key = "F"+minId+"L"+maxId;
			PathEndJunctionEdge pathEndEdge = pathEndEdges.computeIfAbsent(key, v->new PathEndJunctionEdge(minId, maxId));
			log.info("Adding vote to relationship between path ends "+minId+" "+maxId+" edge: "+edge);
			pathEndEdge.addVote(edge.getCost());
		}
		log.info("Found "+pathEndEdges.size()+" candidate edges connecting vertices");
		Map<Integer,List<Integer>> mergedPathIds = findPathsToMerge(pathEndEdges.values(),paths.size());
		log.info("Merged path ids: "+mergedPathIds.values());
		List<List<AssemblyEdge>> answer = new ArrayList<List<AssemblyEdge>>(mergedPathIds.size());
		
		for(List<Integer> idsPath:mergedPathIds.values()) {
			List<AssemblyEdge> nextPath = new ArrayList<AssemblyEdge>();
			int lastEndId = 0;
			int lastInputPathMerged = 0;
			for(int pathEndId:idsPath) {
				if(pathEndId!=-lastEndId) {
					lastEndId = pathEndId;
					continue;
				}
				if(pathEndId > 0) {
					LinkedList<AssemblyEdge> nextInputPath = paths.get(pathEndId-1);
					Collections.reverse(nextInputPath);
					if(!mergePaths(graph, nextPath,nextInputPath, lastInputPathMerged, lastEndId)) {
						answer.add(nextPath);
						nextPath = new ArrayList<AssemblyEdge>();
						nextPath.addAll(nextInputPath);
					}
				} else {
					LinkedList<AssemblyEdge> nextInputPath = paths.get(lastEndId-1);
					if(!mergePaths(graph, nextPath,nextInputPath, lastInputPathMerged, lastEndId)) {
						answer.add(nextPath);
						nextPath = new ArrayList<AssemblyEdge>();
						nextPath.addAll(nextInputPath);
					}
				}
				lastInputPathMerged = pathEndId;
				lastEndId = pathEndId;
			}
			//log.info("Next final path size: "+nextPath.size());
			if(nextPath.size()>0) answer.add(nextPath);
		}
		
		return answer;
	}
	
	
	
	private boolean mergePaths(AssemblyGraph graph, List<AssemblyEdge> path1, List<AssemblyEdge> path2, int pathEndId1, int pathEndId2) {
		int n = path1.size();
		//if(pathEndId2==1) System.out.println("Merging paths. End ids: "+pathEndId1+" "+pathEndId2+" sizes: "+path1.size()+" "+path2.size());
		if(n==0) {
			path1.addAll(path2);
			return true;
		}
		if(path2.size()==0) return true;
		AssemblyEdge lastEdge1 = path1.get(n-1);
		AssemblyEdge secondLastEdge1 = path1.get(n-2);
		AssemblyVertex connecting1 = lastEdge1.getSharedVertex(secondLastEdge1);
		AssemblyVertex lastVertex1 = lastEdge1.getConnectingVertex(connecting1);
		AssemblyVertex secondLastVertex1 = secondLastEdge1.getConnectingVertex(connecting1);
		int lastOverlap = secondLastEdge1.getOverlap();
		int lastCost = secondLastEdge1.getCost();
		log.info("Last vertices path 1: "+secondLastVertex1+" "+lastVertex1+" path ends: "+pathEndId1+" "+pathEndId2);
		
		AssemblyEdge firstEdge2 = path2.get(0);
		AssemblyEdge secondEdge2 = path2.get(1);
		AssemblyVertex connecting2 = firstEdge2.getSharedVertex(secondEdge2);
		AssemblyVertex firstVertex2 = firstEdge2.getConnectingVertex(connecting2);
		AssemblyVertex secondVertex2 = secondEdge2.getConnectingVertex(connecting2);
		int nextOverlap = secondEdge2.getOverlap();
		int nextCost = secondEdge2.getCost();
		log.info("first vertices path 2: "+firstVertex2+" "+secondVertex2+" path ends: "+pathEndId1+" "+pathEndId2);
		
		
		AssemblyEdge c1 = graph.getEdge(lastVertex1, firstVertex2);
		if(c1!=null && c1.getCost()<2*lastCost && c1.getCost()<2*nextCost && c1.getOverlap()>lastOverlap && c1.getOverlap()>nextOverlap) {
			path1.add(c1);
			path1.addAll(path2);
			return true;
		}
		AssemblyEdge c2 = graph.getEdge(lastVertex1, secondVertex2);
		AssemblyEdge c3 = graph.getEdge(secondLastVertex1, firstVertex2);
		if(c2!=null && c3!=null && c2.getOverlap()<c3.getOverlap()) {
			AssemblyEdge t = c2;
			c2 = c3;
			c3 = t;
		}
		log.info("Edges to merge paths: "+c2+" "+c3+" path ends: "+pathEndId1+" "+pathEndId2 +" overlaps: "+lastOverlap+" "+nextOverlap+" costs: "+lastCost+" "+nextCost);
		if(c2!=null && c2.getOverlap()>nextOverlap && c2.getCost()<2*lastCost && c2.getCost()<2*nextCost && (c3==null || c2.getCost()<2*c3.getCost())) {
			log.info("Merging paths with ends: "+pathEndId1+" "+pathEndId2+" Merging edge: "+c2);
			//if(c3==null) {
				path2.remove(0);
				path2.remove(0);
				path1.add(c2);
				path1.addAll(path2);
				//return;
			//}
			return true;
		}
		if(c3!=null && c3.getOverlap()>lastOverlap && c3.getCost()<2*lastCost && c3.getCost()<2*nextCost && (c2==null || c3.getCost()<2*c2.getCost())) {
			log.info("Merging paths with ends: "+pathEndId1+" "+pathEndId2+" Merging edge: "+c3);
			//if(c3==null) {
				path1.remove(n-1);
				path1.remove(n-2);
				path1.add(c3);
				path1.addAll(path2);
				//return;
			//}
			return true;
		}
		return false;
		
	}
	private AssemblyEdge getEdgeBestOverlap(AssemblyGraph graph, AssemblyVertex vertex) {
		List<AssemblyEdge> edges = graph.getEdges(vertex);
		AssemblyEdge maxOverlapEdge = null;
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			if(maxOverlapEdge==null || maxOverlapEdge.getOverlap()<edge.getOverlap()) maxOverlapEdge = edge;
		}
		return maxOverlapEdge;
	}
	private Map<Integer,List<Integer>> findPathsToMerge(Collection<PathEndJunctionEdge> pathEndEdges, int n) {
		List<PathEndJunctionEdge> pathEndEdgesList = new ArrayList<PathEndJunctionEdge>(pathEndEdges);
		Collections.sort(pathEndEdgesList,(e1,e2)->e1.getCost()-e2.getCost());
		Map<Integer,Integer> pathGroups = new HashMap<Integer, Integer>();
		
		Map<Integer,List<Integer>> answer = new TreeMap<Integer,List<Integer>>();
		for(int i=0;i<n;i++) {
			int pathId = i+1;
			pathGroups.put(pathId, pathId);
			pathGroups.put(-pathId, pathId);
			LinkedList<Integer> nextPath = new LinkedList<Integer>();
			nextPath.add(pathId);
			nextPath.add(-pathId);
			answer.put(pathId,nextPath);
		}
		Set<Integer> usedPathEnds = new HashSet<Integer>();
		for(PathEndJunctionEdge pathEdge:pathEndEdgesList) {
			int p1 = pathEdge.getPath1EndId();
			int p2 = pathEdge.getPath2EndId();
			if(usedPathEnds.contains(p1) || usedPathEnds.contains(p2)) continue;
			int g1 = pathGroups.get(p1);
			int g2 = pathGroups.get(p2);
			if(g1==g2) continue;
			//Merge path1 with path2
			
			List<Integer> path1 = answer.get(g1);
			List<Integer> path2 = answer.get(g2);
			if(path2.get(path2.size()-1)==p2) {
				Collections.reverse(path2);
			} else if (path2.get(0)!=p2) {
				continue;
			}
			boolean lastP1 = false;
			if(path1.get(path1.size()-1)==p1) {
				lastP1 = true;
			} else if (path1.get(0)!=p1) {
				continue;
			}
			for(int i:path2) {
				if(lastP1) path1.add(i);
				else path1.add(0, i);
				pathGroups.put(i, g1);
			}
			answer.remove(g2);
			usedPathEnds.add(p1);
			usedPathEnds.add(p2);
		}
		return answer;
	}
	//Other algotrithms not used by now
	private void addEdges3(AssemblyGraph graph, List<LinkedList<AssemblyEdge>> paths, List<AssemblyEdge> pathEdges) {
		AssemblyVertex [] vertices = extractEndVertices(paths);
		int debugSeq = -1;
		int n = vertices.length;
		int [] clusters = new int[n];
		boolean [] used = new boolean[n];
		Arrays.fill(used, false);
		for(int i=0;i<n;i++) {
			clusters[i] = i/2;
		}
		Map<Integer,Integer> verticesMap = new HashMap<Integer, Integer>(n);
		for(int i=0;i<n;i++) verticesMap.put(vertices[i].getUniqueNumber(), i);
		Map<Integer,Set<Integer>> vertexSequences = new HashMap<Integer, Set<Integer>>();
		// Try to expand from vertex i
		for(int i=0;i<vertices.length;i++) {
			if(vertices[i].getSequenceIndex()==debugSeq) System.out.println("AddEdges3. Vertex "+vertices[i]+" used: "+used[i]+" path size: "+paths.get(i/2).size());
			if(used[i]) continue;
			if(paths.get(i/2).size()==1) continue;
			AssemblyVertex vertex = vertices[i];
			List<AssemblyEdge> edges = new ArrayList<AssemblyEdge>();
			for(AssemblyEdge edge:graph.getEdges(vertex)) if(!edge.isSameSequenceEdge()) edges.add(edge);
			if(vertices[i].getSequenceIndex()==debugSeq) System.out.println("AddEdges3. Vertex "+vertices[i]+" edges: "+edges.size());
			if(edges.size()<2) continue;
			Collections.sort(edges,(e1,e2)-> e2.getScore()-e1.getScore());
			int maxScore= edges.get(0).getScore();
			Collections.sort(edges,(e1,e2)-> e2.getOverlap()-e1.getOverlap());
			
			int maxOverlap = edges.get(0).getOverlap();
			if(vertices[i].getSequenceIndex()==debugSeq) System.out.println("AddEdges3. Vertex "+vertices[i]+" max overlap: "+maxOverlap+" max score: "+maxScore);
			for(AssemblyEdge edge:edges) {
				if(edge.isSameSequenceEdge()) continue;
				if(vertices[i].getSequenceIndex()==debugSeq) System.out.println("AddEdges3. Vertex "+vertices[i]+" next edge: "+edge);
				if(edge.getScore()<0.9*maxScore) continue;
				if(edge.getOverlap()<0.9*maxOverlap) break;
				AssemblyVertex v = edge.getConnectingVertex(vertex);
				Integer j = verticesMap.get(v.getUniqueNumber());
				if(j==null || used[j]) continue;
				if(paths.get(j/2).size()>1) continue;
				if(vertices[i].getSequenceIndex()==debugSeq) System.out.println("AddEdges3. Vertex "+vertices[i]+" adding vertex: "+v);
				Set<Integer> selectedSequences = vertexSequences.computeIfAbsent(i, (k)->new HashSet<Integer>());
				selectedSequences.add(v.getSequenceIndex());
			}
			//if(vertices[i].getUniqueNumber()==0) System.out.println("Vertex "+vertices[i]+" sequences: "+vertexSequences.get(0));
		}
		for(int i:vertexSequences.keySet()) {
			Set<Integer> sequencesI = vertexSequences.get(i);
			if(vertices[i].getUniqueNumber()==30 || vertices[i].getUniqueNumber()==-1032) System.out.println("AddEdges3. Vertex "+vertices[i]+" sequences: "+sequencesI+" used: "+used[i]);
			if(used[i]) continue;
			if(sequencesI.size()>5) continue;
			for(int j:vertexSequences.keySet()) {
				if(i>=j) continue;
				if(used[j]) continue;
				if(clusters[i]==clusters[j]) continue;
				Set<Integer> sequencesJ = vertexSequences.get(j);
				if(!sequencesI.equals(sequencesJ)) continue;
				List<Integer> groupIndexes = new ArrayList<Integer>();
				Set<Integer> groupClusters = new HashSet<Integer>();
				groupClusters.add(clusters[j]);
				boolean groupFree=true;
				for(int k:sequencesI) {
					int iv1 = verticesMap.get(graph.getVertex(k, true).getUniqueNumber());
					groupFree = groupFree && !used[iv1] && clusters[iv1]!=clusters[i];
					groupIndexes.add(iv1);
					groupClusters.add(clusters[iv1]);
					int iv2 = verticesMap.get(graph.getVertex(k, false).getUniqueNumber());
					groupFree = groupFree && !used[iv2] && clusters[iv2]!=clusters[i];
					groupIndexes.add(iv2);
					groupClusters.add(clusters[iv2]);
				}
				if(!groupFree) continue;
				List<AssemblyEdge> path = findPathExtensiveSearch(graph,vertices[i],vertices[j],sequencesI);
				if(path.size()==0) continue;
				for(AssemblyEdge edge:path) {
					if(edge.isSameSequenceEdge()) continue;
					pathEdges.add(edge);
				}
				used[i] = used[j] = true;
				for(int k:groupIndexes) used[k] = true;
				int ci = clusters[i];
				int cj = clusters[j];
				for(int c=0;c<n;c++) {
					if(clusters[c]==cj || groupClusters.contains(clusters[c])) clusters[c] = ci;
				}
				System.out.println("Added semisafe edges connecting "+vertices[i]+" and "+vertices[j]+" sequences in the middle: "+sequencesI);
				break;
			}
		}
	}
	
	private List<AssemblyEdge> findPathExtensiveSearch(AssemblyGraph graph, AssemblyVertex v1 ,AssemblyVertex v2, Set<Integer> sequencesPath) {
		LinkedList<List<Integer>> agenda = new LinkedList<List<Integer>>();
		List<Integer> origin = new ArrayList<Integer>();
		origin.add(v1.getUniqueNumber());
		agenda.add(origin);
		while(agenda.size()>0) {
			List<Integer> path = agenda.removeFirst();
			int n =path.size();
			int last = path.get(n-1);
			AssemblyVertex lastV = graph.getVertexByUniqueId(last);
			if(n==2*(sequencesPath.size()+1)) {
				if(last!=v2.getUniqueNumber()) continue;
				//Build list of edges and return
				List<AssemblyEdge> answer = new ArrayList<AssemblyEdge>();
				for(int i=0;i<n-1;i++) {
					answer.add(graph.getEdge(graph.getVertexByUniqueId(path.get(i)),graph.getVertexByUniqueId(path.get(i+1))));
				}
				return answer;
			} else if (last==v2.getUniqueNumber()) continue;
			
			List<AssemblyEdge> edges = graph.getEdges(lastV);
			for(AssemblyEdge edge:edges) {
				AssemblyVertex v = edge.getConnectingVertex(lastV);
				if(path.contains(v.getUniqueNumber())) continue;
				if(v!=v2 && !sequencesPath.contains(v.getSequenceIndex())) continue;
				List<Integer> newPath = new ArrayList<Integer>(n+2);
				newPath.addAll(path);
				newPath.add(v.getUniqueNumber());
				if(v!=v2) newPath.add(graph.getSameSequenceEdge(v).getConnectingVertex(v).getUniqueNumber());
				agenda.add(newPath);
			}
			
		}
		return new ArrayList<AssemblyEdge>();
	}
	private void addEdges2(AssemblyGraph graph, List<LinkedList<AssemblyEdge>> paths, List<AssemblyEdge> pathEdges) {
		AssemblyVertex [] vertices = extractEndVertices(paths);
		int debugSeq = -1;
		int n = vertices.length;
		int [] clusters = new int[n];
		boolean [] used = new boolean[n];
		Arrays.fill(used, false);
		for(int i=0;i<n;i++) {
			clusters[i] = i/2;
		}
		Map<Integer,Integer> verticesMap = new HashMap<Integer, Integer>(n);
		for(int i=0;i<n;i++) verticesMap.put(vertices[i].getUniqueNumber(), i);
		// Try to expand from vertex i
		for(int i=0;i<vertices.length;i++) {
			if(used[i]) continue;
			AssemblyVertex vertex = vertices[i];
			List<AssemblyEdge> edges = graph.getEdges(vertex);
			AssemblyEdge bestOverlap=null;
			AssemblyEdge bestWCSK=null;
			for(AssemblyEdge edge:edges) {
				if(edge.isSameSequenceEdge()) continue;
				Integer j = verticesMap.get(edge.getConnectingVertex(vertex).getUniqueNumber());
				if(j==null || used[j]) continue;
				if(bestOverlap==null || bestOverlap.getOverlap()<edge.getOverlap()) bestOverlap = edge;
				if(bestWCSK==null || bestWCSK.getWeightedCoverageSharedKmers()<edge.getWeightedCoverageSharedKmers()) bestWCSK = edge;
			}
			if(vertex.getSequenceIndex()==debugSeq) System.out.println("AddEdges2. Trying to connect vertex "+vertex+" edges: "+edges.size()+" max overlap: "+bestOverlap+" max WCSK: "+bestWCSK);
			if(bestOverlap == bestWCSK) continue;
			int diffOverlap = Math.abs(bestOverlap.getOverlap()-bestWCSK.getOverlap());
			int diffWCSK = Math.abs(bestOverlap.getWeightedCoverageSharedKmers()-bestWCSK.getWeightedCoverageSharedKmers());
			if(vertex.getSequenceIndex()==debugSeq) System.out.println("AddEdges2. Trying to connect vertex "+vertex+" diff overlap: "+diffOverlap+" diff wcsk: "+diffWCSK);
			if(diffOverlap>0.02*bestOverlap.getOverlap()) continue;
			if(diffWCSK>0.02*bestWCSK.getWeightedCoverageSharedKmers()) continue;
			
			AssemblyEdge thirdOverlap=null;
			AssemblyEdge thirdWCSK=null;
			for(AssemblyEdge edge:edges) {
				if(edge.isSameSequenceEdge()) continue;
				if(edge == bestOverlap) continue;
				if(edge == bestWCSK) continue;
				Integer j = verticesMap.get(edge.getConnectingVertex(vertex).getUniqueNumber());
				if(j==null || used[j]) continue;
				if(thirdOverlap==null || thirdOverlap.getOverlap()<edge.getOverlap()) thirdOverlap = edge;
				if(thirdWCSK==null || thirdWCSK.getWeightedCoverageSharedKmers()<edge.getWeightedCoverageSharedKmers()) thirdWCSK = edge;
			}
			if(vertex.getSequenceIndex()==debugSeq) System.out.println("AddEdges2. Trying to connect vertex "+vertex+" third overlap: "+thirdOverlap+" third WCSK: "+thirdWCSK);
			if(thirdOverlap!=null && (thirdOverlap.getOverlap()>=bestWCSK.getOverlap() || thirdOverlap.getWeightedCoverageSharedKmers()>=bestOverlap.getWeightedCoverageSharedKmers())) continue;
			if(thirdWCSK!=null && (thirdWCSK.getWeightedCoverageSharedKmers()>=bestOverlap.getWeightedCoverageSharedKmers()|| thirdWCSK.getOverlap()>=bestWCSK.getOverlap())) continue;
			
			AssemblyVertex vBestOv = bestOverlap.getConnectingVertex(vertex);
			AssemblyVertex vBestWCSK = bestWCSK.getConnectingVertex(vertex);
			
			int j = verticesMap.get(vBestOv.getUniqueNumber());
			int k = verticesMap.get(vBestWCSK.getUniqueNumber());
			if(paths.get(j/2).size()>1) continue;
			if(paths.get(k/2).size()>1) continue;
			
			AssemblyVertex vBestOv2 = graph.getSameSequenceEdge(vBestOv).getConnectingVertex(vBestOv);
			AssemblyVertex vBestWCSK2 = graph.getSameSequenceEdge(vBestWCSK).getConnectingVertex(vBestWCSK);
			AssemblyEdge edgeOv2 = graph.getEdge(vBestOv2, vBestWCSK);
			AssemblyEdge edgeWCSK2 = graph.getEdge(vBestWCSK2, vBestOv);
			
			if(graph.isBestEdge(vBestOv, bestOverlap) && edgeOv2!=null && edgeWCSK2==null && edgeOv2.getOverlap()>0.8*bestWCSK.getOverlap()) {
				//Best overlap follows transitivity
				int c1 = clusters[i];
				int c2 = clusters[j];
				int c3 = clusters[k];
				
				if(c1!=c2 && c1!=c3 && c2!=c3) {
					pathEdges.add(bestOverlap);
					pathEdges.add(edgeOv2);
					used[i] = used[j] = used[k] = true;
					int j2 = verticesMap.get(vBestOv2.getUniqueNumber());
					used[j2] = true;
					for(int c=0;c<n;c++) {
						if(clusters[c]==c2 || clusters[c]==c3) clusters[c] = c1;
					}
					System.out.println("Added transitivity edges with best overlap "+bestOverlap+" "+edgeOv2);
					System.out.println("Competing edge: "+bestWCSK);
				}
			} else if (graph.isBestEdge(vBestWCSK, bestWCSK) && edgeOv2==null && edgeWCSK2!=null && edgeWCSK2.getOverlap()>0.8*bestOverlap.getOverlap()) {
				//Best WCSK follows transitivity
				int c1 = clusters[i];
				int c2 = clusters[j];
				int c3 = clusters[k];
				
				if(c1!=c2 && c1!=c3 && c2!=c3) {
					pathEdges.add(bestWCSK);
					pathEdges.add(edgeWCSK2);
					used[i] = used[j] = used[k] = true;
					int k2 = verticesMap.get(vBestWCSK2.getUniqueNumber());
					used[k2] = true;
					for(int c=0;c<n;c++) {
						if(clusters[c]==c2 || clusters[c]==c3) clusters[c] = c1;
					}
					System.out.println("Added transitivity edges with best WCSK: "+bestWCSK+" "+edgeWCSK2);
					System.out.println("Competing edge: "+bestOverlap);
				}
			}
		}
	}

}
class VertexPathLocation {
	private int vertexId;
	private int pathId;
	private int pathLength;
	private int pathPosition;
	public VertexPathLocation(int vertexId, int pathId, int pathLength, int pathPosition) {
		super();
		this.vertexId = vertexId;
		this.pathId = pathId;
		this.pathLength = pathLength;
		this.pathPosition = pathPosition;
	}
	public int getVertexId() {
		return vertexId;
	}
	public int getPathId() {
		return pathId;
	}
	public int getPathLength() {
		return pathLength;
	}
	public int getPathPosition() {
		return pathPosition;
	}
	
}
class PathEndJunctionEdge {
	private int path1EndId;
	private int path2EndId;
	private int totalCost=0;
	private int countVotes=0;
	
	
	
	public PathEndJunctionEdge(int path1EndId, int path2EndId) {
		super();
		this.path1EndId = path1EndId;
		this.path2EndId = path2EndId;
	}
	public void addVote (int cost) {
		totalCost+=cost;
		countVotes++;
	}
	public int getPath1EndId() {
		return path1EndId;
	}
	public int getPath2EndId() {
		return path2EndId;
	}
	public int getCost() {
		if(countVotes==0) return Integer.MAX_VALUE;
		return totalCost/countVotes;
	}
	public String toString() {
		return "P1: "+path1EndId+" P2: "+path2EndId+" Total cost: "+totalCost+" Votes: "+countVotes;
	}
	
	
}
