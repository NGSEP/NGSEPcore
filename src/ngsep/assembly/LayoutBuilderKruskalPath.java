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

	private int minPathLength = LayoutBuilder.DEF_MIN_PATH_LENGTH;
	
	private boolean runImprovementAlgorithms = true;
	
	
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
	
	public boolean isRunImprovementAlgorithms() {
		return runImprovementAlgorithms;
	}
	public void setRunImprovementAlgorithms(boolean runImprovementAlgorithms) {
		this.runImprovementAlgorithms = runImprovementAlgorithms;
	}
	@Override
	public void findPaths(AssemblyGraph graph) {
		List<AssemblyEdge> pathEdges = graph.selectSafeEdges();
		log.info("Number of safe edges: "+pathEdges.size());
		
		List<AssemblyPath> safePaths = graph.buildPaths(pathEdges);
		//Map<Integer,Integer> vertexPathIds = calculateVertexPathIdsMap(safePaths);
		log.info("Number of paths safe edges: "+safePaths.size());
		
		addConnectingEdges(graph, safePaths, pathEdges);
		List<AssemblyPath> rawPathsCostsAlgorithm = graph.buildPaths(pathEdges);
		log.info("Paths costs algorithm: "+rawPathsCostsAlgorithm.size());
		List<AssemblyPath> paths = new ArrayList<>();
		for(AssemblyPath path:rawPathsCostsAlgorithm) {
			if(path.getPathLength()>1) paths.add(path);
			//else path.print(System.out);
		}
		log.info("Paths removing unconnected reads: "+paths.size());
		if(runImprovementAlgorithms) {
			Distribution [] distsEdges = calculateDistributions(pathEdges);
			log.info("Average path edge cost: "+distsEdges[0].getAverage());
			distsEdges[0].printDistributionInt(System.out);
			paths = collectAlternativeSmallPaths(graph, paths);
			log.info("Paths after collecting small paths: "+paths.size());
			paths = mergeClosePaths(graph, paths, distsEdges);
			log.info("Paths after first round of merging: "+paths.size());
			expandPathsWithEmbedded(graph, paths, distsEdges);
			log.info("Paths after expanding with embedded: "+paths.size());
			paths = mergeClosePaths(graph, paths, distsEdges);
			log.info("Paths after second round of merging: "+paths.size());
			paths = collectAlternativeSmallPaths(graph, paths);
			log.info("Paths after second round collecting small paths: "+paths.size());
			//paths = mergeClosePaths(graph, paths, distsEdges);
			//log.info("Paths after third round of merging: "+paths.size());
		}
		
		for(AssemblyPath path:paths) {
			if(path.getPathLength()>=minPathLength) graph.addPath(path);
			//else path.print(System.out);
		}
		log.info("Final number of paths: "+graph.getPaths().size());
		System.out.println("Estimated N statistics");
		long [] stats = graph.estimateNStatisticsFromPaths();
		if(stats!=null) NStatisticsCalculator.printNStatistics(stats, System.out);
	}
 
	private AssemblyVertex [] extractEndVertices (List<AssemblyPath> paths) {
		int p = paths.size();
		AssemblyVertex [] vertices = new AssemblyVertex[2*p];
		int v=0;
		for(int i=0;i<p;i++) {
			AssemblyPath path = paths.get(i);
			vertices[v] = path.getVertexLeft();
			v++;
			vertices[v] = path.getVertexRight();
			v++;
		}
		return vertices;
	}
	private void addConnectingEdges(AssemblyGraph graph, List<AssemblyPath> paths, List<AssemblyEdge> pathEdges) {
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
		double limitIKBP = distIKBP.getMean()+15*Math.sqrt(distIKBP.getVariance());
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
			//if(v1.getUniqueNumber()==-97856 || v2.getUniqueNumber()==69473) System.out.println("SelectingEdges. next edge: "+nextEdge+" used: "+used[posV1]+" "+used[posV2]+" clusters: "+clusters[posV1]+" "+clusters[posV2]);
			
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
	private Distribution[] calculateDistributions(List<AssemblyEdge> pathEdges) {
		Distribution costsDistribution = new Distribution(0, 100000, 1000);
		Distribution ikbpDistribution = new Distribution(0, 50, 0.25);
		for(AssemblyEdge edge:pathEdges) {
			if(edge.isSameSequenceEdge()) continue;
			costsDistribution.processDatapoint(edge.getCost());
			ikbpDistribution.processDatapoint(edge.getIndelsPerKbp());
		}
		//log.info("Average cost path edges: "+costsDistribution.getAverage()+" STDEV: "+Math.sqrt(costsDistribution.getVariance()));
		//costsDistribution.printDistribution(System.out);
		Distribution[] answer = {costsDistribution,ikbpDistribution}; 
		return answer;
	}
	private List<AssemblyPath> collectAlternativeSmallPaths(AssemblyGraph graph, List<AssemblyPath> paths) {
		Map<Integer,VertexPathLocation> vertexPositions = getPathPositionsMap(paths);
		Set<Integer> indexesToRemove = new HashSet<>();
		for(int i=0;i<paths.size();i++) {
			AssemblyPath path = paths.get(i);
			//log.info("CollectSmallPaths. Next path: "+(i+1)+" length " + path.getPathLength());
			if(path.getPathLength()>20) continue;
			AssemblyVertex leftVertex = path.getVertexLeft();
			AssemblyVertex rightVertex = path.getVertexRight();
			AssemblyEdge leftEdge = graph.getEdgeMinCost(leftVertex);
			AssemblyEdge rightEdge = graph.getEdgeMinCost(rightVertex);
			if(leftEdge==null || rightEdge == null) continue;
			AssemblyVertex leftConnecting = leftEdge.getConnectingVertex(leftVertex);
			AssemblyVertex rightConnecting = rightEdge.getConnectingVertex(rightVertex);
			if(leftConnecting==null || rightConnecting == null) continue;
			VertexPathLocation leftLocation = vertexPositions.get(leftConnecting.getUniqueNumber());
			VertexPathLocation rightLocation = vertexPositions.get(rightConnecting.getUniqueNumber());
			//log.info("CollectSmallPaths. Next path: "+(i+1)+" length " + path.getPathLength()+" end vertices: "+leftVertex+" "+rightVertex+" connecting "+leftConnecting+" "+rightConnecting+" "+ leftLocation+" "+rightLocation);
			if(leftLocation==null || rightLocation == null) continue;
			if(leftLocation.getPath()==path) continue;
			if(leftLocation.getPath()!=rightLocation.getPath()) continue;
			AssemblyPath hostPath = leftLocation.getPath();
			
			
			if(0.1*hostPath.getPathLength()<path.getPathLength()) continue;
			if(Math.abs(leftLocation.getPathPosition()-rightLocation.getPathPosition())>1.5*path.getPathLength()) continue;
			//log.info("CollectSmallPaths. Integration of path: "+(i+1)+" into path with length: "+hostPath.getPathLength()+" conecting pos: "+leftLocation.getPathPosition()+" "+rightLocation.getPathPosition());
			hostPath.addAlternativeSmallPath(path);
			indexesToRemove.add(i);
		}
		log.info("Internal path ids: "+indexesToRemove.size());
		List<AssemblyPath> answer = new ArrayList<AssemblyPath>();
		for(int i=0;i<paths.size();i++) {
			if(!indexesToRemove.contains(i)) answer.add(paths.get(i));
		}
		return answer;
	}
	private Map<Integer,VertexPathLocation> getPathPositionsMap (List<AssemblyPath> paths) {
		Map<Integer,VertexPathLocation> vertexPositions = new HashMap<Integer,VertexPathLocation>();
		for(int i=0;i<paths.size();i++) {
			AssemblyPath path = paths.get(i);
			List<AssemblyEdge> edges = path.getEdges();
			int n = edges.size();
			AssemblyVertex vertexLeft = path.getVertexLeft();
			vertexPositions.put(vertexLeft.getUniqueNumber(), new VertexPathLocation(vertexLeft, path, n, 0));
			AssemblyVertex lastVertex = vertexLeft;
			int k = 1;
			for(AssemblyEdge edge:edges) {
				AssemblyVertex nextVertex = edge.getConnectingVertex(lastVertex);
				vertexPositions.put(nextVertex.getUniqueNumber(), new VertexPathLocation(nextVertex, path, n, k));
				k++;
				lastVertex = nextVertex;
			}
		}
		return vertexPositions;
	}
	private boolean expandPathsWithEmbedded(AssemblyGraph graph, List<AssemblyPath> paths, Distribution [] dists) {
		double averageCost = dists[0].getAverage();
		double averageIKBP = dists[1].getAverage();
		log.info("Expanding paths with embedded. Average cost: "+averageCost+" average IKBP: "+averageIKBP);
		Set<Integer> usedEmbedded = new HashSet<Integer>();
		for(AssemblyPath path:paths) {
			
			//System.out.println("Expanding paths with embedded. Next start: "+endVertex);
			while(true) {
				AssemblyVertex endVertex = path.getVertexLeft();
				AssemblyEdge edge = graph.getEdgeMinCost(endVertex);
				if(edge==null) break;
				if(edge.getCost()>2*averageCost || edge.getIndelsPerKbp()>averageIKBP) break;
				AssemblyVertex next = edge.getConnectingVertex(endVertex);
				if(!usedEmbedded.contains(next.getSequenceIndex()) && graph.isEmbedded(next.getSequenceIndex())) {
					//log.info("Expanding path starting with: "+path.getVertexLeft()+" with embedded. Edge: "+edge);
					path.connectEdgeLeft(graph, edge);
					endVertex = path.getVertexLeft();
					usedEmbedded.add(endVertex.getSequenceIndex());
				} else break;
			}
			//System.out.println("Expanding paths with embedded. Next end: "+endVertex);
			while(true) {
				AssemblyVertex endVertex = path.getVertexRight();
				AssemblyEdge edge = graph.getEdgeMinCost(endVertex);
				//if(endVertex.getUniqueNumber()==24728) System.out.println("Expanding paths with embedded. Best overlap: "+edge);
				if(edge==null) break;
				if(edge.getCost()>2*averageCost || edge.getIndelsPerKbp()>averageIKBP) break;
				AssemblyVertex next = edge.getConnectingVertex(endVertex);
				if(!usedEmbedded.contains(next.getSequenceIndex()) && graph.isEmbedded(next.getSequenceIndex())) {
					//log.info("Expanding path ending with: "+path.getVertexRight()+" with embedded. Edge: "+edge);
					path.connectEdgeRight(graph, edge);
					endVertex = path.getVertexRight();
					usedEmbedded.add(endVertex.getSequenceIndex());
				} else break;
			}
		}
		return usedEmbedded.size()>0;
	}
	private List<AssemblyPath> mergeClosePaths (AssemblyGraph graph, List<AssemblyPath> paths, Distribution [] dists) {
		Map<String,PathEndJunctionEdge> pathEndEdges = buildPathJunctionEdgesByCloseEdges(graph, paths, dists);
		log.info("Found "+pathEndEdges.size()+" candidate edges connecting vertices");
		//for(PathEndJunctionEdge edge:pathEndEdges.values()) {
			//if(edge.getPath1EndId()<200 && edge.getPath2EndId()<200) System.out.println(edge);
		//}
		Map<Integer,List<Integer>> mergedPathIds = findPathsToMerge(pathEndEdges.values(),paths.size());
		log.info("Merged path ids");
		//for(List<Integer> nextPath:mergedPathIds.values()) {
			//if(nextPath.size()>2) System.out.println(nextPath);
		//}
		//return paths;
		return mergePaths(graph, paths, mergedPathIds, pathEndEdges);
	}
	
	private Map<String,PathEndJunctionEdge> buildPathJunctionEdgesByCloseEdges(AssemblyGraph graph, List<AssemblyPath> paths, Distribution[] dists) {
		log.info("Merging paths from "+paths.size()+" paths");
		boolean debug = false;
		Map<Integer,VertexPathLocation> vertexPositions = getPathPositionsMap(paths);
		//Find edges connecting paths
		Map<String,PathEndJunctionEdge> pathEndEdges = new HashMap<String, PathEndJunctionEdge>();

		Distribution costsDistribution = dists[0];
		double limitCost = costsDistribution.getAverage()+4*Math.sqrt(costsDistribution.getVariance());
		log.info("Cost limit to identify possible merging edges: "+limitCost);
		for(int i=0;i<paths.size();i++) {
			AssemblyPath path = paths.get(i);
			int pathId = i+1;
			path.setPathId(pathId);
			if(debug) System.out.println("next input path. id: "+pathId+" left: "+path.getVertexLeft()+" right: "+path.getVertexRight()+" edges: "+path.getEdges().size());
		}

    	//Build connections graph
    	for (Map.Entry<Integer, VertexPathLocation> entry:vertexPositions.entrySet()) {
    		VertexPathLocation loc1 = entry.getValue();
    		Integer pathEndV1 = loc1.getPathEnd();
    		if(pathEndV1==null) continue;
    		if(debug && pathEndV1==-2) System.out.println("Vertex: "+loc1.getVertex()+" path: "+loc1.getPath().getPathId()+" position: "+loc1.getPathPosition()+" end: "+pathEndV1);
    		
    		AssemblyVertex v1 = loc1.getVertex();
    		List<AssemblyEdge> edges = graph.getEdges(v1);
    		for(AssemblyEdge edge:edges) {
    			if(edge.isSameSequenceEdge()) continue;
    			//TODO. Improve cost
    			if(edge.getCost()>limitCost) continue;
    			AssemblyVertex v2 = edge.getConnectingVertex(v1);
    			VertexPathLocation loc2 = vertexPositions.get(v2.getUniqueNumber());
    			if(loc2==null) continue;
    			if(loc1.getPath()==loc2.getPath()) continue;
    			Integer pathEndV2 = loc2.getPathEnd();
    			if(debug && pathEndV1==-2) System.out.println("Vertex 2: "+v2+" cost: "+edge.getCost()+" pathEnd2: "+pathEndV2);
    			if(pathEndV2 == null) continue;
    			addVote(pathEndV1, pathEndV2, edge, pathEndEdges);
    		}
    		List<AssemblyEmbedded> embeddedList = graph.getEmbeddedByHostId(v1.getSequenceIndex());
    		for(AssemblyEmbedded embedded:embeddedList) {
    			if(embedded.getCost()>limitCost) continue;
    			AssemblyVertex v2 = graph.getVertex(embedded.getSequenceId(), true);
    			VertexPathLocation loc2 = vertexPositions.get(v2.getUniqueNumber());
    			if(loc2==null) continue;
    			if(loc1.getPath()==loc2.getPath()) continue;
    			Integer pathEndV2 = loc2.getPathEnd();
    			if(pathEndV2 == null) continue;
    			addVote(pathEndV1, pathEndV2, embedded, pathEndEdges);
    		}
    	}
    	return pathEndEdges;
	}
	
	private void addVote (int pathEndV1, int pathEndV2, AssemblySequencesRelationship rel,  Map<String,PathEndJunctionEdge> pathEndEdges) {
		
		int minId = Math.min(pathEndV1, pathEndV2);
		int maxId = Math.max(pathEndV1, pathEndV2);
		String key = "F"+minId+"L"+maxId;
		PathEndJunctionEdge pathEndEdge = pathEndEdges.computeIfAbsent(key, v->new PathEndJunctionEdge(minId, maxId));
		pathEndEdge.addVote(rel);
	}
	
	private Map<Integer,List<Integer>> findPathsToMerge(Collection<PathEndJunctionEdge> pathEndEdges, int n) {
		boolean debug = false;
		List<PathEndJunctionEdge> pathEndEdgesList = new ArrayList<PathEndJunctionEdge>(pathEndEdges);
		Collections.sort(pathEndEdgesList,(e1,e2)->e1.getCost()-e2.getCost());
		Map<Integer,Integer> pathGroups = new HashMap<Integer, Integer>();
		
		Map<Integer,List<Integer>> answer = new TreeMap<Integer,List<Integer>>();
		Set<Integer> usedPathEnds = new HashSet<Integer>();
		for(int i=0;i<n;i++) {
			int pathId = i+1;
			pathGroups.put(pathId, pathId);
			pathGroups.put(-pathId, pathId);
			LinkedList<Integer> nextPath = new LinkedList<Integer>();
			nextPath.add(pathId);
			nextPath.add(-pathId);
			answer.put(pathId,nextPath);
		}
		
		for(PathEndJunctionEdge pathEdge:pathEndEdgesList) {
			if(debug) System.out.println("Path edge: "+pathEdge.getPath1EndId()+" "+pathEdge.getPath2EndId()+" cost: "+pathEdge.getCost());
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
	
	private List<AssemblyPath> mergePaths(AssemblyGraph graph, List<AssemblyPath> paths, Map<Integer, List<Integer>> mergedPathIds, Map<String,PathEndJunctionEdge> pathEndEdges) {
		boolean debug = false;
		List<AssemblyPath> answer = new ArrayList<>();
		for(List<Integer> idsPath:mergedPathIds.values()) {
			AssemblyPath nextPath = null;
			int lastEndId = 0;
			PathEndJunctionEdge je = null;
			for(int pathEndId:idsPath) {
				if(pathEndId!=-lastEndId) {
					int minId = Math.min(lastEndId, pathEndId);
					int maxId = Math.max(lastEndId, pathEndId);
					String key = "F"+minId+"L"+maxId;
					je = pathEndEdges.get(key);
					lastEndId = pathEndId;
					continue;
				}
				if(pathEndId > 0) {
					AssemblyPath nextInputPath = paths.get(pathEndId-1);
					if(nextPath==null) {
						nextInputPath.reverse();
						nextPath = nextInputPath;
					} else {
						List<AssemblyPath> leftoverpaths = nextPath.connectPathRight(graph, nextInputPath, true, je.getRelationships()); 
						if(leftoverpaths==null) {
							answer.add(nextPath);
							nextPath = nextInputPath;
						} else {
							if(debug) System.out.println("merged path: "+nextPath.getPathId() +" with "+nextInputPath.getPathId()+". leftover paths: "+leftoverpaths.size());
							answer.addAll(leftoverpaths);
						}
					}
				} else {
					AssemblyPath nextInputPath = paths.get(lastEndId-1);
					if(nextPath==null) {
						nextPath = nextInputPath;
					} else {
						List<AssemblyPath> leftoverpaths = nextPath.connectPathRight(graph, nextInputPath, false, je.getRelationships());
						if(leftoverpaths==null) {
							answer.add(nextPath);
							nextPath = nextInputPath;
						} else {
							if(debug) System.out.println("merged path: "+nextPath.getPathId() +" with "+nextInputPath.getPathId()+" leftover paths: "+leftoverpaths.size());
							answer.addAll(leftoverpaths);
						}
						
					}
				}
				lastEndId = pathEndId;
			}
			//log.info("Next final path size: "+nextPath.size());
			if(nextPath!=null) answer.add(nextPath);
		}
		return answer;
	}
	
	

}
class VertexPathLocation {
	private AssemblyVertex vertex;
	private AssemblyPath path;
	private int pathLength;
	private int pathPosition;
	public VertexPathLocation(AssemblyVertex vertex, AssemblyPath path, int pathLength, int pathPosition) {
		super();
		this.vertex = vertex;
		this.path = path;
		this.pathLength = pathLength;
		this.pathPosition = pathPosition;
	}
	
	public Integer getPathEnd() {
		int distanceEnd = pathLength-pathPosition;
		boolean outwardsLeft=pathPosition%2==0;
		
		if(pathPosition<20 && outwardsLeft ) return path.getPathId();
		if(distanceEnd<20 && !outwardsLeft ) return -path.getPathId();
		return null;
	}

	public AssemblyVertex getVertex() {
		return vertex;
	}

	public AssemblyPath getPath() {
		return path;
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
	private List<AssemblySequencesRelationship> relationships = new ArrayList<>();
	
	
	public PathEndJunctionEdge(int path1EndId, int path2EndId) {
		super();
		this.path1EndId = path1EndId;
		this.path2EndId = path2EndId;
	}
	public void addVote (AssemblySequencesRelationship rel) {
		relationships.add(rel);
		totalCost+=rel.getCost();
		countVotes++;
	}
	public int getPath1EndId() {
		return path1EndId;
	}
	public int getPath2EndId() {
		return path2EndId;
	}
	
	public int getTotalCost() {
		return totalCost;
	}
	public int getCountVotes() {
		return countVotes;
	}
	public int getCost() {
		if(countVotes==0) return Integer.MAX_VALUE;
		return totalCost/countVotes;
	}
	
	public List<AssemblySequencesRelationship> getRelationships() {
		return relationships;
	}
	public String toString() {
		return "P1: "+path1EndId+" P2: "+path2EndId+" Total cost: "+totalCost+" Votes: "+countVotes+" rels: "+relationships;
	}
	
	
}
