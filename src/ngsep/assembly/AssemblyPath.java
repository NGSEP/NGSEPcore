package ngsep.assembly;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class AssemblyPath {
	private LinkedList<AssemblyEdge> edges = new LinkedList<AssemblyEdge>();
	private AssemblyVertex vertexLeft;
	private AssemblyVertex vertexRight;
	private int pathId;
	private String sequenceName;
	private List<AssemblyPath> alternativeSmallPaths=new ArrayList<>();
	private String consensus;
	private Map<Integer, Integer> pathVerticesConsensusEnds;
	
	public AssemblyPath(AssemblyEdge edge) {
		edges.add(edge);
		vertexLeft = edge.getVertex1();
		vertexRight = edge.getVertex2();
	}
	public void connectEdgeLeft (AssemblyGraph graph, AssemblyEdge edge) {
		AssemblyVertex newSecond = edge.getConnectingVertex(vertexLeft);
		if(newSecond==null && edges.size()==1) {
			newSecond = edge.getConnectingVertex(vertexRight);
			if(newSecond!=null) reverse();
		}
		if(newSecond==null) throw new RuntimeException("Invalid edge to add as left edge. Current left vertex: "+vertexLeft+" incorrect edge. "+edge);
		AssemblyEdge newLeftEdge = graph.getSameSequenceEdge(newSecond);
		if(newLeftEdge==null) throw new RuntimeException("Invalid edge to add as left edge. Following same sequence edge not found. Current left vertex: "+vertexLeft+" edge: "+edge);
		AssemblyVertex newLeft = newLeftEdge.getConnectingVertex(newSecond);
		if(newLeft==null) throw new RuntimeException("Error adding edge. Out vertex not found for self sequence edge. Second vertex: "+newSecond+" self sequence edge: "+newLeftEdge );
		edges.add(0, edge);
		edges.add(0,newLeftEdge);
		vertexLeft = newLeft;
	}
	public void connectEdgeRight (AssemblyGraph graph, AssemblyEdge edge) {
		AssemblyVertex newSecond = edge.getConnectingVertex(vertexRight);
		if(newSecond==null && edges.size()==1) {
			newSecond = edge.getConnectingVertex(vertexLeft);
			if(newSecond!=null) reverse();
		}
		if(newSecond==null) throw new RuntimeException("Invalid edge to add as right edge. Current vertex right: "+vertexRight+" incorrect edge. "+edge);
		AssemblyEdge newRightEdge = graph.getSameSequenceEdge(newSecond);
		if(newRightEdge==null) throw new RuntimeException("Invalid edge to add as right edge. Following same sequence edge not found. Current vertex right: "+vertexRight+" edge: "+edge);
		AssemblyVertex newRight = newRightEdge.getConnectingVertex(newSecond);
		if(newRight==null) throw new RuntimeException("Error adding edge. Out vertex not found for self sequence edge. Second right vertex: "+newSecond+" self sequence edge: "+newRightEdge );
		edges.add(edge);
		edges.add(newRightEdge);
		vertexRight = newRight;
	}
	public int getPathLength() {
		return edges.size();
	}
	public List<AssemblyEdge> getEdges() {
		return edges;
	}
	public AssemblyVertex getVertexLeft() {
		return vertexLeft;
	}
	public AssemblyVertex getVertexRight() {
		return vertexRight;
	}
	public int getPathId() {
		return pathId;
	}
	public void setPathId(int pathId) {
		this.pathId = pathId;
	}
	public String getSequenceName() {
		return sequenceName;
	}
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
	public List<AssemblyVertex> extractVertices() {
		List<AssemblyVertex> answer = new ArrayList<AssemblyVertex>(edges.size()+1);
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) {
				answer.add(edge.getVertex1());
				answer.add(edge.getVertex2());
			}
		}
		return answer;
	}
	public boolean connectPathRight(AssemblyGraph graph, AssemblyPath path2, boolean reverse) {
		int n1 = edges.size();
		AssemblyEdge lastEdge1 = edges.get(n1-1);
		AssemblyVertex lastVertex1 = getVertexRight();
		AssemblyVertex secondLastVertex1 = null;
		int lastOverlap = 0;
		int lastCost = 0;
		if(n1>1) {
			AssemblyEdge secondLastEdge1 = edges.get(n1-2);
			AssemblyVertex connecting1 = lastEdge1.getSharedVertex(secondLastEdge1);
			secondLastVertex1 = secondLastEdge1.getConnectingVertex(connecting1);
			lastOverlap = secondLastEdge1.getOverlap();
			lastCost = secondLastEdge1.getCost();
		}
		
		List<AssemblyEdge> path2Edges = new ArrayList<AssemblyEdge>(path2.edges);
		if(reverse) Collections.reverse(path2Edges);
		int n2 = path2Edges.size();
		AssemblyEdge firstEdge2 = path2Edges.get(0);
		AssemblyVertex firstVertex2 = (reverse?path2.getVertexRight():path2.getVertexLeft());
		AssemblyVertex secondVertex2 = null;
		int nextOverlap = 0;
		int nextCost = 0;
		if(n2>1) {
			AssemblyEdge secondEdge2 = path2.edges.get(1);
			AssemblyVertex connecting2 = firstEdge2.getSharedVertex(secondEdge2);
			secondVertex2 = secondEdge2.getConnectingVertex(connecting2);
			nextOverlap = secondEdge2.getOverlap();
			nextCost = secondEdge2.getCost();
		}
		
		
		AssemblyEdge c1 = graph.getEdge(lastVertex1, firstVertex2);
		//System.err.println("Trying to merge paths with ends: "+lastVertex1+" "+firstVertex2+" Merging direct edge: "+c1);
		if(c1!=null && (n1==1 || c1.getCost()<2*lastCost) && c1.getCost()<2*nextCost && c1.getOverlap()>0.5*lastOverlap && (n2==1 || c1.getOverlap()>0.5*nextOverlap)) {
			System.err.println("Merging paths with ends: "+lastVertex1+" "+firstVertex2+" Merging direct edge: "+c1);
			connectEdgeRight(graph, c1);
			for(AssemblyEdge edge:path2Edges) {
				if(!edge.isSameSequenceEdge()) connectEdgeRight(graph, edge);
			}
			this.alternativeSmallPaths.addAll(path2.alternativeSmallPaths);
			return true;
		}
		if(secondLastVertex1 == null || secondVertex2==null) return false; 
		AssemblyEdge c2 = graph.getEdge(lastVertex1, secondVertex2);
		AssemblyEdge c3 = graph.getEdge(secondLastVertex1, firstVertex2);
		//System.err.println("Edges to merge paths: "+c2+" "+c3+" overlaps: "+lastOverlap+" "+nextOverlap+" costs: "+lastCost+" "+nextCost);
		boolean validC2 = c2!=null && c2.getOverlap()>0.8*nextOverlap && c2.getCost()<1.5*lastCost && c2.getCost()<1.5*nextCost && (c3==null || c2.getCost()<2*c3.getCost());
		if(!validC2) c2= null;
		boolean validC3 = c3!=null && c3.getOverlap()>0.8*lastOverlap && c3.getCost()<1.5*lastCost && c3.getCost()<1.5*nextCost && (c2==null || c3.getCost()<2*c2.getCost());
		if(!validC3) c3= null;
		
		if(c2!=null && (c3==null || c2.getOverlap()>c3.getOverlap())) {
			System.err.println("Merging path with ends: "+secondLastVertex1+" "+lastVertex1+" with "+firstVertex2+" "+ secondVertex2+" Merging edge: "+c2);
			connectEdgeRight(graph, c2);
			for(int i=3;i<n2;i++) {
				AssemblyEdge edge = path2Edges.get(i);
				if(!edge.isSameSequenceEdge()) connectEdgeRight(graph, edge);
			}
			this.alternativeSmallPaths.addAll(path2.alternativeSmallPaths);
			return true;
		}
		if(c3!=null && (c2==null || c3.getOverlap()>c2.getOverlap())) {
			System.err.println("Merging path with ends: "+secondLastVertex1+" "+lastVertex1+" with "+firstVertex2+" "+ secondVertex2+" Merging edge: "+c3);
			edges.removeLast();
			edges.removeLast();
			vertexRight = secondLastVertex1;
			connectEdgeRight(graph, c3);
			for(int i=1;i<n2;i++) {
				AssemblyEdge edge = path2Edges.get(i);
				if(!edge.isSameSequenceEdge()) connectEdgeRight(graph, edge);
			}
			this.alternativeSmallPaths.addAll(path2.alternativeSmallPaths);
			return true;
		}
		return false;
		
	}
	
	public List<AssemblyPath> connectPathRight(AssemblyGraph graph, AssemblyPath path2, boolean reverse, List<AssemblySequencesRelationship> junctionRels) {
		List<AssemblyEdge> path1Edges = new ArrayList<AssemblyEdge>(edges);
		int n1 = path1Edges.size();
		List<AssemblyEdge> path2Edges = new ArrayList<AssemblyEdge>(path2.edges);
		int n2 = path2Edges.size();
		AssemblyVertex nextVertex2 = path2.vertexLeft;
		if(reverse) {
			nextVertex2 = path2.vertexRight;
			Collections.reverse(path2Edges);
		}
		Map<Integer,Integer> lastVerticesLocations1 = new HashMap<>();
		AssemblyVertex nextVertex1 = vertexRight;
		for(int i=n1;i>0 && i>=n1-20;i--) {
			if(i%2==1) lastVerticesLocations1.put(nextVertex1.getUniqueNumber(), i);
			AssemblyEdge edge = path1Edges.get(i-1);
			nextVertex1 = edge.getConnectingVertex(nextVertex1);
		}
		Map<Integer,Integer> lastVerticesLocations2 = new HashMap<>();
		
		for(int i=0;i<n2 && i<=20;i++) {
			if(i%2==0) lastVerticesLocations2.put(nextVertex2.getUniqueNumber(), i);
			AssemblyEdge edge = path2Edges.get(i);
			nextVertex2 = edge.getConnectingVertex(nextVertex2);
		}
		
		Collections.sort(junctionRels,(r1,r2)->r1.getCost()-r2.getCost());
		//System.out.println("Trying to merge path: "+pathId+" with "+path2.getPathId()+" reverse: "+reverse+" rels: "+junctionRels.size());
		//for(AssemblySequencesRelationship rel: junctionRels) System.out.println(rel);
		//Try relationship with lowest cost
		AssemblySequencesRelationship rel = junctionRels.get(0);
		if(rel instanceof AssemblyEdge) {
			AssemblyEdge edge = (AssemblyEdge)rel;
			AssemblyVertex vertexRightAfterRemove = edge.getVertex1(); 
			Integer loc1 = lastVerticesLocations1.get(edge.getVertex1().getUniqueNumber());
			Integer loc2 = lastVerticesLocations2.get(edge.getVertex2().getUniqueNumber());
			if (loc1==null || loc2==null) {
				vertexRightAfterRemove = edge.getVertex2();
				loc1 = lastVerticesLocations1.get(edge.getVertex2().getUniqueNumber());
				loc2 = lastVerticesLocations2.get(edge.getVertex1().getUniqueNumber());
			}
			//System.out.println("Loc1: "+loc1+" loc2: "+loc2+" min cost edge: "+edge);
			if (loc1==null || loc2==null) return null;
			int newPathLength = loc1+(n2-loc2);
			//The new path length should be larger than the original path lengths and the sum of the leftovers
			if(newPathLength<n1 || newPathLength<n2 || newPathLength<=(n1-loc1)+loc2) return null;
			List<AssemblyPath> leftoverPaths = new ArrayList<>();
			AssemblyPath leftover1 = null;
			while(edges.size()>loc1) {
				AssemblyEdge edgeRemoved = edges.removeLast();
				if(leftover1==null) leftover1 = new AssemblyPath(edgeRemoved);
				else if(!edgeRemoved.isSameSequenceEdge()) leftover1.connectEdgeLeft(graph, edgeRemoved);
			}
			if(leftover1!=null && leftover1.getPathLength()>1) leftoverPaths.add(leftover1);
			AssemblyPath leftover2 = null;
			for(int i=0;i<loc2;i++) {
				AssemblyEdge edge2 = path2Edges.get(i);
				if(leftover2==null) leftover2 = new AssemblyPath(edge2);
				else if(!edge2.isSameSequenceEdge()) leftover2.connectEdgeRight(graph, edge2);
			}
			if(leftover2!=null && leftover2.getPathLength()>1) leftoverPaths.add(leftover2);
			vertexRight = vertexRightAfterRemove;
			connectEdgeRight(graph, edge);
			for(int i=loc2;i<n2;i++) {
				AssemblyEdge edge2 = path2Edges.get(i);
				if(!edge2.isSameSequenceEdge()) connectEdgeRight(graph, edge2);
			}
			this.alternativeSmallPaths.addAll(path2.alternativeSmallPaths);
			return leftoverPaths;
		}
		
		
		
		
		
		return null;
	}
	public void reverse() {
		Collections.reverse(edges);
		AssemblyVertex t = vertexLeft;
		vertexLeft = vertexRight;
		vertexRight = t;
		
	}
	public void addAlternativeSmallPath(AssemblyPath path, AssemblyEdge leftEdge, AssemblyEdge rightEdge) {
		//int idx = Math.min(locLeft.getPathPosition(), locRight.getPathPosition());
		//if(idx == locRight.getPathPosition()) path.reverse();
		alternativeSmallPaths.add(path);
	}
	public List<AssemblyPath> getAlternativeSmallPaths() {
		return alternativeSmallPaths;
	}
	public String getConsensus() {
		return consensus;
	}
	
	public void setConsensus(String consensus) {
		this.consensus = consensus;
	}
	public Map<Integer, Integer> getPathVerticesConsensusEnds() {
		return pathVerticesConsensusEnds;
	}
	public void setPathVerticesConsensusEnds(Map<Integer, Integer> pathVerticesConsensusEnds) {
		this.pathVerticesConsensusEnds = pathVerticesConsensusEnds;
	}
	public void print(PrintStream out) {
		out.print("Path: "+pathId);
		for(AssemblyEdge edge:edges) {
			out.print(edge);
		}
		System.out.println();
	}
	public void integrateAlternativeSmallPaths() {
		Map<Integer,Integer> verticesLocations = new HashMap<Integer, Integer>();
		AssemblyVertex v = vertexLeft;
		verticesLocations.put(v.getUniqueNumber(), 0);
		int i=1;
		
		for(AssemblyEdge edge:edges) {
			v = edge.getConnectingVertex(v);
			verticesLocations.put(v.getUniqueNumber(), i);
			System.out.println("Added vertex "+v+" with id: "+v.getUniqueNumber()+" pos "+i);
			i++;
		}
		int j=0;
		Set<Integer> idxsToRemove = new HashSet<Integer>();
		Map<Integer,AssemblyPath> pathsToIntegrateMap = new TreeMap<Integer, AssemblyPath>();
		System.out.println("Path start: "+vertexLeft+" length: "+edges.size()+" internal paths: "+alternativeSmallPaths.size());
		for(AssemblyPath path:alternativeSmallPaths) {
			AssemblyVertex l1 = path.getVertexLeft();
			AssemblyVertex l2 = path.getVertexRight();
			Integer p1 = verticesLocations.get(l1.getUniqueNumber());
			Integer p2 = verticesLocations.get(l2.getUniqueNumber());
			System.out.println("Next small Path borders: "+l1+" "+l2+" ids: "+l1.getUniqueNumber()+" "+l2.getUniqueNumber()+" pos: "+p1+" "+p2);
			if(p1!=null && p2!=null && Math.abs(p1-p2)==1) {
				if(p1>p2)path.reverse();
				int pos = Math.min(p1, p2);
				if(!pathsToIntegrateMap.containsKey(pos)) {
					pathsToIntegrateMap.put(pos, path);
					idxsToRemove.add(j);
				}
			}
			j++;
		}
		if(pathsToIntegrateMap.size()==0) return;
		v = vertexLeft;
		LinkedList<AssemblyEdge> newEdges = new LinkedList<AssemblyEdge>();
		List<AssemblyPath> pathsToIntegrate = new ArrayList<AssemblyPath>(pathsToIntegrateMap.values());
		i=0;
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) {
				newEdges.add(edge);
				v = edge.getConnectingVertex(v);
				continue;
			}
			
			AssemblyPath nextPath = pathsToIntegrate.get(i);
			if(v.getUniqueNumber()==nextPath.getVertexLeft().getUniqueNumber()) {
				AssemblyVertex v2 = edge.getConnectingVertex(v);
				if(v2.getUniqueNumber()==nextPath.getVertexRight().getUniqueNumber()) {
					newEdges.addAll(nextPath.edges);
				} else {
					System.err.println("Inconsistent end vertex to integrate alternative path. Expected: "+v2+" given: "+nextPath.getVertexRight());
					newEdges.add(edge);
				}
				i++;
			} else {
				newEdges.add(edge);
				
			}
			v = edge.getConnectingVertex(v);
		}
		edges = newEdges;
		List<AssemblyPath> newAlternativeSmallPaths=new ArrayList<AssemblyPath>();
		for(int k=0;k<alternativeSmallPaths.size();k++) {
			if(!idxsToRemove.contains(k)) newAlternativeSmallPaths.add(alternativeSmallPaths.get(k));
		}
		alternativeSmallPaths = newAlternativeSmallPaths;
	}
	public List<AssemblyPath> extractUnphasedPaths(AssemblyGraph graph) {
		List<AssemblyPath> answer = new ArrayList<AssemblyPath>();
		AssemblyPath unphasedPath = null;
		AssemblyVertex v = vertexLeft;
		
		for(AssemblyEdge edge:edges) {
			AssemblyVertex v2 = edge.getConnectingVertex(v);
			if(v.isInHomozygousRegion() && v2.isInHomozygousRegion()) {
				if(unphasedPath==null && edge.isSameSequenceEdge()) {
					unphasedPath = new AssemblyPath(edge);
				} else if (unphasedPath!=null && !edge.isSameSequenceEdge()) {
					unphasedPath.connectEdgeRight(graph, edge);
				}
			} else if(unphasedPath!=null) {
				answer.add(unphasedPath);
				if(getPathLength()>=5) System.out.println("Found unphased path within path starting from: "+vertexLeft+" length: "+unphasedPath.getPathLength()+" first read: "+graph.getSequence(unphasedPath.vertexLeft.getSequenceIndex()).getName()+" last: "+graph.getSequence(unphasedPath.vertexRight.getSequenceIndex()).getName());
				unphasedPath = null;
			}
			v=v2;
		}
		if(unphasedPath!=null) {
			answer.add(unphasedPath);
			if(getPathLength()>=5) System.out.println("Found unphased path within path starting from: "+vertexLeft+" length: "+unphasedPath.getPathLength()+" first read: "+graph.getSequence(unphasedPath.vertexLeft.getSequenceIndex()).getName()+" last: "+graph.getSequence(unphasedPath.vertexRight.getSequenceIndex()).getName());
		}
		return answer;
	}
}
