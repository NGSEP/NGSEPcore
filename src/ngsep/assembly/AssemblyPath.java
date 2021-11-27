package ngsep.assembly;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class AssemblyPath {
	private LinkedList<AssemblyEdge> edges = new LinkedList<AssemblyEdge>();
	private AssemblyVertex vertexLeft;
	private AssemblyVertex vertexRight;
	private int pathId;
	private String sequenceName;
	private List<AssemblyPath> alternativeSmallPaths=new ArrayList<AssemblyPath>();
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
		System.err.println("Trying to merge paths with ends: "+lastVertex1+" "+firstVertex2+" Merging direct edge: "+c1);
		if(c1!=null && (n1==1 || c1.getCost()<2*lastCost) && c1.getCost()<2*nextCost && c1.getOverlap()>0.5*lastOverlap && (n2==1 || c1.getOverlap()>0.5*nextOverlap)) {
			System.err.println("Merging paths with ends: "+lastVertex1+" "+firstVertex2+" Merging direct edge: "+c1);
			connectEdgeRight(graph, c1);
			for(AssemblyEdge edge:path2Edges) {
				if(!edge.isSameSequenceEdge()) connectEdgeRight(graph, edge);
			}
			return true;
		}
		if(secondLastVertex1 == null || secondVertex2==null) return false; 
		AssemblyEdge c2 = graph.getEdge(lastVertex1, secondVertex2);
		AssemblyEdge c3 = graph.getEdge(secondLastVertex1, firstVertex2);
		System.err.println("Edges to merge paths: "+c2+" "+c3+" overlaps: "+lastOverlap+" "+nextOverlap+" costs: "+lastCost+" "+nextCost);
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
			return true;
		}
		return false;
		
	}
	public void reverse() {
		Collections.reverse(edges);
		AssemblyVertex t = vertexLeft;
		vertexLeft = vertexRight;
		vertexRight = t;
		
	}
	public void addAlternativeSmallPath(AssemblyPath path) {
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
		out.println("Path: "+pathId);
		for(AssemblyEdge edge:edges) {
			out.println(edge);
		}
	}
	
	
}
