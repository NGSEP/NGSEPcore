package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

public class LayoutBuilderGreedy implements LayoutBuilder 
{
	private Comparator<AssemblyEdge> edgesComparator;
	@Override
	public void findPaths(AssemblyGraph graph) {
		Set<Integer> sequencesInPaths = new HashSet<>();
		List<AssemblyVertex> vertices = graph.getVertices();
		for(int i=0;i<vertices.size();i++) {
			AssemblyVertex vertex = findNextUncoveredVertex(graph, vertices,sequencesInPaths);
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
				if(nextEdgeRight==null || (nextEdgeLeft!=null && edgesComparator.compare(nextEdgeLeft, nextEdgeRight)<0)) {
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

	private AssemblyVertex findNextUncoveredVertex(AssemblyGraph graph, List<AssemblyVertex> vertices, Set<Integer> sequencesInPaths) {
		for(AssemblyVertex vertex:vertices) {
			if(!graph.isEmbedded(vertex.getSequenceIndex()) && !sequencesInPaths.contains(vertex.getSequenceIndex())) return vertex;
		}
		return null;
	}
	private AssemblyEdge findBestUncoveredEdge(AssemblyGraph graph, AssemblyVertex vertex, Set<Integer> sequencesInPaths) {
		List<AssemblyEdge> edgesVertex = new ArrayList<AssemblyEdge>(); 
		edgesVertex.addAll(graph.getEdges(vertex));
		Collections.sort(edgesVertex, edgesComparator);
		for(AssemblyEdge edge:edgesVertex) {
			AssemblyVertex connectingVertex = edge.getConnectingVertex(vertex);
			if(!sequencesInPaths.contains(connectingVertex.getSequenceIndex())&& isSafe(graph, connectingVertex,edge)) {
				return edge;
			}
		}
		return null;
	}
	private boolean isSafe(AssemblyGraph graph, AssemblyVertex vertex, AssemblyEdge bestEdge1) {
		List<AssemblyEdge> edges = graph.getEdges(vertex);
		AssemblyEdge bestEdge2 = null;
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			if(bestEdge2==null || edgesComparator.compare(edge, bestEdge2)<0) {
				bestEdge2 = edge;
			}
		}
		return bestEdge1==bestEdge2;
	}

	public void setComparator(Comparator<AssemblyEdge> comparator) {
		this.edgesComparator = comparator;
	}
}
