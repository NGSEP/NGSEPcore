package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class LayoutBuilderGreedy implements LayoutBuilder 
{
	private int minPathLength = LayoutBuilder.DEF_MIN_PATH_LENGTH;
	
	public int getMinPathLength() {
		return minPathLength;
	}
	public void setMinPathLength(int minPathLength) {
		this.minPathLength = minPathLength;
	}
	
	private Comparator<AssemblyEdge> edgesComparator;
	@Override
	public void findPaths(AssemblyGraph graph) {
		Set<Integer> sequencesInPaths = new HashSet<>();
		List<AssemblyVertex> vertices = graph.getVertices();
		for(int i=0;i<vertices.size();i++) {
			AssemblyVertex vertex = findNextUncoveredVertex(graph, vertices,sequencesInPaths);
			if(vertex==null) return;
			
			
			
			AssemblyEdge nextEdgeLeft = graph.getSameSequenceEdge (vertex);
			AssemblyEdge nextEdgeRight = nextEdgeLeft;
			AssemblyPath currentPath = new AssemblyPath(nextEdgeRight);
			sequencesInPaths.add(vertex.getSequenceIndex());
			nextEdgeLeft = findBestUncoveredEdge(graph, currentPath.getVertexLeft(), sequencesInPaths);
			nextEdgeRight = findBestUncoveredEdge(graph, currentPath.getVertexRight(), sequencesInPaths);
			while(nextEdgeLeft!=null || nextEdgeRight!=null) {
				if(nextEdgeRight==null || (nextEdgeLeft!=null && edgesComparator.compare(nextEdgeLeft, nextEdgeRight)<0)) {
					currentPath.connectEdgeLeft(graph, nextEdgeLeft);
					sequencesInPaths.add(currentPath.getVertexLeft().getSequenceIndex());
				} else if(nextEdgeRight!=null) {
					currentPath.connectEdgeRight(graph, nextEdgeRight);
					sequencesInPaths.add(currentPath.getVertexRight().getSequenceIndex());
				}
				nextEdgeLeft = findBestUncoveredEdge(graph, currentPath.getVertexLeft(), sequencesInPaths);
				nextEdgeRight = findBestUncoveredEdge(graph, currentPath.getVertexRight(), sequencesInPaths);
				//System.out.println("Vertex left "+vertexLeft.getIndex()+" vertex right: "+vertexRight.getIndex());
			}
			if(currentPath.getPathLength()>1) {
				System.err.println("Found path of length "+currentPath.getPathLength());
				graph.addPath(currentPath);
			}
		}
		
		
		
		
	}

	private AssemblyVertex findNextUncoveredVertex(AssemblyGraph graph, List<AssemblyVertex> vertices, Set<Integer> sequencesInPaths) {
		for(AssemblyVertex vertex:vertices) {
			if(!graph.isEmbedded(vertex.getSequenceIndex()) && graph.getEdges(vertex).size()>0 && !sequencesInPaths.contains(vertex.getSequenceIndex())) return vertex;
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
