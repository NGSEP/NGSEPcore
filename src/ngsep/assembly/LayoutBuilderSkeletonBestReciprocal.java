package ngsep.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import ngsep.math.Distribution;

public class LayoutBuilderSkeletonBestReciprocal implements LayoutBuilder {

	@Override
	public void findPaths(AssemblyGraph graph) {
		List<AssemblyEdge> safeReciprocalEdges = selectSafeEdges(graph);
		System.out.println("Number of safe edges: "+safeReciprocalEdges.size());
		Distribution [] edgesStats = calculateStatistics(safeReciprocalEdges);
		System.out.println("Average overlap: "+edgesStats[0].getAverage()+" average coverage shared kmers: "+edgesStats[1].getAverage());
		List<List<AssemblyEdge>> safePaths = buildPaths(graph, safeReciprocalEdges);
		System.out.println("Number of initial paths: "+safePaths.size());
		List<List<AssemblyEdge>> mergedPaths = mergePaths(graph, safePaths, edgesStats);
		for(List<AssemblyEdge> path:mergedPaths) {
			graph.addPath(path);
		}

	}

	private List<AssemblyEdge> selectSafeEdges(AssemblyGraph graph) {
		List<AssemblyEdge> allEdges = graph.getEdges();
		List<AssemblyEdge> safeEdges = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge edge:allEdges) {
			if(!edge.isSameSequenceEdge() && isRecipocalBest(graph, edge)) {
				safeEdges.add(edge);
			}
		}
		return safeEdges;
	}

	private boolean isRecipocalBest(AssemblyGraph graph, AssemblyEdge edge) {
		AssemblyVertex v1 = edge.getVertex1();
		List<AssemblyEdge> edgesV1 = graph.getEdges(v1);
		for(AssemblyEdge edgeV1:edgesV1) {
			if(edgeV1.isSameSequenceEdge() || edge==edgeV1) continue;
			//if(edgeV1.getOverlap()>=edge.getOverlap() || edgeV1.getCoverageSharedKmers()>=edge.getCoverageSharedKmers()) return false;
			if(edgeV1.getOverlap()>=edge.getOverlap()) return false;
		}
		AssemblyVertex v2 = edge.getVertex2();
		List<AssemblyEdge> edgesV2 = graph.getEdges(v2);
		for(AssemblyEdge edgeV2:edgesV2) {
			if(edgeV2.isSameSequenceEdge() || edge==edgeV2) continue;
			//if(edgeV2.getOverlap()>=edge.getOverlap() || edgeV2.getCoverageSharedKmers()>=edge.getCoverageSharedKmers()) return false;
			if(edgeV2.getOverlap()>=edge.getOverlap()) return false;
		}
		return true;
	}

	private Distribution[] calculateStatistics(List<AssemblyEdge> edges) {
		Distribution overlapDistribution = new Distribution(0, 100000, 1000);
		Distribution kmerHitCoverageDistribution = new Distribution(0, 100000, 1000);
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			overlapDistribution.processDatapoint(edge.getOverlap());
			kmerHitCoverageDistribution.processDatapoint(edge.getCoverageSharedKmers());
		}
		Distribution [] answer = {overlapDistribution, kmerHitCoverageDistribution};
		return answer;
	}

	private List<List<AssemblyEdge>> buildPaths(AssemblyGraph graph, List<AssemblyEdge> edges) {
		Map<Integer,AssemblyEdge> bestEdgesByVertex = new HashMap<Integer, AssemblyEdge>();
		Set<Integer> sequencesWithEdges = new TreeSet<Integer>(); 
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			AssemblyVertex v1 = edge.getVertex1();
			AssemblyVertex v2 = edge.getVertex2();
			bestEdgesByVertex.put(v1.getUniqueNumber(), edge);
			bestEdgesByVertex.put(v2.getUniqueNumber(), edge);
			sequencesWithEdges.add(v1.getSequenceIndex());
			sequencesWithEdges.add(v2.getSequenceIndex());
		}
		List<List<AssemblyEdge>> paths = new ArrayList<List<AssemblyEdge>>();
		Set<Integer> sequencesInPaths = new HashSet<>();
		for(int i:sequencesWithEdges) {
			if(graph.isEmbedded(i)) continue;
			if(sequencesInPaths.contains(i)) continue;
			AssemblyVertex nextVertex = graph.getVertex(i, true);
			List<AssemblyEdge> nextPath = new LinkedList<AssemblyEdge>();
			nextPath.add(graph.getSameSequenceEdge(nextVertex));
			//Expand v1
			while(nextVertex!=null) {
				AssemblyEdge nextEdge = bestEdgesByVertex.get(nextVertex.getUniqueNumber());
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
				AssemblyEdge nextEdge = bestEdgesByVertex.get(nextVertex.getUniqueNumber());
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

	private AssemblyEdge selectNextUncoveredEdge(List<AssemblyEdge> edges, Set<Integer> sequencesInPaths) {
		for(AssemblyEdge edge:edges) {
			if(!sequencesInPaths.contains(edge.getVertex1().getSequenceIndex()) && !sequencesInPaths.contains(edge.getVertex2().getSequenceIndex())) return edge;
		}
		return null;
	}

	private List<List<AssemblyEdge>> mergePaths(AssemblyGraph graph, List<List<AssemblyEdge>> paths, Distribution[] edgesStats) {
		int p = paths.size();
		AssemblyVertex [] vertices = new AssemblyVertex[2*p];
		int v=0;
		for(int i=0;i<p;i++) {
			List<AssemblyEdge> path = paths.get(i);
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
		List<AssemblyEdgePathEnd> candidateEdges = new ArrayList<AssemblyEdgePathEnd>();
		for(int i=0;i<vertices.length;i++) {
			for(int j=i+1;j<vertices.length;j++) {
				if(i%2==0 && j==i+1) continue;
				AssemblyEdge edge = graph.getEdge(vertices[i], vertices[j]);
				if(edge !=null && passFilters(edge,edgesStats)) candidateEdges.add(new AssemblyEdgePathEnd(edge, i, j));
			}
		}
		Map<Integer,AssemblyEdgePathEnd> pathEdgesByVertexId = selectEdgesToMergePaths(candidateEdges,vertices.length);
		List<AssemblyEdge> newPathsEdges = new ArrayList<AssemblyEdge>();
		for(List<AssemblyEdge> path:paths) newPathsEdges.addAll(path);
		for(AssemblyEdgePathEnd edgePathEnd:pathEdgesByVertexId.values()) newPathsEdges.add(edgePathEnd.getEdgeAssemblyGraph());
		return buildPaths(graph, newPathsEdges);
		/*List<List<AssemblyEdge>> finalPaths = new ArrayList<List<AssemblyEdge>>();
		Set<Integer> mergedPaths = new HashSet<Integer>();
		for(int i=0;i<p;i++) {
			if(mergedPaths.contains(i)) continue;
			List<AssemblyEdge> nextPath = new LinkedList<AssemblyEdge>();
			nextPath.addAll(paths.get(i));
			mergedPaths.add(i);
			//Augment left
			int v1 = 2*i;
			AssemblyEdgePathEnd nextEdge = pathEdgesByVertexId.get(v1);
			while (nextEdge!=null) {
				int v2 = nextEdge.getVertex1();
				if(v2==v1) v2 = nextEdge.getVertex2();
				if (mergedPaths.contains(v2/2)) {
					System.err.println("WARN: cycle detected merging paths");
					break;
				}
				List<AssemblyEdge> pathToMerge = paths.get(v2/2);
				if(v2%2==0) Collections.reverse(pathToMerge);
				nextPath.addAll(0, pathToMerge);
				mergedPaths.add(v2/2);
				v1=v2%2==0?v2+1:v2-1;
				nextEdge = pathEdgesByVertexId.get(v1);
			}
			//Augment right
			v1 = 2*i+1;
			nextEdge = pathEdgesByVertexId.get(v1);
			while (nextEdge!=null) {
				int v2 = nextEdge.getVertex1();
				if(v2==v1) v2 = nextEdge.getVertex2();
				if (mergedPaths.contains(v2/2)) {
					System.err.println("WARN: cycle detected merging paths");
					break;
				}
				List<AssemblyEdge> pathToMerge = paths.get(v2/2);
				if(v2%2==1) Collections.reverse(pathToMerge);
				nextPath.addAll(pathToMerge);
				mergedPaths.add(v2/2);
				v1=v2%2==0?v2+1:v2-1;
				nextEdge = pathEdgesByVertexId.get(v1);
			}
			finalPaths.add(nextPath);
		}
		return finalPaths; */
	}

	private Map<Integer,AssemblyEdgePathEnd> selectEdgesToMergePaths(List<AssemblyEdgePathEnd> candidateEdges, int length) {
		Collections.sort(candidateEdges,(e1,e2)->e2.getOverlap()-e1.getOverlap());
		int [] clusters = new int[length];
		boolean [] used = new boolean[length];
		Arrays.fill(used, false);
		for(int i=0;i<length;i++) {
			clusters[i] = i/2;
		}
		Map<Integer,AssemblyEdgePathEnd> answer = new HashMap<Integer, AssemblyEdgePathEnd>();
		for(AssemblyEdgePathEnd nextEdge:candidateEdges) {
			int v1 = nextEdge.getVertex1();
			if(used[v1]) continue;
			int v2 = nextEdge.getVertex2();
			if(used[v2]) continue;
			int c1 = clusters[v1];
			int c2 = clusters[v2];
			if(c1!=c2) {
				answer.put(v1, nextEdge);
				answer.put(v2, nextEdge);
				used[v1] =used[v2] = true;
				for(int i=0;i<length;i++) {
					if(clusters[i]==c2) clusters[i] = c1;
				}
			}
		}
		return answer;
	}

	private boolean passFilters(AssemblyEdge edge, Distribution[] edgesStats) {
		// TODO Auto-generated method stub
		return true;
	}

}
class AssemblyEdgePathEnd  {

	private AssemblyEdge edgeAssemblyGraph;	
	private int vertex1 = 0;
	private int vertex2 = 0;
	public AssemblyEdgePathEnd(AssemblyEdge edge, int v1, int v2 ) {
		edgeAssemblyGraph = edge;
		//Numbers encoding path index and start / end status
		vertex1 = v1;
		vertex2 = v2;
	}
	public AssemblyEdge getEdgeAssemblyGraph() {
		return edgeAssemblyGraph;
	}
	public int getVertex1() {
		return vertex1;
	}
	public int getVertex2() {
		return vertex2;
	}
	public int getOverlap() {
		return edgeAssemblyGraph.getOverlap();
	}
	public int getCoverageSharedKmers() {
		return edgeAssemblyGraph.getCoverageSharedKmers();
	}
	
	

}
