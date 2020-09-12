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

import JSci.maths.statistics.NormalDistribution;
import ngsep.math.Distribution;
import ngsep.math.PhredScoreHelper;

public class LayoutBuilderKruskalPath implements LayoutBuilder {

	@Override
	public void findPaths(AssemblyGraph graph) {
		Distribution initialDegreesDist = new Distribution(0, 10000, 1);
		List<AssemblyVertex> vertices = graph.getVertices();
		List<AssemblyEdge> allEdges = graph.getEdges();
		for(AssemblyVertex v:vertices) {
			if(!graph.isEmbedded(v.getSequenceIndex())) initialDegreesDist.processDatapoint(v.getDegreeUnfilteredGraph());
		}
		System.out.println("Degree average: "+initialDegreesDist.getAverage()+" variance "+initialDegreesDist.getVariance());
		NormalDistribution distDegrees = new NormalDistribution(initialDegreesDist.getAverage(), initialDegreesDist.getVariance());
		
		Set<Integer> repetitiveVertices = new HashSet<Integer>();
		for(AssemblyVertex v:vertices) {
			if(isRepetivive(v, distDegrees)) repetitiveVertices.add(v.getUniqueNumber());
		}
		System.out.println("Total vertices for layout: "+initialDegreesDist.getCount()+" Number of repetitive vertices: "+repetitiveVertices.size());
		List<AssemblyEdge> pathEdges = selectSafeEdges(graph, repetitiveVertices);
		System.out.println("Number of safe edges: "+pathEdges.size());
		
		List<LinkedList<AssemblyEdge>> safePaths = buildPaths(graph, pathEdges);
		//Map<Integer,Integer> vertexPathIds = calculateVertexPathIdsMap(safePaths);
		System.out.println("Number of initial paths: "+safePaths.size());
		Distribution [] edgesStats = calculateStatistics(pathEdges, allEdges, repetitiveVertices);
		System.out.println("Average overlap TP: "+edgesStats[0].getAverage()+" SD: "+Math.sqrt(edgesStats[0].getVariance())+ " Total: "+edgesStats[0].getCount());
		System.out.println("Average coverage shared kmers TP: "+edgesStats[1].getAverage()+" SD: "+Math.sqrt(edgesStats[1].getVariance())+ " Total: "+edgesStats[1].getCount());
		System.out.println("Average weighted coverage shared kmers TP: "+edgesStats[2].getAverage()+" SD: "+Math.sqrt(edgesStats[2].getVariance())+ " Total: "+edgesStats[2].getCount());
		System.out.println("Average weighted coverage proportion TP: "+edgesStats[3].getAverage()+" SD: "+Math.sqrt(edgesStats[3].getVariance())+ " Total: "+edgesStats[3].getCount());
		addConnectingEdges(graph, safePaths, pathEdges, edgesStats);
		List<LinkedList<AssemblyEdge>> paths = buildPaths(graph, pathEdges);
		for(LinkedList<AssemblyEdge> path:paths) {
			if(path.size()<=5) continue;
			graph.addPath(path);
		}

	}
	private boolean isRepetivive(AssemblyVertex vertex, NormalDistribution distDegrees) {
		int degree = vertex.getDegreeUnfilteredGraph();
		double pValue = distDegrees.cumulative(degree);
		boolean repetitive = pValue>0.999;
		//if(repetitive) System.out.println("Repetitive vertex: "+vertex+" initial degree: "+vertex.getDegreeUnfilteredGraph()+" average: "+distDegrees.getMean()+" sd "+Math.sqrt(distDegrees.getVariance()));
		return repetitive;
	}
	
	private List<AssemblyEdge> selectSafeEdges(AssemblyGraph graph, Set<Integer> repetitiveVertices ) {
		List<AssemblyEdge> allEdges = graph.getEdges();
		List<AssemblyEdge> safeEdges = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge edge:allEdges) {
			if(edge.isSameSequenceEdge()) continue;
			boolean r1 = repetitiveVertices.contains(edge.getVertex1().getUniqueNumber());
			boolean r2 = repetitiveVertices.contains(edge.getVertex2().getUniqueNumber()); 
			int d1 = graph.getEdges(edge.getVertex1()).size();
			int d2 = graph.getEdges(edge.getVertex2()).size();
			if(logEdge(edge)) System.out.println("Select safe edges. repetitive: "+r1+" "+r2+" initial degrees "+edge.getVertex1().getDegreeUnfilteredGraph()+" "+edge.getVertex2().getDegreeUnfilteredGraph()+" current degrees "+d1+" "+d2+" reciprocal best: "+isRecipocalBest(graph, edge)+" edge: "+edge);
			if((d1==2 && d2==2) || (!r1 && !r2 && isRecipocalBest(graph, edge))) {
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
			if(edgeV1.getOverlap()>=edge.getOverlap() || edgeV1.getWeightedCoverageSharedKmers()>=edge.getWeightedCoverageSharedKmers()) return false;
			//if(edgeV1.getOverlap()>=edge.getOverlap()) return false;
		}
		AssemblyVertex v2 = edge.getVertex2();
		List<AssemblyEdge> edgesV2 = graph.getEdges(v2);
		for(AssemblyEdge edgeV2:edgesV2) {
			if(edgeV2.isSameSequenceEdge() || edge==edgeV2) continue;
			if(edgeV2.getOverlap()>=edge.getOverlap() || edgeV2.getWeightedCoverageSharedKmers()>=edge.getWeightedCoverageSharedKmers()) return false;
			//if(edgeV2.getOverlap()>=edge.getOverlap()) return false;
		}
		return true;
	}
	private Map<Integer, Integer> calculateVertexPathIdsMap(List<List<AssemblyEdge>> paths) {
		Map<Integer, Integer> answer = new HashMap<Integer, Integer>();
		for(int i=0;i<paths.size();i++) {
			List<AssemblyEdge> path = paths.get(i);
			for(AssemblyEdge edge:path) {
				answer.put(edge.getVertex1().getUniqueNumber(), i);
				answer.put(edge.getVertex2().getUniqueNumber(), i);
			}
		}
		return answer;
	}
	private Distribution[] calculateStatistics(List<AssemblyEdge> safeEdges, List<AssemblyEdge> allEdges, Set<Integer> repetitiveVertices) {
		Distribution overlapDistributionTP = new Distribution(0, 100000, 1);
		Distribution kmerHitCoverageDistributionTP = new Distribution(0, 100000, 1);
		Distribution kmerHitWCovDistributionTP = new Distribution(0, 100000, 1);
		Distribution coverageProportionDistributionTP = new Distribution(0, 1.5, 0.01);
		for(AssemblyEdge edge:safeEdges) {
			overlapDistributionTP.processDatapoint(edge.getOverlap());
			kmerHitCoverageDistributionTP.processDatapoint(edge.getCoverageSharedKmers());
			kmerHitWCovDistributionTP.processDatapoint(edge.getWeightedCoverageSharedKmers());
			coverageProportionDistributionTP.processDatapoint((double)edge.getWeightedCoverageSharedKmers()/edge.getOverlap());
		}
		Distribution [] answer = {overlapDistributionTP, kmerHitCoverageDistributionTP, kmerHitWCovDistributionTP, coverageProportionDistributionTP};
		return answer;
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
			if(graph.isEmbedded(i)) continue;
			if(graph.getVertex(i, true)==null || graph.getVertex(i, false)==null) continue;
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

	private void addConnectingEdges(AssemblyGraph graph, List<LinkedList<AssemblyEdge>> paths, List<AssemblyEdge> pathEdges, Distribution[] edgesStats) {
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
		List<AssemblyEdgePathEnd> candidateEdges = new ArrayList<AssemblyEdgePathEnd>();
		for(int i=0;i<vertices.length;i++) {
			for(int j=i+1;j<vertices.length;j++) {
				if(i%2==0 && j==i+1) continue;
				AssemblyEdge edge = graph.getEdge(vertices[i], vertices[j]);
				if(edge !=null && passFilters(edge,edgesStats)) candidateEdges.add(new AssemblyEdgePathEnd(edge, i, j));
			}
		}
		Collections.sort(candidateEdges,(e1,e2)->calculateCost(e1,edgesStats)-calculateCost(e2,edgesStats));
		List<AssemblyEdgePathEnd> pathEdgesByVertexId = selectEdgesToMergePaths(candidateEdges,vertices.length, edgesStats);
		for(AssemblyEdgePathEnd edgeP:pathEdgesByVertexId) pathEdges.add(edgeP.getEdgeAssemblyGraph());
	}

	private int calculateCost(AssemblyEdgePathEnd edge, Distribution[] edgesStats) {
		Distribution overlapTP = edgesStats[0];
		Distribution covTP = edgesStats[1];
		Distribution wCovTP = edgesStats[2];
		Distribution wCovPropTP = edgesStats[3];
		//Distribution covPropTP = edgesStats[2];
		//double prop = (double)edge.getCoverageSharedKmers()/edge.getOverlap();
		//int cost = (int)Math.round(1000*(1.5-prop));
		//return cost;
		NormalDistribution noTP = new NormalDistribution(overlapTP.getAverage(),overlapTP.getVariance());
		NormalDistribution ncTP = new NormalDistribution(covTP.getAverage(),covTP.getVariance());
		NormalDistribution nwcTP = new NormalDistribution(wCovTP.getAverage(),wCovTP.getVariance());
		NormalDistribution nwcpTP = new NormalDistribution(wCovPropTP.getAverage(),wCovPropTP.getVariance());
		double pValueOTP = noTP.cumulative(edge.getOverlap());
		int cost1 = PhredScoreHelper.calculatePhredScore(pValueOTP);
		double pValueCTP = ncTP.cumulative(edge.getCoverageSharedKmers());
		int cost2 = PhredScoreHelper.calculatePhredScore(pValueCTP);
		double pValueWCTP = nwcTP.cumulative(edge.getEdgeAssemblyGraph().getWeightedCoverageSharedKmers());
		int cost3 = PhredScoreHelper.calculatePhredScore(pValueWCTP);
		double pValueWCPTP = nwcpTP.cumulative((double)edge.getEdgeAssemblyGraph().getWeightedCoverageSharedKmers()/edge.getOverlap());
		int cost4 = PhredScoreHelper.calculatePhredScore(pValueWCPTP);
		if( logEdge(edge.getEdgeAssemblyGraph())) System.out.println("CalculateCost. Pvalues "+pValueOTP+" "+pValueCTP+" "+pValueWCTP+" "+pValueWCPTP+" costs: "+cost1+" "+cost2+" "+cost3+" "+cost4+" sum: " +(cost1+cost2+cost3)+ " Edge: "+edge.getEdgeAssemblyGraph());
		//return cost1+cost3+cost4;
		return cost1+cost3;
		//return cost3;
		//return cost1+cost2;
	}
	private boolean logEdge(AssemblyEdge edge) {
		int n = -1;
		return edge.getVertex1().getSequenceIndex()==n || edge.getVertex2().getSequenceIndex()==n;
		//return false;
	}
	private List<AssemblyEdgePathEnd> selectEdgesToMergePaths(List<AssemblyEdgePathEnd> candidateEdges, int length, Distribution[] edgesStats) {
		
		int [] clusters = new int[length];
		boolean [] used = new boolean[length];
		Arrays.fill(used, false);
		for(int i=0;i<length;i++) {
			clusters[i] = i/2;
		}
		List<AssemblyEdgePathEnd> answer = new ArrayList<AssemblyEdgePathEnd>();
		for(AssemblyEdgePathEnd nextEdge:candidateEdges) {
			//if(nextEdge.getEdgeAssemblyGraph().getVertex1().getUniqueNumber()==2854) System.out.println("score: "+calculateCost(nextEdge, edgesStats)+" edge: "+nextEdge.getEdgeAssemblyGraph());
			//if(answer.size()<20) System.out.println("score: "+calculateCost(nextEdge, edgesStats)+ " edge: "+nextEdge.getEdgeAssemblyGraph());
			int v1 = nextEdge.getVertex1();
			if(used[v1]) continue;
			int v2 = nextEdge.getVertex2();
			if(used[v2]) continue;
			int c1 = clusters[v1];
			int c2 = clusters[v2];
			if(c1!=c2) {
				answer.add(nextEdge);
				used[v1] =used[v2] = true;
				for(int i=0;i<length;i++) {
					if(clusters[i]==c2) clusters[i] = c1;
				}
			}
		}
		return answer;
	}

	private boolean passFilters(AssemblyEdge edge, Distribution[] edgesStats) {
		/*Distribution overlapTP = edgesStats[0];
		Distribution covTP = edgesStats[1];
		Distribution overlapFP = edgesStats[2];
		Distribution covFP = edgesStats[3];
		if(overlapTP.getCount()==0 || overlapFP.getCount()==0) return true;
		NormalDistribution noTP = new NormalDistribution(overlapTP.getAverage(),overlapTP.getVariance());
		NormalDistribution ncTP = new NormalDistribution(covTP.getAverage(),covTP.getVariance());
		NormalDistribution noFP = new NormalDistribution(overlapFP.getAverage(),overlapFP.getVariance());
		NormalDistribution ncFP = new NormalDistribution(covFP.getAverage()-covTP.getAverage()/2,covFP.getVariance());
		double pValueOTP = noTP.cumulative(edge.getOverlap());
		if(pValueOTP>0.5) pValueOTP=1-pValueOTP;
		double pValueCTP = ncTP.cumulative(edge.getCoverageSharedKmers());
		if(pValueCTP>0.5) pValueCTP=1-pValueCTP;
		double pValueOFP = noFP.cumulative(edge.getOverlap());
		if(pValueOFP>0.5) pValueOFP=1-pValueOFP;
		double pValueCFP = ncFP.cumulative(edge.getCoverageSharedKmers());
		if(pValueCFP>0.5) pValueCFP=1-pValueCFP;
		boolean passFilter = pValueOTP*pValueCTP/(pValueOFP*pValueCFP) > 2;
		if(!passFilter) System.out.println("False positive edge "+edge+" p values: "+pValueOTP+" "+pValueCTP+" "+pValueOFP+" "+pValueCFP);
		//return pValueTP>=0.01;
		//double pValueFP = overlapFP.getEmpiricalPvalue(edge.getOverlap())*covFP.getEmpiricalPvalue(edge.getCoverageSharedKmers());
		//return pValueTP/pValueFP>=2;
		return passFilter;*/
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
