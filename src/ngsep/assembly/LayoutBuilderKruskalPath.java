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

	private int minPathLength = 6;
	private boolean useIndels = false;
	
	public int getMinPathLength() {
		return minPathLength;
	}
	public void setMinPathLength(int minPathLength) {
		this.minPathLength = minPathLength;
	}
	
	
	public boolean isUseIndels() {
		return useIndels;
	}
	public void setUseIndels(boolean useIndels) {
		this.useIndels = useIndels;
	}
	@Override
	public void findPaths(AssemblyGraph graph) {
		filterEdgesAndEmbedded(graph, 0.3);
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
		System.out.println("Average Evidence proportion TP: "+edgesStats[4].getAverage()+" SD: "+Math.sqrt(edgesStats[4].getVariance())+ " Total: "+edgesStats[4].getCount());
		System.out.println("Average indels kbp TP: "+edgesStats[5].getAverage()+" SD: "+Math.sqrt(edgesStats[5].getVariance())+ " Total: "+edgesStats[5].getCount());
		System.out.println("Distribution indels kbp");
		edgesStats[5].printDistribution(System.out);
		//Debug
		for(AssemblyEdge edge: graph.getEdges()) calculateCost(edge, edgesStats);
		//Algorithms to resolve conflicts between almost safe close edges
		addEdges2(graph, safePaths, pathEdges);
		System.out.println("Number of edges after second algorithm: "+pathEdges.size());
		safePaths = buildPaths(graph, pathEdges);
		System.out.println("Updated number of paths: "+safePaths.size());
		
		addEdges3(graph, safePaths, pathEdges);
		System.out.println("Number of edges after third algorithm: "+pathEdges.size());
		safePaths = buildPaths(graph, pathEdges);
		System.out.println("Updated number of paths: "+safePaths.size());
		
		addConnectingEdges(graph, safePaths, pathEdges, edgesStats);
		List<LinkedList<AssemblyEdge>> paths = buildPaths(graph, pathEdges);
		for(LinkedList<AssemblyEdge> path:paths) {
			if(path.size()<minPathLength) continue;
			graph.addPath(path);
		}
		System.out.println("Final number of paths: "+graph.getPaths().size());
		System.out.println("Estimated N statistics");
		long [] stats = graph.estimateNStatisticsFromPaths();
		if(stats!=null) NStatisticsCalculator.printNStatistics(stats, System.out);

	}
	
	//Initial edge filtering
	public void filterEdgesAndEmbedded(AssemblyGraph graph, double minScoreProportionEdges) {
		Distribution lengthsDistribution = new Distribution(0, graph.getSequenceLength(0), 1);
		int n = graph.getNumSequences();
		for(int i=0;i<n;i++) lengthsDistribution.processDatapoint(graph.getSequenceLength(i));
		/*Distribution evidenceProportionEmbedded = new Distribution(0, 1, 0.01);
		Distribution cskProportionSelfEmbedded = new Distribution(0, 1, 0.01);
		Distribution wcskProportionSelfEmbedded = new Distribution(0, 1, 0.01);
		for(int seqId:embeddedMapBySequence.keySet()) {
			int length = getSequenceLength(seqId);
			AssemblyEdge edge = getSameSequenceEdge(seqId);
			if(edge == null) continue;
			int selfCSK = edge.getCoverageSharedKmers();
			int selfWCSK = edge.getWeightedCoverageSharedKmers();
			List<AssemblyEmbedded> relations = embeddedMapBySequence.get(seqId);
			for(AssemblyEmbedded embedded:relations) {
				double evidenceLength = embedded.getHostEvidenceEnd()-embedded.getHostEvidenceStart();
				evidenceProportionEmbedded.processDatapoint(evidenceLength/length);
				cskProportionSelfEmbedded.processDatapoint((double)embedded.getCoverageSharedKmers()/selfCSK);
				wcskProportionSelfEmbedded.processDatapoint((double)embedded.getCoverageSharedKmers()/selfWCSK);
			}
		}
		System.out.println("Proportion of evidence vs read length for embedded relationships");
		evidenceProportionEmbedded.printDistribution(System.out);
		System.out.println("Proportion of CSK vs self CSK for embedded relationships");
		cskProportionSelfEmbedded.printDistribution(System.out);
		System.out.println("Proportion of WCSK vs self WCSK for embedded relationships");
		wcskProportionSelfEmbedded.printDistribution(System.out);*/
		int medianLength = graph.getMedianLength();
		System.out.println("Median read length: "+medianLength);
		int numEmbedded = 0;
		for (int seqId = n-1; seqId >=0; seqId--) {
			if(filterEdgesAndEmbedded(graph, seqId,medianLength, minScoreProportionEdges)) numEmbedded++;
		}
		System.out.println("Filtered edges and embedded. Final number of embedded sequences: "+numEmbedded);
		graph.pruneEmbeddedSequences();
		System.out.println("Prunned embedded sequences");
		//filterEdgesCloseRelationships();
	}

	private boolean filterEdgesAndEmbedded(AssemblyGraph graph, int sequenceId,int medianLength, double minScoreProportionEdges) {
		int debugIdx = -1;
		int sequenceLength = graph.getSequenceLength(sequenceId);
		AssemblyVertex vS = graph.getVertex(sequenceId, true);
		AssemblyVertex vE = graph.getVertex(sequenceId, false);
		if(vS==null || vE==null) return false;
		if(sequenceId == debugIdx) System.out.println("Filtered edges with abnormal features");
		List<AssemblyEdge> edgesS = new ArrayList<AssemblyEdge>();
		if(vS!=null) edgesS.addAll(graph.getEdges(vS));
		
		double maxScoreSE = 0;
		double maxScoreSF = 0;
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			maxScoreSF = Math.max(maxScoreSF, calculateScoreForEdgeFiltering(edge));
			//int connectingLength = getSequenceLength(edge.getConnectingVertex(vS).getSequenceIndex());
			/*if(connectingLength<1.2*sequenceLength && edge.getOverlap() >0.8*sequenceLength)*/ maxScoreSE = Math.max(maxScoreSE, calculateScoreForEmbedded(edge));
		}
		List<AssemblyEdge> edgesE = new ArrayList<AssemblyEdge>();
		if(vE!=null) edgesE.addAll(graph.getEdges(vE));
		double maxScoreEE = 0;
		double maxScoreEF = 0;			
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			maxScoreEF = Math.max(maxScoreEF, calculateScoreForEdgeFiltering(edge));
			//int connectingLength = getSequenceLength(edge.getConnectingVertex(vE).getSequenceIndex());
			/*if(connectingLength<1.5*sequenceLength && edge.getOverlap() >0.8*sequenceLength)*/ maxScoreEE = Math.max(maxScoreEE, calculateScoreForEmbedded(edge));
		}
		double minScoreEdges = minScoreProportionEdges*Math.max(maxScoreSF, maxScoreEF);
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScoreForEdgeFiltering(edge);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next edge start "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber()+" overlap: "+edge.getOverlap()+" score: "+score+" Max score start: "+maxScoreSF+" limit: "+minScoreEdges);
			if(score <500 || (score < maxScoreSF && score < minScoreEdges)) {
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				graph.removeEdge(edge);
			}
		}
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Initial edges start "+edgesS.size()+" Max scores start: "+maxScoreSF+" "+maxScoreSE);
		
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScoreForEdgeFiltering(edge);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next edge end "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber()+" overlap: "+edge.getOverlap()+" score: "+score+" Max score end: "+maxScoreEF+" limit: "+minScoreEdges);
			if(score < 500 || (score < maxScoreEF && score < minScoreEdges)) {
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				graph.removeEdge(edge);
			}
		}
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Initial edges end "+edgesE.size()+" Max scores end: "+maxScoreEF+" "+maxScoreEE);
		
		double medianRelationship = 1.0*sequenceLength/(double)medianLength;
		
		//double minScoreProportionEmbedded = 0.8;
		//double cumulative = lengthsDistribution.getCumulativeCount(sequenceLength)/lengthsDistribution.getCount();
		double minScoreProportionEmbedded = Math.min(0.8, 0.5*medianRelationship);
		//double minScoreProportionEmbedded = 0.8*cumulative;
		if(minScoreProportionEmbedded<0.5) minScoreProportionEmbedded = 0.5;
		//if(medianRelationship>1 && minScoreProportionEmbedded<0.7) minScoreProportionEmbedded = 0.7;
		
		double maxScoreFilterEmbedded = minScoreProportionEmbedded*Math.max(maxScoreSE, maxScoreEE);
	
		AssemblyEdge sameSequenceEdge = graph.getSameSequenceEdge(sequenceId);
		double sameSeqCSK = sameSequenceEdge.getCoverageSharedKmers();
		double sameSeqWCSK = sameSequenceEdge.getCoverageSharedKmers();
		List<AssemblyEmbedded> embeddedList= new ArrayList<AssemblyEmbedded>();
		embeddedList.addAll(graph.getEmbeddedBySequenceId(sequenceId));
		if(embeddedList.size()==0) return false;
		double maxEvidencePropEmbedded = 0;
		AssemblyEmbedded embeddedMax = null;
		double maxScoreEmbedded = -1;
		int countPass = 0;
		for(AssemblyEmbedded embedded:embeddedList) {
			maxEvidencePropEmbedded = Math.max(maxEvidencePropEmbedded, embedded.calculateEvidenceProportion());
			double CSKprop = (double)embedded.getCoverageSharedKmers()/sameSeqCSK;
			double WCSKprop = (double)embedded.getWeightedCoverageSharedKmers()/sameSeqWCSK;
			double score = calculateScore(embedded);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next embedded "+embedded+" score: "+score+" evProp: "+embedded.calculateEvidenceProportion()+" CSK prop "+CSKprop+" WCSK prop: "+WCSKprop+" Indels: "+embedded.getNumIndels()+" IKBP: "+embedded.getIndelsPerKbp());
			
			if(embeddedMax==null || maxScoreEmbedded<score) {
				maxScoreEmbedded = score;
				embeddedMax = embedded;
			}
			
			//if(evidenceProp*CSKprop >=0.25) countPass++;
			if(embedded.calculateEvidenceProportion() >=0.95) countPass++;
		}
			
		//Score proportion filter calculation
		double maxScorePropEmbedded = maxScoreEmbedded/calculateScoreForEmbedded(sameSequenceEdge);
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Median relationship: "+medianRelationship+" evidence proportion "+ maxEvidencePropEmbedded +" max score embedded "+maxScoreEmbedded+" same seq score: "+calculateScoreForEmbedded(sameSequenceEdge)+ " max score prop self: "+maxScorePropEmbedded+" count pass: "+countPass);
		
		//if(countPass==0) {
		if(maxScoreEmbedded<maxScoreFilterEmbedded) {
		//if(maxEvidencePropEmbedded<evidenceProportionThreshold) {
		//if(maxEvidencePropEmbedded<minProportionEmbedded || maxScorePropEmbedded<0.5*minProportionEmbedded) {
		//if(maxScoreEmbedded<maxScoreFilterEmbedded || maxScoreEmbedded < 0.2*calculateScoreForEmbedded(sameSequenceEdge)) {
		//if(maxScoreEmbedded < 0.2*calculateScoreForEmbedded(sameSequenceEdge)) {
			//Replace embedded relationships with edges to make the sequence not embedded
			for(AssemblyEmbedded embedded:embeddedList) {
				graph.removeEmbedded(embedded);
				if(sequenceId == debugIdx) System.out.println("Adding edge replacing embedded "+embedded.getHostId()+" limits: "+embedded.getHostStart()+" "+embedded.getHostEnd()+" host length: "+graph.getSequenceLength(embedded.getHostId())+"score: "+calculateScore(embedded));
				addEdgeFromEmbedded(graph, embedded);
			}
			return false;
		} else {
			if(sequenceId == debugIdx) System.out.println("Sequence is embedded ");
			for(AssemblyEmbedded embedded:embeddedList) {
				if(embedded!=embeddedMax) {
					graph.removeEmbedded(embedded);
					//if(sequenceId == debugIdx) System.out.println("Assembly graph. Removed embedded host: "+embedded.getHostId()+" Embedded relations: "+embeddedMapBySequence.get(sequenceId)+" is embedded: "+isEmbedded(sequenceId));
				}
			}
			return true;
		}
	}

	private void addEdgeFromEmbedded(AssemblyGraph graph, AssemblyEmbedded embedded) {
		int distanceStart = embedded.getHostStart();
		int distanceEnd = graph.getSequenceLength(embedded.getHostId())-embedded.getHostEnd();
		AssemblyVertex vertexHost=null;
		AssemblyVertex vertexEmbedded=null;
		if(distanceStart<0.5*distanceEnd) {
			vertexHost = graph.getVertex(embedded.getHostId(), true);
			vertexEmbedded = graph.getVertex(embedded.getSequenceId(), embedded.isReverse());
		} else if (distanceEnd<0.5*distanceStart) {
			vertexHost = graph.getVertex(embedded.getHostId(), false);
			vertexEmbedded = graph.getVertex(embedded.getSequenceId(), !embedded.isReverse());
		}
		if(vertexHost==null || vertexEmbedded==null) return;
		int overlap = Math.max(embedded.getCoverageSharedKmers(), embedded.getSequenceEvidenceEnd()-embedded.getSequenceEvidenceStart());
		AssemblyEdge edge = new AssemblyEdge(vertexHost, vertexEmbedded, overlap);
		edge.setOverlapStandardDeviation(embedded.getHostStartStandardDeviation());
		edge.setWeightedCoverageSharedKmers(embedded.getWeightedCoverageSharedKmers());
		edge.setCoverageSharedKmers(embedded.getCoverageSharedKmers());
		edge.setNumSharedKmers(embedded.getNumSharedKmers());
		edge.setRawKmerHits(embedded.getRawKmerHits());
		edge.setRawKmerHitsSubjectStartSD(embedded.getRawKmerHitsSubjectStartSD());
		edge.setVertex1EvidenceStart(embedded.getHostEvidenceStart());
		edge.setVertex1EvidenceEnd(embedded.getHostEvidenceEnd());
		edge.setVertex2EvidenceStart(embedded.getSequenceEvidenceStart());
		edge.setVertex2EvidenceEnd(embedded.getSequenceEvidenceEnd());
		edge.setNumIndels(embedded.getNumIndels());
		edge.setNumMismatches(embedded.getNumMismatches());
		graph.addEdge(edge);
	}

	private double calculateScore(AssemblyEmbedded embedded) {
		double evProp = embedded.calculateEvidenceProportion();
		double indelsScore = Math.max(1,embedded.getIndelsPerKbp()-15);
		//return embedded.getCoverageSharedKmers();
		//return embedded.getWeightedCoverageSharedKmers();
		//return embedded.getRawKmerHits();
		return embedded.getWeightedCoverageSharedKmers()*evProp*evProp;
		//return embedded.getWeightedCoverageSharedKmers()*evProp*evProp/indelsScore;
	}

	private double calculateScoreForEmbedded(AssemblyEdge edge) {
		double evProp = edge.calculateEvidenceProportion();
		double indelsScore = Math.max(1,edge.getIndelsPerKbp()-15);
		//return edge.getCoverageSharedKmers();
		//return edge.getWeightedCoverageSharedKmers();
		//return edge.getRawKmerHits();
		return edge.getWeightedCoverageSharedKmers()*evProp*evProp;
		//return edge.getWeightedCoverageSharedKmers()*evProp*evProp/indelsScore;
	}
	
	private double calculateScoreForEdgeFiltering(AssemblyEdge edge) {
		double evProp = edge.calculateEvidenceProportion();
		//return edge.getCoverageSharedKmers();
		//return edge.getRawKmerHits();
		return edge.getWeightedCoverageSharedKmers()*evProp*evProp;
	}
	
	private boolean isRepetivive(AssemblyVertex vertex, NormalDistribution distDegrees) {
		int degree = vertex.getDegreeUnfilteredGraph();
		double pValue = distDegrees.cumulative(degree);
		boolean repetitive = pValue>0.999;
		//if(repetitive) System.out.println("Repetitive vertex: "+vertex+" initial degree: "+vertex.getDegreeUnfilteredGraph()+" average: "+distDegrees.getMean()+" sd "+Math.sqrt(distDegrees.getVariance()));
		return repetitive;
	}
	private int calculateScoreForEdgeLabeling(AssemblyEdge edge) {
		double evProp = edge.calculateEvidenceProportion();
		double indelsScore = Math.max(1,edge.getIndelsPerKbp()-5);
		double score = edge.getWeightedCoverageSharedKmers()*evProp*evProp;
		//double score = edge.getWeightedCoverageSharedKmers()*evProp*evProp/indelsScore;
		return (int)Math.round(score);
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
		if(!isBestEdge(graph, v1, edge)) return false;
		AssemblyVertex v2 = edge.getVertex2();
		if(!isBestEdge(graph, v2, edge)) return false;
		return true;
	}
	private boolean isBestEdge (AssemblyGraph graph, AssemblyVertex v, AssemblyEdge edgeT) {
		List<AssemblyEdge> edges = graph.getEdges(v);
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge() || edge==edgeT) continue;
			if(edge.getOverlap()>=edgeT.getOverlap()) return false;
			if(calculateScoreForEdgeLabeling(edge)>=calculateScoreForEdgeLabeling(edgeT)) return false;
		}
		return true;
	}
	/*private Map<Integer, Integer> calculateVertexPathIdsMap(List<List<AssemblyEdge>> paths) {
		Map<Integer, Integer> answer = new HashMap<Integer, Integer>();
		for(int i=0;i<paths.size();i++) {
			List<AssemblyEdge> path = paths.get(i);
			for(AssemblyEdge edge:path) {
				answer.put(edge.getVertex1().getUniqueNumber(), i);
				answer.put(edge.getVertex2().getUniqueNumber(), i);
			}
		}
		return answer;
	}*/
	private Distribution[] calculateStatistics(List<AssemblyEdge> safeEdges, List<AssemblyEdge> allEdges, Set<Integer> repetitiveVertices) {
		Distribution overlapDistributionTP = new Distribution(0, 100000, 1);
		Distribution kmerHitCoverageDistributionTP = new Distribution(0, 100000, 1);
		Distribution kmerHitWCovDistributionTP = new Distribution(0, 100000, 1);
		Distribution coverageProportionDistributionTP = new Distribution(0, 1.5, 0.01);
		Distribution evidenceProportionDistributionTP = new Distribution(0, 1.1, 0.01);
		Distribution indelsKbpDistributionTP = new Distribution(0, 300, 1);
		
		for(AssemblyEdge edge:safeEdges) {
			if (edge.isSameSequenceEdge()) continue;
			double overlap = edge.getOverlap();
			overlapDistributionTP.processDatapoint(overlap);
			kmerHitCoverageDistributionTP.processDatapoint(edge.getCoverageSharedKmers());
			kmerHitWCovDistributionTP.processDatapoint(edge.getWeightedCoverageSharedKmers());
			coverageProportionDistributionTP.processDatapoint((double)edge.getWeightedCoverageSharedKmers()/overlap);
			evidenceProportionDistributionTP.processDatapoint(edge.calculateEvidenceProportion());
			indelsKbpDistributionTP.processDatapoint(edge.getIndelsPerKbp());
		}
		Distribution [] answer = {overlapDistributionTP, kmerHitCoverageDistributionTP, kmerHitWCovDistributionTP, coverageProportionDistributionTP,evidenceProportionDistributionTP,indelsKbpDistributionTP};
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
			Collections.sort(edges,(e1,e2)-> calculateScoreForEdgeLabeling(e2)-calculateScoreForEdgeLabeling(e1));
			int maxScore= calculateScoreForEdgeLabeling(edges.get(0));
			Collections.sort(edges,(e1,e2)-> e2.getOverlap()-e1.getOverlap());
			
			int maxOverlap = edges.get(0).getOverlap();
			if(vertices[i].getSequenceIndex()==debugSeq) System.out.println("AddEdges3. Vertex "+vertices[i]+" max overlap: "+maxOverlap+" max score: "+maxScore);
			for(AssemblyEdge edge:edges) {
				if(edge.isSameSequenceEdge()) continue;
				if(vertices[i].getSequenceIndex()==debugSeq) System.out.println("AddEdges3. Vertex "+vertices[i]+" next edge: "+edge);
				if(calculateScoreForEdgeLabeling(edge)<0.9*maxScore) continue;
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
			AssemblyEdge bestScore=null;
			for(AssemblyEdge edge:edges) {
				if(edge.isSameSequenceEdge()) continue;
				Integer j = verticesMap.get(edge.getConnectingVertex(vertex).getUniqueNumber());
				if(j==null || used[j]) continue;
				if(bestOverlap==null || bestOverlap.getOverlap()<edge.getOverlap()) bestOverlap = edge;
				if(bestScore==null || calculateScoreForEdgeLabeling(bestScore)<calculateScoreForEdgeLabeling(edge)) bestScore = edge;
			}
			if(vertex.getSequenceIndex()==debugSeq) System.out.println("AddEdges2. Trying to connect vertex "+vertex+" edges: "+edges.size()+" max overlap: "+bestOverlap+" max score: "+bestScore);
			if(bestOverlap == bestScore) continue;
			int diffOverlap = Math.abs(bestOverlap.getOverlap()-bestScore.getOverlap());
			int diffScore = Math.abs(calculateScoreForEdgeLabeling(bestOverlap)-calculateScoreForEdgeLabeling(bestScore));
			if(vertex.getSequenceIndex()==debugSeq) System.out.println("AddEdges2. Trying to connect vertex "+vertex+" diff overlap: "+diffOverlap+" diff score: "+diffScore);
			if(diffOverlap>0.05*bestOverlap.getOverlap()) continue;
			if(diffScore>0.05*calculateScoreForEdgeLabeling(bestScore)) continue;
			
			AssemblyEdge thirdOverlap=null;
			AssemblyEdge thirdScore=null;
			for(AssemblyEdge edge:edges) {
				if(edge.isSameSequenceEdge()) continue;
				if(edge == bestOverlap) continue;
				if(edge == bestScore) continue;
				Integer j = verticesMap.get(edge.getConnectingVertex(vertex).getUniqueNumber());
				if(j==null || used[j]) continue;
				if(thirdOverlap==null || thirdOverlap.getOverlap()<edge.getOverlap()) thirdOverlap = edge;
				if(thirdScore==null || calculateScoreForEdgeLabeling(thirdScore)<calculateScoreForEdgeLabeling(edge)) thirdScore = edge;
			}
			if(vertex.getSequenceIndex()==debugSeq) System.out.println("AddEdges2. Trying to connect vertex "+vertex+" third overlap: "+thirdOverlap+" third score: "+thirdScore);
			if(thirdOverlap!=null && thirdOverlap.getOverlap()>0.9*bestScore.getOverlap()) continue;
			if(thirdScore!=null && calculateScoreForEdgeLabeling(thirdScore)>0.9*calculateScoreForEdgeLabeling(bestOverlap)) continue;
			
			AssemblyVertex vBestOv = bestOverlap.getConnectingVertex(vertex);
			AssemblyVertex vBestScore = bestScore.getConnectingVertex(vertex);
			
			int j = verticesMap.get(vBestOv.getUniqueNumber());
			int k = verticesMap.get(vBestScore.getUniqueNumber());
			if(paths.get(j/2).size()>1) continue;
			if(paths.get(k/2).size()>1) continue;
			
			AssemblyVertex vBestOv2 = graph.getSameSequenceEdge(vBestOv).getConnectingVertex(vBestOv);
			AssemblyVertex vBestScore2 = graph.getSameSequenceEdge(vBestScore).getConnectingVertex(vBestScore);
			AssemblyEdge edgeOv2 = graph.getEdge(vBestOv2, vBestScore);
			AssemblyEdge edgeScore2 = graph.getEdge(vBestScore2, vBestOv);
			
			if(isBestEdge(graph, vBestOv, bestOverlap) && edgeOv2!=null && edgeScore2==null && edgeOv2.getOverlap()>0.8*bestScore.getOverlap()) {
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
					System.out.println("Competing edge: "+bestScore);
				}
			} else if (isBestEdge(graph, vBestScore, bestScore) && edgeOv2==null && edgeScore2!=null && edgeScore2.getOverlap()>0.8*bestScore.getOverlap()) {
				//Best score follows transitivity
				int c1 = clusters[i];
				int c2 = clusters[j];
				int c3 = clusters[k];
				
				if(c1!=c2 && c1!=c3 && c2!=c3) {
					pathEdges.add(bestScore);
					pathEdges.add(edgeScore2);
					used[i] = used[j] = used[k] = true;
					int k2 = verticesMap.get(vBestScore2.getUniqueNumber());
					used[k2] = true;
					for(int c=0;c<n;c++) {
						if(clusters[c]==c2 || clusters[c]==c3) clusters[c] = c1;
					}
					System.out.println("Added transitivity edges with best score: "+bestScore+" "+edgeScore2);
					System.out.println("Competing edge: "+bestOverlap);
				}
			}
		}
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
	private void addConnectingEdges(AssemblyGraph graph, List<LinkedList<AssemblyEdge>> paths, List<AssemblyEdge> pathEdges, Distribution[] edgesStats) {
		AssemblyVertex [] vertices = extractEndVertices(paths);
		List<AssemblyEdgePathEnd> candidateEdges = new ArrayList<AssemblyEdgePathEnd>();
		for(int i=0;i<vertices.length;i++) {
			for(int j=i+1;j<vertices.length;j++) {
				if(i%2==0 && j==i+1) continue;
				AssemblyEdge edge = graph.getEdge(vertices[i], vertices[j]);
				if(edge !=null) candidateEdges.add(new AssemblyEdgePathEnd(edge, i, j));
			}
		}
		Collections.sort(candidateEdges,(e1,e2)->calculateCost(e1.getEdgeAssemblyGraph(),edgesStats)-calculateCost(e2.getEdgeAssemblyGraph(),edgesStats));
		List<AssemblyEdgePathEnd> pathEdgesByVertexId = selectEdgesToMergePaths(candidateEdges,vertices.length, edgesStats);
		for(AssemblyEdgePathEnd edgeP:pathEdgesByVertexId) pathEdges.add(edgeP.getEdgeAssemblyGraph());
	}

	private int calculateCost(AssemblyEdge edge, Distribution[] edgesStats) {
		Distribution overlapTP = edgesStats[0];
		Distribution covTP = edgesStats[1];
		Distribution wCovTP = edgesStats[2];
		Distribution wCovPropTP = edgesStats[3];
		Distribution evPropTP = edgesStats[4];
		Distribution indelsKbpTP = edgesStats[5];
		//double prop = (double)edge.getCoverageSharedKmers()/edge.getOverlap();
		//int cost = (int)Math.round(1000*(1.5-prop));
		//return cost;
		NormalDistribution noTP = new NormalDistribution(overlapTP.getAverage(),overlapTP.getVariance()+1);
		NormalDistribution ncTP = new NormalDistribution(covTP.getAverage(),covTP.getVariance()+1);
		NormalDistribution nwcTP = new NormalDistribution(wCovTP.getAverage(),wCovTP.getVariance()+1);
		NormalDistribution nwcpTP = new NormalDistribution(wCovPropTP.getAverage(),wCovPropTP.getVariance()+0.0001);
		NormalDistribution evpTP = new NormalDistribution(evPropTP.getAverage(),evPropTP.getVariance()+0.0001);
		NormalDistribution ikbpTP = new NormalDistribution(indelsKbpTP.getAverage(),indelsKbpTP.getVariance()+1);
		//NormalDistribution niTP = new NormalDistribution(200,40000);
		double pValueOTP = noTP.cumulative(edge.getOverlap());
		//if(pValueOTP>0.5) pValueOTP = 1- pValueOTP;
		int cost1 = PhredScoreHelper.calculatePhredScore(pValueOTP);
		double pValueCTP = ncTP.cumulative(edge.getCoverageSharedKmers());
		int cost2 = PhredScoreHelper.calculatePhredScore(pValueCTP);
		double pValueWCTP = nwcTP.cumulative(edge.getWeightedCoverageSharedKmers());
		int cost3 = PhredScoreHelper.calculatePhredScore(pValueWCTP);
		double pValueWCPTP = nwcpTP.cumulative((double)edge.getWeightedCoverageSharedKmers()/edge.getOverlap());
		int cost4 = PhredScoreHelper.calculatePhredScore(pValueWCPTP);
		double pValueEvProp = evpTP.cumulative(edge.calculateEvidenceProportion());
		int cost5 = PhredScoreHelper.calculatePhredScore(pValueEvProp);
		double pValueIKBP = 1-ikbpTP.cumulative(edge.getIndelsPerKbp());
		int cost6 = PhredScoreHelper.calculatePhredScore(pValueIKBP);
		double costD = 0;
		costD+=cost1;
		//cost += cost2;
		costD += cost3;
		costD += cost5;
		if(useIndels) costD += cost6;
		
		costD*=1000;
		int cost = (int)Math.min(1000000000, costD);
		
		//cost+= (int) (1000000*(1-pValueOTP)*(1-pValueCTP));
		cost+= (int) (1000*(1-pValueOTP)*(1-pValueWCTP));

		if( logEdge(edge)) System.out.println("CalculateCost. Pvalues "+pValueOTP+" "+pValueCTP+" "+pValueWCTP+" "+pValueWCPTP+" "+pValueEvProp+" "+pValueIKBP+" costs: "+cost1+" "+cost2+" "+cost3+" "+cost4+" "+cost5+" "+cost6+" cost: " +cost+ " Edge: "+edge);
		return cost;
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
