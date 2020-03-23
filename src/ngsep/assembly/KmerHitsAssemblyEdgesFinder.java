package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.alignments.LongReadsAligner;
import ngsep.alignments.ReadAlignment;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;

public class KmerHitsAssemblyEdgesFinder {

	private AssemblyGraph graph;
	private int minKmerPercentage;
	
	private LongReadsAligner aligner = new LongReadsAligner();
	
	private int idxDebug = -1;
	
	
	
	public KmerHitsAssemblyEdgesFinder(AssemblyGraph graph, int minKmerPercentage) {
		this.graph = graph;
		this.minKmerPercentage = minKmerPercentage;
	}
	
	public AssemblyGraph getGraph() {
		return graph;
	}

	public int getMinKmerPercentage() {
		return minKmerPercentage;
	}



	public void updateGraphWithKmerHits(int querySequenceId, CharSequence query, boolean queryRC, List<UngappedSearchHit> kmerHitsList, int kmersCount, double averageHits) {
		// Cluster hits by target region
		int minKmers = (int) (0.5*minKmerPercentage*kmersCount/100.0);
		List<KmerHitsCluster> clusteredKmerAlns = clusterKmerHits(querySequenceId, query, kmerHitsList, Math.max(10, minKmers));
		if(querySequenceId==idxDebug) System.out.println("Query id: "+querySequenceId+" RC: "+queryRC+" kmers: "+kmersCount+" Clusters: "+clusteredKmerAlns.size());
		
		updateGraphWithKmerClusters(querySequenceId, queryRC, kmersCount, averageHits, clusteredKmerAlns);
	}

	public void updateGraphWithKmerClusters(int querySequenceId, boolean queryRC, int kmersCount, double averageHits, List<KmerHitsCluster> clusteredKmerAlns) {
		//Process clusters
		Collections.sort(clusteredKmerAlns, (o1,o2)-> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
		for (int i=0;i<clusteredKmerAlns.size() && i<10;i++) {
			KmerHitsCluster cluster = clusteredKmerAlns.get(i);
			cluster.summarize(averageHits, kmersCount);
			double pct = 100.0*cluster.getProportionKmers();
			if(querySequenceId==idxDebug) System.out.println("Processing cluster. QueryStart: "+cluster.getQueryStart()+" query end: "+cluster.getQueryEnd()+" Subject: "+cluster.getSequenceIdx()+" first: "+cluster.getFirst()+" last: "+cluster.getLast()+" plain count: "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount()+" pct: "+pct+" coverage: "+cluster.getQueryCoverage());
			if(pct<minKmerPercentage) break;
			synchronized (graph) {
				processAlignment(graph, querySequenceId, queryRC, cluster);
			}
		}
	}

	private List<KmerHitsCluster> clusterKmerHits(int querySequenceId, CharSequence query, List<UngappedSearchHit> kmerHits, int minKmers) {
		List<KmerHitsCluster> clusters = new ArrayList<>();
		Map<Integer,List<UngappedSearchHit>> hitsByTargetSequence = new HashMap<>();
		for(UngappedSearchHit kmerHit: kmerHits) {
			int targetSequenceId = kmerHit.getSequenceIdx();
			List<UngappedSearchHit> targetHits = hitsByTargetSequence.computeIfAbsent(targetSequenceId, k-> new ArrayList<UngappedSearchHit>());
			targetHits.add(kmerHit);
		}
		
		for(int targetIdx:hitsByTargetSequence.keySet()) {
			List<UngappedSearchHit> targetHits = hitsByTargetSequence.get(targetIdx);
			// TODO: choose better min coverage
			if(targetHits.size()>=minKmers) clusters.addAll(LongReadsAligner.clusterSequenceKmerAlns(querySequenceId, query, targetHits, 0));
		}
		return clusters;
	}
	
	
	public void printTargetHits(List<UngappedSearchHit> targetHits) {
		for(UngappedSearchHit hit:targetHits) {
			System.out.println(hit.getQueryIdx()+" "+hit.getSequenceIdx()+":"+hit.getStart());
		}
		
	}

	private void processAlignment(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int targetSeqIdx = cluster.getSequenceIdx();
		int targetLength = graph.getSequenceLength(targetSeqIdx);
		//Zero based limits
		int firstTarget = cluster.getFirst();
		int lastTarget = cluster.getLast();
		//System.out.println("Processing cluster. Query: "+querySequenceId+" length: "+queryLength+ " target: "+cluster.getSequenceIdx()+" length: "+targetLength);
		if(firstTarget>0 && lastTarget<=targetLength) {
			addEmbedded(graph, querySequenceId, queryRC, cluster);
		} else if (firstTarget>0) {
			addQueryAfterTargetEdge(graph, querySequenceId, queryRC, cluster);
		} else if (lastTarget<=targetLength) {
			addQueryBeforeTargetEdge(graph, querySequenceId, queryRC, cluster);
		} else {
			ReadAlignment aln = aligner.alignRead(graph.getSequence(targetSeqIdx), graph.getSequence(querySequenceId), 0, targetLength, "Subject", 0.5);
			if(aln!=null) {
				int firstQueryMatch = 0;
				int newFirst = targetLength;
				int newLast = -1;
				for(int i=0;i<queryLength;i++) {
					int subjectPos = aln.getReferencePosition(i);
					if(subjectPos>0) {
						firstQueryMatch = i;
						newFirst = subjectPos-1-i;
						break;
					}
				}
				for(int i=queryLength-1;i>=firstQueryMatch;i--) {
					int subjectPos = aln.getReferencePosition(i);
					if(subjectPos>0) {
						newLast= subjectPos-1+(queryLength-1-i);
						break;
					}
				}
				if(newFirst==targetLength || newLast==-1) return;
				cluster.setFirst(newFirst);
				cluster.setLast(newLast);
				if(newFirst>0 && newLast<=targetLength) {
					addEmbedded(graph, querySequenceId, queryRC, cluster);
				} else if (newFirst>0) {
					addQueryAfterTargetEdge(graph, querySequenceId, queryRC, cluster);
				} else if (newLast<=targetLength) {
					addQueryBeforeTargetEdge(graph, querySequenceId, queryRC, cluster);
				}
			}
		}
	}
	private void addEmbedded(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerHitsCluster cluster) {
		int firstTarget = cluster.getFirst()-1;
		int targetSeqIdx = cluster.getSequenceIdx();
		double pct = 100.0*cluster.getProportionKmers();
		double queryCoverage = cluster.getQueryCoverage();
		//TODO: improve rules for embedded sequences
		if(pct>=2*minKmerPercentage && queryCoverage>=0.5) {
			AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(querySequenceId, graph.getSequence(querySequenceId), queryRC, targetSeqIdx, firstTarget);
			embeddedEvent.setEvidence(cluster);
			graph.addEmbedded(embeddedEvent);
			if (querySequenceId==idxDebug) System.out.println("Query: "+querySequenceId+" embedded in "+targetSeqIdx);
		}
	}
	private void addQueryAfterTargetEdge(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int targetSeqIdx = cluster.getSequenceIdx();
		int targetLength = graph.getSequenceLength(targetSeqIdx);
		int queryRegionLength = cluster.getQueryEnd()-cluster.getQueryStart();
		int startTarget = cluster.getFirst()-1;
		AssemblyVertex vertexTarget = graph.getVertex(targetSeqIdx, false);
		AssemblyVertex vertexQuery;
		if(queryRC) {
			vertexQuery = graph.getVertex(querySequenceId, false); 
		} else {
			vertexQuery = graph.getVertex(querySequenceId, true);
		}
		int overlap = targetLength-startTarget;
		if (queryRegionLength > 0.5*overlap) {
			int cost = targetLength + queryLength - overlap;
			AssemblyEdge edge = new AssemblyEdge(vertexTarget, vertexQuery, cost, overlap);
			edge.setEvidence(cluster);
			graph.addEdge(edge);
		}
		//System.out.println("Edge between target: "+targetSeqIdx+" and query "+querySequenceId+" overlap: "+overlap+" weight: "+weight);
	}
	private void addQueryBeforeTargetEdge(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int targetSeqIdx = cluster.getSequenceIdx();
		int targetLength = graph.getSequenceLength(targetSeqIdx);
		int queryRegionLength = cluster.getQueryEnd()-cluster.getQueryStart();
		int endTarget = cluster.getLast();
		AssemblyVertex vertexTarget = graph.getVertex(targetSeqIdx, true);
		AssemblyVertex vertexQuery;
		if(queryRC) {
			vertexQuery = graph.getVertex(querySequenceId, true); 
		} else {
			vertexQuery = graph.getVertex(querySequenceId, false);
		}
		int overlap = endTarget;
		if (queryRegionLength > 0.5*overlap) {
			int cost = targetLength + queryLength -overlap;
			AssemblyEdge edge = new AssemblyEdge(vertexQuery, vertexTarget, cost, overlap);
			edge.setEvidence(cluster);
			graph.addEdge(edge);
		}
		// System.out.println("Edge between query: "+querySequenceId+" and target "+targetSeqIdx+" overlap: "+overlap+" weight: "+weight);
	}
}
