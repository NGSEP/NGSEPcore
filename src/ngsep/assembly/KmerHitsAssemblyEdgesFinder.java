package ngsep.assembly;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.alignments.LongReadsAligner;
import ngsep.alignments.ReadAlignment;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;

public class KmerHitsAssemblyEdgesFinder {

	private AssemblyGraph graph;
	
	private int minKmerPercentage=10;
	
	private double minProportionOverlap = 0.1;
	
	private double minProportionEvidence = 0.5;
	
	private int meanDepth = 10;
	
	
	
	
	private LongReadsAligner aligner = new LongReadsAligner();
	
	private int idxDebug = -1;
	
	
	
	public KmerHitsAssemblyEdgesFinder(AssemblyGraph graph) {
		this.graph = graph;
	}
	
	public AssemblyGraph getGraph() {
		return graph;
	}
	
	public int getMeanDepth() {
		return meanDepth;
	}

	public void setMeanDepth(int meanDepth) {
		this.meanDepth = meanDepth;
	}
	

	public int getMinKmerPercentage() {
		return minKmerPercentage;
	}

	public void setMinKmerPercentage(int minKmerPercentage) {
		this.minKmerPercentage = minKmerPercentage;
	}

	public double getMinProportionOverlap() {
		return minProportionOverlap;
	}

	public void setMinProportionOverlap(double minProportionOverlap) {
		this.minProportionOverlap = minProportionOverlap;
	}

	public void updateGraphWithKmerHitsMap(int queryIdx, CharSequence query, boolean queryRC, int selfHitsCount, Map<Integer, List<UngappedSearchHit>> hitsBySubjectIdx) {
		Set<Integer> subjectIdxs = new HashSet<Integer>();
		for(int subjectIdx:hitsBySubjectIdx.keySet()) {
			int subjectCount = hitsBySubjectIdx.get(subjectIdx).size();
			if(subjectIdx< queryIdx) {
				if (queryIdx == idxDebug && subjectCount>10) System.out.println("EdgesFinder. Query: "+queryIdx+" "+queryRC+" Subject sequence: "+subjectIdx+" hits: "+subjectCount+" self hits: "+selfHitsCount);
				subjectIdxs.add(subjectIdx);
			}
		}
	
		//COmbined query min coverage and percentage of kmers
		int minCount = (int) (minProportionOverlap*minKmerPercentage*selfHitsCount/100);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" "+queryRC+" Subject sequences: "+subjectIdxs.size());
		for(int subjectIdx:subjectIdxs) {
			List<UngappedSearchHit> hits = hitsBySubjectIdx.get(subjectIdx);
			if(hits.size()<minCount) continue;
			List<KmerHitsCluster> subjectClusters = LongReadsAligner.clusterRegionKmerAlns(queryIdx, query, hits, 0);
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" "+queryRC+" Subject idx: "+subjectIdx+" hits: "+hits.size()+" clusters: "+subjectClusters.size());
			updateGraphWithKmerClusters(queryIdx, query.length(), queryRC, selfHitsCount, subjectClusters);
		}
	}
	public void updateGraphWithKmerClusters(int querySequenceId, int queryLength,  boolean queryRC, int selfHitsCount, List<KmerHitsCluster> clusteredKmerAlns) {
		//Process clusters
		Collections.sort(clusteredKmerAlns, (o1,o2)-> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
		for (int i=0;i<clusteredKmerAlns.size() && i<10;i++) {
			KmerHitsCluster cluster = clusteredKmerAlns.get(i);
			int targetSeqIdx = cluster.getSequenceIdx();
			int targetLength = graph.getSequenceLength(targetSeqIdx);
			cluster.summarize(meanDepth);
			double overlap = estimateOverlap (cluster, queryLength, targetLength);
			double regionSelfCount = overlap*selfHitsCount/queryLength;
			double pct = 100.0*cluster.getNumDifferentKmers()/regionSelfCount;
			int queryEvidenceLength = cluster.getQueryEnd()-cluster.getQueryStart();
			int subjectEvidenceLength = cluster.getSubjectEvidenceEnd() - cluster.getSubjectEvidenceStart();
			if(querySequenceId==idxDebug) System.out.println("EdgesFinder. Processing cluster. query length "+queryLength+" target length: "+targetLength+" QueryStart: "+cluster.getQueryStart()+" query end: "+cluster.getQueryEnd()+" overlap "+overlap+" Subject: "+cluster.getSequenceIdx()+" predicted limits: "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd()+" evidence limits: "+cluster.getSubjectEvidenceStart()+" - "+cluster.getSubjectEvidenceEnd() +" plain count: "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount()+" pct: "+pct+" coverage: "+cluster.getQueryCoverage());
			if(overlap < minProportionOverlap*queryLength || overlap < minProportionOverlap*targetLength) continue;
			if(queryEvidenceLength < minProportionEvidence*overlap || subjectEvidenceLength < minProportionEvidence*overlap) continue;
			if(pct<minKmerPercentage) continue;
			synchronized (graph) {
				processAlignment(graph, querySequenceId, queryRC, cluster);
			}
		}
	}

	private double estimateOverlap(KmerHitsCluster cluster, int queryLength, int targetLength) {
		double overlap = Math.min(targetLength,cluster.getSubjectPredictedEnd())-Math.max(0, cluster.getSubjectPredictedStart());
		return Math.min(overlap, queryLength);
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
		int startTarget = cluster.getSubjectPredictedStart();
		int endTarget = cluster.getSubjectPredictedEnd();
		if(querySequenceId==idxDebug) System.out.println("Processing cluster. Query: "+querySequenceId+" length: "+queryLength+ " target: "+cluster.getSequenceIdx()+" length: "+targetLength);
		if(startTarget>=0 && endTarget<=targetLength) {
			addEmbedded(graph, querySequenceId, queryRC, cluster);
		} else if (startTarget>=0) {
			addQueryAfterTargetEdge(graph, querySequenceId, queryRC, cluster);
		} else if (endTarget<=targetLength) {
			addQueryBeforeTargetEdge(graph, querySequenceId, queryRC, cluster);
		} else {
			ReadAlignment aln = aligner.alignRead(graph.getSequence(targetSeqIdx), graph.getSequence(querySequenceId), 0, targetLength, "Subject", 0.5);
			if(aln!=null) {
				int firstQueryMatch = 0;
				int newStart = targetLength;
				int newEnd = -1;
				for(int i=0;i<queryLength;i++) {
					int subjectPos = aln.getReferencePosition(i);
					if(subjectPos>0) {
						firstQueryMatch = i;
						newStart = subjectPos-1-i;
						break;
					}
				}
				for(int i=queryLength-1;i>firstQueryMatch;i--) {
					int subjectPos = aln.getReferencePosition(i);
					if(subjectPos>0) {
						newEnd= subjectPos-1+(queryLength-i);
						break;
					}
				}
				if(newStart==targetLength || newEnd==-1) return;
				cluster.setSubjectPredictedLimits(newStart, newEnd);
				if(newStart>=0 && newEnd<=targetLength) {
					addEmbedded(graph, querySequenceId, queryRC, cluster);
				} else if (newStart>=0) {
					addQueryAfterTargetEdge(graph, querySequenceId, queryRC, cluster);
				} else if (newEnd<=targetLength) {
					addQueryBeforeTargetEdge(graph, querySequenceId, queryRC, cluster);
				}
			}
		}
	}
	private void addEmbedded(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerHitsCluster cluster) {
		int startTarget = cluster.getSubjectPredictedStart();
		int targetSeqIdx = cluster.getSequenceIdx();
		AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(querySequenceId, graph.getSequence(querySequenceId), queryRC, targetSeqIdx, startTarget);
		embeddedEvent.setEvidence(cluster);
		graph.addEmbedded(embeddedEvent);
		if (querySequenceId==idxDebug) System.out.println("Query: "+querySequenceId+" embedded in "+targetSeqIdx);
	}
	private void addQueryAfterTargetEdge(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int targetSeqIdx = cluster.getSequenceIdx();
		int targetLength = graph.getSequenceLength(targetSeqIdx);
		AssemblyVertex vertexTarget = graph.getVertex(targetSeqIdx, false);
		AssemblyVertex vertexQuery;
		if(queryRC) {
			vertexQuery = graph.getVertex(querySequenceId, false); 
		} else {
			vertexQuery = graph.getVertex(querySequenceId, true);
		}
		int overlap = cluster.getPredictedOverlap();
		int cost = targetLength + queryLength - overlap;
		AssemblyEdge edge = new AssemblyEdge(vertexTarget, vertexQuery, cost, overlap);
		edge.setEvidence(cluster);
		graph.addEdge(edge);
		if(querySequenceId==idxDebug) System.out.println("Edge between target: "+vertexTarget.getUniqueNumber()+" and query "+vertexQuery.getUniqueNumber()+" overlap: "+overlap+" cost: "+cost);
	}
	private void addQueryBeforeTargetEdge(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int targetSeqIdx = cluster.getSequenceIdx();
		int targetLength = graph.getSequenceLength(targetSeqIdx);
		AssemblyVertex vertexTarget = graph.getVertex(targetSeqIdx, true);
		AssemblyVertex vertexQuery;
		if(queryRC) {
			vertexQuery = graph.getVertex(querySequenceId, true); 
		} else {
			vertexQuery = graph.getVertex(querySequenceId, false);
		}
		int overlap = cluster.getPredictedOverlap();
		int cost = targetLength + queryLength -overlap;
		AssemblyEdge edge = new AssemblyEdge(vertexQuery, vertexTarget, cost, overlap);
		edge.setEvidence(cluster);
		graph.addEdge(edge);
		if(querySequenceId==idxDebug) System.out.println("Edge between query: "+vertexQuery.getUniqueNumber()+" and target "+vertexTarget.getUniqueNumber()+" overlap: "+overlap+" cost: "+cost);
	}
}
