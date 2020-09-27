package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;

public class KmerHitsAssemblyEdgesFinder {

	private AssemblyGraph graph;
	
	public static final int DEF_MIN_HITS = 25;
	
	private double minProportionOverlap = 0.1;
	
	private double minProportionEvidence = 0;
	
	private int idxDebug = -1;
	
	
	
	public KmerHitsAssemblyEdgesFinder(AssemblyGraph graph) {
		this.graph = graph;
	}
	
	public AssemblyGraph getGraph() {
		return graph;
	}

	public double getMinProportionOverlap() {
		return minProportionOverlap;
	}

	public void setMinProportionOverlap(double minProportionOverlap) {
		this.minProportionOverlap = minProportionOverlap;
	}

	public void updateGraphWithKmerHitsMap(int queryIdx, CharSequence query, boolean queryRC, double compressionFactor, int selfHitsCount, int kmerLength, Map<Integer, List<UngappedSearchHit>> hitsBySubjectIdx) {
		Map<Integer,Integer> subjectNormalizedCounts = new HashMap<Integer,Integer>();
		List<Integer> subjectIdxs = new ArrayList<Integer>();
		//TODO: Make parameter
		int minHits = (int) Math.max(selfHitsCount*minProportionOverlap,DEF_MIN_HITS);
		//int minHits = DEF_MIN_HITS;
		for(int subjectIdx:hitsBySubjectIdx.keySet()) {
			if(subjectIdx>= queryIdx) continue;
			int subjectCount = hitsBySubjectIdx.get(subjectIdx).size();
			//Calculated over the query to avoid missing embedded sequences
			int subjectLength = graph.getSequenceLength(subjectIdx);
			int subjectNormalizedCount = subjectCount*query.length()*kmerLength/subjectLength;
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" "+queryRC+" Subject sequence: "+subjectIdx+" hits: "+subjectCount+" normalized: "+subjectNormalizedCount+" self hits: "+selfHitsCount+" min hits: "+minHits);
			if(subjectCount<minHits) continue;
			
			subjectNormalizedCounts.put(subjectIdx,subjectNormalizedCount);
			subjectIdxs.add(subjectIdx);
		}
		if(subjectNormalizedCounts.size()==0) return;
		Collections.sort(subjectIdxs,(i1,i2)-> subjectNormalizedCounts.get(i2)-subjectNormalizedCounts.get(i1));
		int maxNormalizedCount = 0;
		//Combined query min coverage and percentage of kmers
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" "+queryRC+" Subject sequences: "+subjectIdxs.size());
		Map<Integer, Long> queryCodes = KmersExtractor.extractDNAKmerCodes(query.toString(), kmerLength, 0, query.length());
		for(int i=0;i<subjectIdxs.size();i++) {
			int subjectIdx = subjectIdxs.get(i);
			List<UngappedSearchHit> hits = hitsBySubjectIdx.get(subjectIdx);
			QualifiedSequence subjectSequence = graph.getSequence(subjectIdx);
			int subjectLength = subjectSequence.getLength();
			List<KmerHitsCluster> subjectClusters = KmerHitsCluster.clusterRegionKmerAlns(query.length(), subjectLength, hits, 0);
			if(subjectClusters.size()==0) continue;
			Collections.sort(subjectClusters, (o1,o2)-> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
			KmerHitsCluster subjectCluster = subjectClusters.get(0);
			int numKmers = subjectCluster.getNumDifferentKmers();
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" "+queryRC+" Subject idx: "+subjectIdx+" hits: "+hits.size()+" clusters: "+subjectClusters.size()+" hits best cluster: "+numKmers);
			if(numKmers<minHits) continue;
			subjectCluster.completeMissingHits(subjectSequence.getCharacters().toString(),queryCodes);
			int normalizedCount = query.length()*kmerLength/subjectLength;
			numKmers = subjectCluster.getNumDifferentKmers();
			//if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" "+queryRC+" Query length: "+query.length()+" subject length: "+subjectLength+" kmer length: "+kmerLength+" numkmers: "+numKmers+" normalized count "+normalizedCount);
			normalizedCount*=numKmers;
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" "+queryRC+" Subject idx: "+subjectIdx+" new hits subject cluster: "+numKmers+" normalized count "+normalizedCount+" max: "+maxNormalizedCount);
			//if(normalizedCount<0.3*maxNormalizedCount) return;
			maxNormalizedCount = Math.max(maxNormalizedCount, normalizedCount);
			updateGraphWithKmerCluster(queryIdx, query, queryRC, compressionFactor, subjectCluster);
		}
	}
	private void updateGraphWithKmerCluster(int querySequenceId, CharSequence query,  boolean queryRC, double compressionFactor, KmerHitsCluster cluster) {
		//Process cluster
		cluster.summarize();
		if(passFilters(querySequenceId, query.length(), cluster)) {
			processCluster(querySequenceId, query, queryRC, compressionFactor, cluster);
		}
		cluster.disposeHits();
	}
	
	private boolean passFilters (int querySequenceId, int queryLength, KmerHitsCluster cluster) {
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		double overlap = cluster.getPredictedOverlap();
		int queryEvidenceLength = cluster.getQueryEvidenceEnd()-cluster.getQueryEvidenceStart();
		int subjectEvidenceLength = cluster.getSubjectEvidenceEnd() - cluster.getSubjectEvidenceStart();
		//if(querySequenceId==idxDebug) System.out.println("EdgesFinder. Evaluating cluster. qlen "+queryLength+" QPred: "+cluster.getQueryPredictedStart()+" - "+cluster.getQueryPredictedEnd()+" QEv: "+cluster.getQueryEvidenceStart()+" - "+cluster.getQueryEvidenceEnd()+" subject len: "+subjectLength+" Subject: "+cluster.getSequenceIdx()+" sPred: "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd()+" sEv: "+cluster.getSubjectEvidenceStart()+" - "+cluster.getSubjectEvidenceEnd()+" overlap1 "+overlap+" overlap2: "+cluster.getPredictedOverlap() +" plain count: "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount()+" pct: "+pct);
		if(overlap < minProportionOverlap*queryLength) return false;
		if(overlap < minProportionOverlap*subjectLength) return false;
		if(queryEvidenceLength < minProportionEvidence*overlap) return false;
		if(subjectEvidenceLength < minProportionEvidence*overlap) return false;
		
		return true;
	}

	public void printSubjectHits(List<UngappedSearchHit> subjectHits) {
		for(UngappedSearchHit hit:subjectHits) {
			System.out.println(hit.getQueryIdx()+" "+hit.getSequenceIdx()+":"+hit.getStart());
		}
		
	}

	private void processCluster(int querySequenceId, CharSequence query, boolean queryRC, double compressionFactor, KmerHitsCluster cluster) {
		int queryLength = query.length();
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		//Zero based limits
		int startSubject = cluster.getSubjectPredictedStart();
		int endSubject = cluster.getSubjectPredictedEnd();
		if(querySequenceId==idxDebug) System.out.println("Processing cluster. Query: "+querySequenceId+" length: "+queryLength+ " subject: "+cluster.getSubjectIdx()+" length: "+subjectLength);
		if(startSubject>=0 && endSubject<=subjectLength) {
			addEmbedded(querySequenceId, query, queryRC, cluster);
		} else if (startSubject>=0) {
			addQueryAfterSubjectEdge(querySequenceId, query, queryRC, compressionFactor, cluster);
		} else if (endSubject<=subjectLength) {
			addQueryBeforeSubjectEdge(querySequenceId, query, queryRC, compressionFactor, cluster);
		} else {
			// Similar sequences. Add possible embedded
			addEmbedded(querySequenceId, query, queryRC, cluster);
		}
	}
	private void addEmbedded(int querySequenceId, CharSequence query, boolean queryRC, KmerHitsCluster cluster) {
		int startSubject = cluster.getSubjectPredictedStart();
		int endSubject = cluster.getSubjectPredictedEnd();
		int subjectSeqIdx = cluster.getSubjectIdx();
		QualifiedSequence subjectSequence = graph.getSequence(subjectSeqIdx);
		int subjectLength = subjectSequence.getLength();
		AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(querySequenceId, graph.getSequence(querySequenceId), queryRC, subjectSeqIdx, startSubject, endSubject);
		embeddedEvent.setHostEvidenceStart(cluster.getSubjectEvidenceStart());
		embeddedEvent.setHostEvidenceEnd(cluster.getSubjectEvidenceEnd());
		embeddedEvent.setNumSharedKmers(cluster.getNumDifferentKmers());
		int [] alnData = MinimizersTableReadAlignmentAlgorithm.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, query.length(), cluster);
		embeddedEvent.setCoverageSharedKmers(alnData[0]);
		embeddedEvent.setWeightedCoverageSharedKmers(alnData[1]);
		embeddedEvent.setMismatches(alnData[2]);
		synchronized (graph) {
			graph.addEmbedded(embeddedEvent);
		}
		if (querySequenceId==idxDebug) System.out.println("Query: "+querySequenceId+" embedded in "+subjectSeqIdx);
	}
	private void addQueryAfterSubjectEdge(int querySequenceId, CharSequence query, boolean queryRC, double compressionFactor, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		AssemblyVertex vertexSubject = graph.getVertex(subjectSeqIdx, false);
		AssemblyVertex vertexQuery = graph.getVertex(querySequenceId, !queryRC);
		int overlap = (int) ((double)cluster.getPredictedOverlap()/compressionFactor);
		AssemblyEdge edge = new AssemblyEdge(vertexSubject, vertexQuery, overlap);
		//ReadAlignment aln = aligner.buildCompleteAlignment(subjectSeqIdx, graph.getSequence(subjectSeqIdx).getCharacters(), query, cluster);
		//int mismatches = overlap;
		//if(aln!=null) mismatches = aln.getNumMismatches();
		int [] alnData = MinimizersTableReadAlignmentAlgorithm.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, queryLength, cluster);
		edge.setCoverageSharedKmers(alnData[0]);
		edge.setWeightedCoverageSharedKmers(alnData[1]);
		edge.setMismatches(alnData[2]);
		edge.setNumSharedKmers(cluster.getNumDifferentKmers());
		edge.setOverlapStandardDeviation(cluster.getPredictedOverlapSD());
		synchronized (graph) {
			graph.addEdge(edge);
		}
		if(querySequenceId==idxDebug) System.out.println("Edge between subject: "+vertexSubject.getUniqueNumber()+" and query "+vertexQuery.getUniqueNumber()+" overlap: "+overlap+" mismatches: "+edge.getMismatches()+" cost: "+edge.getCost());
	}
	private void addQueryBeforeSubjectEdge(int querySequenceId, CharSequence query, boolean queryRC, double compressionFactor, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		AssemblyVertex vertexSubject = graph.getVertex(subjectSeqIdx, true);
		AssemblyVertex vertexQuery = graph.getVertex(querySequenceId, queryRC);
		int overlap = (int) ((double)cluster.getPredictedOverlap()/compressionFactor);
		AssemblyEdge edge = new AssemblyEdge(vertexQuery, vertexSubject, overlap);
		//ReadAlignment aln = aligner.buildCompleteAlignment(subjectSeqIdx, graph.getSequence(subjectSeqIdx).getCharacters(), query, cluster);
		//int mismatches = overlap;
		//if(aln!=null) mismatches = aln.getNumMismatches();
		int [] alnData = MinimizersTableReadAlignmentAlgorithm.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, queryLength, cluster);
		edge.setCoverageSharedKmers(alnData[0]);
		edge.setWeightedCoverageSharedKmers(alnData[1]);
		edge.setMismatches(alnData[2]);
		edge.setNumSharedKmers(cluster.getNumDifferentKmers());
		edge.setOverlapStandardDeviation(cluster.getPredictedOverlapSD());
		synchronized (graph) {
			graph.addEdge(edge);
		}
		if(querySequenceId==idxDebug) System.out.println("Edge between query: "+vertexQuery.getUniqueNumber()+" and subject "+vertexSubject.getUniqueNumber()+" overlap: "+overlap+" mismatches: "+edge.getMismatches()+" cost: "+edge.getCost());
	}
}
