package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
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
	
	private double minProportionOverlap = 0.05;
	
	private double minProportionEvidence = 0;
	
	private int countRawHits = 0;
	private int countCompletedHits = 0;
	
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

	public int getCountRawHits() {
		return countRawHits;
	}

	public int getCountCompletedHits() {
		return countCompletedHits;
	}

	public void updateGraphWithKmerHitsMap(int queryIdx, CharSequence queryForward, Map<Integer, List<UngappedSearchHit>> hitsForward, CharSequence queryReverse, Map<Integer, List<UngappedSearchHit>> hitsReverse, double compressionFactor, int kmerLength) {
		List<UngappedSearchHit> selfHits = hitsForward.get(queryIdx);
		int selfHitsCount = (selfHits!=null)?selfHits.size():1;
		int minHits = (int) Math.max(selfHitsCount*minProportionOverlap,DEF_MIN_HITS);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" self hits: "+selfHitsCount+" min hits: "+minHits);
		//Initial selection based on raw hit counts
		List<Integer> subjectIdxsF = filterSubjectIds(queryIdx, hitsForward, minHits);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Filtered subject idxs forward: "+subjectIdxsF.size());
		List<Integer> subjectIdxsR = filterSubjectIds(queryIdx, hitsReverse, minHits);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Filtered subject idxs reverse: "+subjectIdxsR.size());
		
		if(subjectIdxsF.size()==0 && subjectIdxsR.size()==0) {
			System.out.println("Query "+queryIdx+" had zero subject ids after initial filtering. self hits: "+selfHitsCount+" min hits: "+minHits);
			return;
		}
		//Build initial clusters
		List<KmerHitsCluster> clustersForward = createClusters(queryIdx, queryForward.length(), hitsForward, subjectIdxsF);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Clusters forward: "+clustersForward.size());
		List<KmerHitsCluster> clustersReverse = createClusters(queryIdx, queryReverse.length(), hitsReverse, subjectIdxsR);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Clusters reverse: "+clustersReverse.size());
		//Combined query min coverage and percentage of kmers
		boolean completeHits = evaluateClusters (queryIdx, clustersForward, clustersReverse, minHits);
		Map<Integer, Long> queryCodesF = completeHits?KmersExtractor.extractDNAKmerCodes(queryForward.toString(), kmerLength, 0, queryForward.length()):null;
		Map<Integer, Long> queryCodesR = completeHits?KmersExtractor.extractDNAKmerCodes(queryReverse.toString(), kmerLength, 0, queryReverse.length()):null;
		processClusters(queryIdx, queryForward, false, queryCodesF, clustersForward, compressionFactor, kmerLength);
		processClusters(queryIdx, queryReverse, true, queryCodesR, clustersReverse, compressionFactor, kmerLength);
		synchronized (this) {
			if(completeHits) countCompletedHits++;
			else countRawHits++;
		}
	}

	private List<Integer> filterSubjectIds(int queryIdx, Map<Integer, List<UngappedSearchHit>> hits, int minHits) {
		List<Integer> subjectIdxs = new ArrayList<Integer>();
		for(int subjectIdx:hits.keySet()) {
			if(subjectIdx>= queryIdx) continue;
			int subjectCount = hits.get(subjectIdx).size();
			if (queryIdx == idxDebug && subjectCount>DEF_MIN_HITS) System.out.println("EdgesFinder. Query: "+queryIdx+" Subject sequence: "+subjectIdx+" hits: "+subjectCount+" min hits: "+minHits);
			if(subjectCount<minHits) continue;
			subjectIdxs.add(subjectIdx);
		}
		return subjectIdxs;
	}
	
	private List<KmerHitsCluster> createClusters(int queryIdx, int queryLength, Map<Integer, List<UngappedSearchHit>> hitsMap, List<Integer> subjectIdxs) {
		List<KmerHitsCluster> clusters = new ArrayList<KmerHitsCluster>(subjectIdxs.size());
		for(int subjectIdx:subjectIdxs) {
			List<UngappedSearchHit> hits = hitsMap.get(subjectIdx);
			QualifiedSequence subjectSequence = graph.getSequence(subjectIdx);
			int subjectLength = subjectSequence.getLength();
			List<KmerHitsCluster> subjectClusters = KmerHitsCluster.clusterRegionKmerAlns(queryLength, subjectLength, hits, 0);
			if(subjectClusters.size()==0) continue;
			Collections.sort(subjectClusters, (o1,o2)-> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
			KmerHitsCluster subjectCluster = subjectClusters.get(0);
			int numKmers = subjectCluster.getNumDifferentKmers();
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Subject idx: "+subjectIdx+" hits: "+hits.size()+" clusters: "+subjectClusters.size()+" hits best cluster: "+numKmers);
			clusters.add(subjectCluster);
		}
		return clusters;
	}
	private boolean evaluateClusters(int queryIdx, List<KmerHitsCluster> clustersForward, List<KmerHitsCluster> clustersReverse, int minHits) {
		/*List<KmerHitsCluster> allClusters = new ArrayList<KmerHitsCluster>(clustersForward.size()+clustersReverse.size());
		allClusters.addAll(clustersForward);
		allClusters.addAll(clustersReverse);
		int passCount = 0;
		int maxCount = 0;
		for(KmerHitsCluster cluster:allClusters) {
			int count = cluster.getNumDifferentKmers();
			maxCount = Math.max(maxCount, count);
			if(count>=minHits) passCount++;
		}
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Min hits: "+minHits+" number of passing clusters: "+passCount+" max count: "+maxCount);
		return passCount>3 && maxCount>=2*minHits;
		*/
		return true;
	}
	private void processClusters(int queryIdx, CharSequence query, boolean queryRC, Map<Integer, Long> queryCodes, List<KmerHitsCluster> clusters, double compressionFactor, int kmerLength) {
		int queryLength = query.length();
		for(KmerHitsCluster cluster:clusters) {
			//if(numKmers<minHits) continue;
			QualifiedSequence subjectSequence = graph.getSequence(cluster.getSubjectIdx());
			if(queryCodes!=null) cluster.completeMissingHits(subjectSequence.getCharacters().toString(),queryCodes);
			long normalizedCount = queryLength*kmerLength/subjectSequence.getLength();
			int numKmers = cluster.getNumDifferentKmers();
			normalizedCount*=numKmers;
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Query length: "+queryLength+" subject idx: "+cluster.getSubjectIdx()+" subject length: "+subjectSequence.getLength()+" kmer length: "+kmerLength+" numkmers: "+numKmers+" normalized count "+normalizedCount);
			updateGraphWithKmerCluster(queryIdx, query, queryRC, compressionFactor, cluster);
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
		embeddedEvent.setHostStartStandardDeviation((int) Math.round(cluster.getSubjectStartSD()));
		int [] alnData = MinimizersTableReadAlignmentAlgorithm.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, query.length(), cluster);
		embeddedEvent.setCoverageSharedKmers(alnData[0]);
		embeddedEvent.setWeightedCoverageSharedKmers(alnData[1]);
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
		edge.setNumSharedKmers(cluster.getNumDifferentKmers());
		edge.setOverlapStandardDeviation((int) Math.round(cluster.getPredictedOverlapSD()));
		synchronized (graph) {
			graph.addEdge(edge);
		}
		if(querySequenceId==idxDebug) System.out.println("New edge: "+edge);
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
		edge.setNumSharedKmers(cluster.getNumDifferentKmers());
		edge.setOverlapStandardDeviation((int) Math.round(cluster.getPredictedOverlapSD()));
		synchronized (graph) {
			graph.addEdge(edge);
		}
		if(querySequenceId==idxDebug) System.out.println("New edge: "+edge);
	}
}
