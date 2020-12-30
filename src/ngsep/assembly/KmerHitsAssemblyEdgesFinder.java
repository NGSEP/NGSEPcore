package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.math.Distribution;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.QualifiedSequence;

public class KmerHitsAssemblyEdgesFinder {

	private AssemblyGraph graph;
	
	public static final int DEF_MIN_HITS = 50;
	
	private double minProportionOverlap = 0.05;
	
	private double minProportionEvidence = 0;
	
	private long expectedAssemblyLength = 0;
	
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
	
	public long getExpectedAssemblyLength() {
		return expectedAssemblyLength;
	}

	public void setExpectedAssemblyLength(long expectedAssemblyLength) {
		this.expectedAssemblyLength = expectedAssemblyLength;
	}

	public void updateGraphWithKmerHitsMap(int queryIdx, int queryLength, Map<Integer, Long> queryCodesF, Map<Integer, Long> queryCodesR, Map<Integer, List<UngappedSearchHit>> hitsForward, Map<Integer, List<UngappedSearchHit>> hitsReverse, double compressionFactor, int kmerLength ) {
		List<UngappedSearchHit> selfHits = hitsForward.get(queryIdx);
		int selfHitsCount = (selfHits!=null)?selfHits.size():1;
		int minHits = (int) Math.max(selfHitsCount*minProportionOverlap,DEF_MIN_HITS);
		List<KmerHitsCluster> queryClusters = KmerHitsCluster.clusterRegionKmerAlns(queryLength, queryLength, selfHits, 0);
		int kmersSelfCluster = 0;
		if(queryClusters.size()==0) {
			System.err.println("WARN: Self hits for sequence: "+queryIdx+" dd not make clusters");
		} else {
			Collections.sort(queryClusters, (o1,o2)-> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
			KmerHitsCluster cluster = queryClusters.get(0);
			AssemblyEdge edge = graph.getSameSequenceEdge(queryIdx);
			int [] alnData = MinimizersTableReadAlignmentAlgorithm.simulateAlignment(queryIdx, queryLength, queryIdx, queryLength, cluster);
			edge.setCoverageSharedKmers(alnData[0]);
			edge.setWeightedCoverageSharedKmers(alnData[1]);
			edge.setNumSharedKmers(cluster.getNumDifferentKmers());
			kmersSelfCluster = cluster.getNumDifferentKmers();
			edge.setOverlapStandardDeviation((int) Math.round(cluster.getPredictedOverlapSD()));
			edge.setRawKmerHits(cluster.getRawKmerHits());
			edge.setRawKmerHitsSubjectStartSD((int)Math.round(cluster.getRawKmerHitsSubjectStartSD()));
			minHits = (int)Math.min(minHits, minProportionOverlap*cluster.getNumDifferentKmers());
			minHits = (int) Math.max(minHits,DEF_MIN_HITS);
		}
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" self raw hits: "+selfHitsCount+" kmersSelfCluster: "+kmersSelfCluster+" min hits: "+minHits);
		//Initial selection based on raw hit counts
		List<Integer> subjectIdxsF = filterSubjectIds(queryIdx, hitsForward, minHits);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Selected subject idxs forward: "+subjectIdxsF.size());
		List<Integer> subjectIdxsR = filterSubjectIds(queryIdx, hitsReverse, minHits);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Selected subject idxs reverse: "+subjectIdxsR.size());
		
		if(subjectIdxsF.size()==0 && subjectIdxsR.size()==0) {
			//System.out.println("Query "+queryIdx+" had zero subject ids after initial filtering. self hits: "+selfHitsCount+" min hits: "+minHits);
			return;
		}
		//Build initial clusters
		List<KmerHitsCluster> clustersForward = createClusters(queryIdx, queryLength, hitsForward, subjectIdxsF);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Clusters forward: "+clustersForward.size());
		List<KmerHitsCluster> clustersReverse = createClusters(queryIdx, queryLength, hitsReverse, subjectIdxsR);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Clusters reverse: "+clustersReverse.size());
		//Combined query min coverage and percentage of kmers
		int minClusterSize = calculateMinimumClusterSize(queryIdx, clustersForward, clustersReverse, minHits);
		processClusters(queryIdx, queryLength, false, queryCodesF, clustersForward, compressionFactor, kmerLength, minClusterSize);
		processClusters(queryIdx, queryLength, true, queryCodesR, clustersReverse, compressionFactor, kmerLength, minClusterSize);
		/*synchronized (this) {
			if(completeHits) countCompletedHits++;
			else countRawHits++;
		}*/
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
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" name "+graph.getSequence(queryIdx).getName()+" Subject: "+subjectSequence.getName()+" hits: "+hits.size()+" clusters: "+subjectClusters.size()+" hits best cluster: "+numKmers);
			clusters.add(subjectCluster);
		}
		return clusters;
	}
	private int calculateMinimumClusterSize(int queryIdx, List<KmerHitsCluster> clustersForward, List<KmerHitsCluster> clustersReverse, int minHits) {
		List<KmerHitsCluster> allClusters = new ArrayList<KmerHitsCluster>(clustersForward.size()+clustersReverse.size());
		allClusters.addAll(clustersForward);
		allClusters.addAll(clustersReverse);
		int maxCount = 0;
		int passCount = 0;
		for(KmerHitsCluster cluster:allClusters) {
			int count = cluster.getNumDifferentKmers();
			maxCount = Math.max(maxCount, count);
			if(count >=minHits) passCount++;
		}
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Min hits: "+minHits+" number of passing clusters: "+passCount+" max count: "+maxCount);
		return Math.max(minHits, maxCount/5);
		
		//return true;
	}
	//By now not used
	private Map<Integer, Map<Integer,Long>> codesCache = new HashMap<Integer, Map<Integer,Long>>();
	//private long processedAssemblyLength = 0;
	private Distribution candidateClusterSizesDist = new Distribution(0, 100, 1);
	private void processClusters(int queryIdx, int queryLength, boolean queryRC, Map<Integer, Long> queryCodes, List<KmerHitsCluster> clusters, double compressionFactor, int kmerLength, int minClusterSize) {
		for(KmerHitsCluster cluster:clusters) {
			int numKmers = cluster.getNumDifferentKmers();
			QualifiedSequence subjectSequence = graph.getSequence(cluster.getSubjectIdx());
			int subjectLength = subjectSequence.getLength();
			//if (queryIdx == idxDebug || clusters.size()>80) System.out.println("EdgesFinder. Query: "+queryIdx+" name "+graph.getSequence(queryIdx).getName()+" length: "+queryLength+" Subject: "+subjectSequence.getName()+" idx: "+cluster.getSubjectIdx()+" length: "+subjectLength+" numkmers: "+numKmers+" minCluster size "+minClusterSize);
			if(numKmers<minClusterSize) continue;
			/*Map<Integer,Long> subjectCodes = codesCache.get(cluster.getSubjectIdx());
			
			if(subjectCodes==null) {
				subjectCodes = KmersExtractor.extractDNAKmerCodes(subjectSequence.getCharacters().toString(), kmerLength, cluster.getSubjectEvidenceStart(), cluster.getSubjectEvidenceEnd());
			} else {
				Map<Integer,Long> selectedCodes = new HashMap<Integer, Long>();
				for(int i:subjectCodes.keySet()) {
					if(i>=cluster.getSubjectEvidenceStart() && i<=cluster.getSubjectEvidenceEnd()) selectedCodes.put(i, subjectCodes.get(i));
				}
				subjectCodes = selectedCodes;
			}
			
			if(queryCodes!=null) cluster.completeMissingHits(subjectCodes,queryCodes);*/
			long normalizedCount = queryLength*kmerLength/subjectLength;
			
			normalizedCount*=numKmers;
			//if (queryIdx == idxDebug || clusters.size()>80) System.out.println("EdgesFinder. Query: "+queryIdx+" name "+graph.getSequence(queryIdx).getName()+" Query length: "+queryLength+" Subject: "+subjectSequence.getName()+" subject idx: "+cluster.getSubjectIdx()+" subject length: "+subjectLength+" kmer length: "+kmerLength+" numkmers: "+numKmers+" normalized count "+normalizedCount);
			updateGraphWithKmerCluster(queryIdx, queryLength, queryRC, compressionFactor, cluster);
		}
		/*candidateClusterSizesDist.processDatapoint(clusters.size());
		if((queryIdx+1)%1000==0 && !queryRC) {
			candidateClusterSizesDist.printDistributionInt(System.out);
		}*/
		/*if(!queryRC && processedAssemblyLength<expectedAssemblyLength) {
			synchronized (codesCache) {
				codesCache.put(queryIdx, queryCodes);
				processedAssemblyLength+=queryLength;
			}
		}*/
	}
	private void updateGraphWithKmerCluster(int querySequenceId, int queryLength,  boolean queryRC, double compressionFactor, KmerHitsCluster cluster) {
		//Process cluster
		if(passFilters(querySequenceId, queryLength, cluster)) {
			processCluster(querySequenceId, queryLength, queryRC, compressionFactor, cluster);
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

	private void processCluster(int querySequenceId, int queryLength, boolean queryRC, double compressionFactor, KmerHitsCluster cluster) {
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		//Zero based limits
		int startSubject = cluster.getSubjectPredictedStart();
		int endSubject = cluster.getSubjectPredictedEnd();
		if(querySequenceId==idxDebug) System.out.println("Processing cluster. Query: "+querySequenceId+" length: "+queryLength+ " subject: "+cluster.getSubjectIdx()+" length: "+subjectLength);
		if(startSubject>=0 && endSubject<=subjectLength) {
			addEmbedded(querySequenceId, queryLength, queryRC, cluster);
		} else if (startSubject>=0) {
			addQueryAfterSubjectEdge(querySequenceId, queryRC, compressionFactor, cluster);
		} else if (endSubject<=subjectLength) {
			addQueryBeforeSubjectEdge(querySequenceId, queryRC, compressionFactor, cluster);
		} else {
			// Similar sequences. Add possible embedded
			addEmbedded(querySequenceId, queryLength, queryRC, cluster);
		}
	}
	private void addEmbedded(int querySequenceId, int queryLength, boolean queryRC, KmerHitsCluster cluster) {
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
		int [] alnData = MinimizersTableReadAlignmentAlgorithm.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, queryLength, cluster);
		embeddedEvent.setCoverageSharedKmers(alnData[0]);
		embeddedEvent.setWeightedCoverageSharedKmers(alnData[1]);
		embeddedEvent.setRawKmerHits(cluster.getRawKmerHits());
		embeddedEvent.setRawKmerHitsSubjectStartSD((int)Math.round(cluster.getRawKmerHitsSubjectStartSD()));
		synchronized (graph) {
			graph.addEmbedded(embeddedEvent);
		}
		if (querySequenceId==idxDebug) System.out.println("Query: "+querySequenceId+" embedded in "+subjectSeqIdx);
	}
	private void addQueryAfterSubjectEdge(int querySequenceId, boolean queryRC, double compressionFactor, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		AssemblyVertex vertexSubject = graph.getVertex(subjectSeqIdx, false);
		AssemblyVertex vertexQuery = graph.getVertex(querySequenceId, !queryRC);
		int overlap = (int) ((double)cluster.getPredictedOverlap()/compressionFactor);
		AssemblyEdge edge = new AssemblyEdge(vertexSubject, vertexQuery, overlap);
		edge.setAverageOverlap(cluster.getAveragePredictedOverlap());
		edge.setMedianOverlap(cluster.getMedianPredictedOverlap());
		edge.setFromLimitsOverlap(cluster.getFromLimitsPredictedOverlap());
		//ReadAlignment aln = aligner.buildCompleteAlignment(subjectSeqIdx, graph.getSequence(subjectSeqIdx).getCharacters(), query, cluster);
		//int mismatches = overlap;
		//if(aln!=null) mismatches = aln.getNumMismatches();
		int [] alnData = MinimizersTableReadAlignmentAlgorithm.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, queryLength, cluster);
		edge.setCoverageSharedKmers(alnData[0]);
		edge.setWeightedCoverageSharedKmers(alnData[1]);
		edge.setNumSharedKmers(cluster.getNumDifferentKmers());
		edge.setOverlapStandardDeviation((int) Math.round(cluster.getPredictedOverlapSD()));
		edge.setRawKmerHits(cluster.getRawKmerHits());
		edge.setRawKmerHitsSubjectStartSD((int)Math.round(cluster.getRawKmerHitsSubjectStartSD()));
		synchronized (graph) {
			graph.addEdge(edge);
		}
		if(querySequenceId==idxDebug) System.out.println("New edge: "+edge);
	}
	private void addQueryBeforeSubjectEdge(int querySequenceId, boolean queryRC, double compressionFactor, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		AssemblyVertex vertexSubject = graph.getVertex(subjectSeqIdx, true);
		AssemblyVertex vertexQuery = graph.getVertex(querySequenceId, queryRC);
		int overlap = (int) ((double)cluster.getPredictedOverlap()/compressionFactor);
		AssemblyEdge edge = new AssemblyEdge(vertexQuery, vertexSubject, overlap);
		edge.setAverageOverlap(cluster.getAveragePredictedOverlap());
		edge.setMedianOverlap(cluster.getMedianPredictedOverlap());
		edge.setFromLimitsOverlap(cluster.getFromLimitsPredictedOverlap());
		//ReadAlignment aln = aligner.buildCompleteAlignment(subjectSeqIdx, graph.getSequence(subjectSeqIdx).getCharacters(), query, cluster);
		//int mismatches = overlap;
		//if(aln!=null) mismatches = aln.getNumMismatches();
		int [] alnData = MinimizersTableReadAlignmentAlgorithm.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, queryLength, cluster);
		edge.setCoverageSharedKmers(alnData[0]);
		edge.setWeightedCoverageSharedKmers(alnData[1]);
		edge.setNumSharedKmers(cluster.getNumDifferentKmers());
		edge.setOverlapStandardDeviation((int) Math.round(cluster.getPredictedOverlapSD()));
		edge.setRawKmerHits(cluster.getRawKmerHits());
		edge.setRawKmerHitsSubjectStartSD((int)Math.round(cluster.getRawKmerHitsSubjectStartSD()));
		synchronized (graph) {
			graph.addEdge(edge);
		}
		if(querySequenceId==idxDebug) System.out.println("New edge: "+edge);
	}
}
