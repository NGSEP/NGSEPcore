package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.PairwiseAlignerDynamicKmers;
import ngsep.alignments.UngappedSearchHitsCluster;
import ngsep.alignments.UngappedSearchHitsClusterBuilder;
import ngsep.alignments.LongReadsUngappedSearchHitsClusterAligner;
import ngsep.alignments.ReadAlignment;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerSearchResultsCompressedTable;
import ngsep.sequences.QualifiedSequence;

public class KmerHitsAssemblyEdgesFinder {
	
	private Logger log = Logger.getLogger(KmerHitsAssemblyEdgesFinder.class.getName());

	private AssemblyGraph graph;
	
	private UngappedSearchHitsClusterBuilder clustersBuilder = new UngappedSearchHitsClusterBuilder();
	
	public static final int DEF_MIN_HITS = 50;
	
	private double minProportionOverlap = 0;
	
	private double minProportionEvidence = 0;
	
	private int idxDebug = -1;
	
	private boolean extensiveSearch = false;
	
	private boolean completeAlignment = false;
	
	
	
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
	
	public boolean isExtensiveSearch() {
		return extensiveSearch;
	}

	public void setExtensiveSearch(boolean extensiveSearch) {
		this.extensiveSearch = extensiveSearch;
	}

	public boolean isCompleteAlignment() {
		return completeAlignment;
	}

	public void setCompleteAlignment(boolean completeAlignment) {
		this.completeAlignment = completeAlignment;
	}

	public List<AssemblySequencesRelationship> inferRelationshipsFromKmerHits(int queryIdx, CharSequence queryF, CharSequence queryR, KmerSearchResultsCompressedTable hitsForward, KmerSearchResultsCompressedTable hitsReverse) {
		List<AssemblySequencesRelationship> relationships = new ArrayList<AssemblySequencesRelationship>();
		int queryLength = queryF.length();
		
		int numSelfHits = hitsForward.getKmerHitCount(queryIdx); 
		int distinctKmersCount = hitsForward.countDistinctKmerHits(queryIdx);
		if(distinctKmersCount==0) distinctKmersCount = queryLength/10;
		
		int minHits = (int) Math.max(distinctKmersCount/8,DEF_MIN_HITS);
		//minHits = calculateMinimumHitsFromTotal(queryIdx, hitsForward, hitsReverse, minHits);
		if(!extensiveSearch) minHits*=2;
		if (queryIdx == idxDebug) log.info("EdgesFinder. Query: "+queryIdx+" self hits: "+numSelfHits+" self distinct kmers: "+distinctKmersCount+" min hits: "+minHits+" Memory: "+calculateMemoryGbp());
		//Initial selection based on raw hit counts
		List<Integer> subjectIdxsF = filterAndSortSubjectIds(queryIdx, hitsForward, minHits);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Selected subject idxs forward: "+subjectIdxsF.size()+" Memory: "+calculateMemoryGbp());
		List<Integer> subjectIdxsR = filterAndSortSubjectIds(queryIdx, hitsReverse, minHits);
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Selected subject idxs reverse: "+subjectIdxsR.size()+" Memory: "+calculateMemoryGbp());
		//List<AssemblySequencesRelationship> relationships = new ArrayList<AssemblySequencesRelationship>(subjectIdxsF.size()+subjectIdxsR.size());
		if(subjectIdxsF.size()==0 && subjectIdxsR.size()==0) {
			//System.out.println("Query "+queryIdx+" had zero subject ids after initial filtering. self hits: "+selfHitsCount+" min hits: "+minHits);
			return relationships;
		}
		
		if(extensiveSearch) {
			//Build initial clusters
			List<UngappedSearchHitsCluster> clustersForward = createClusters(queryIdx, queryLength, hitsForward, subjectIdxsF);
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Clusters forward: "+clustersForward.size()+" Memory: "+calculateMemoryGbp());
			List<UngappedSearchHitsCluster> clustersReverse = createClusters(queryIdx, queryLength, hitsReverse, subjectIdxsR);
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" Clusters reverse: "+clustersReverse.size()+" Memory: "+calculateMemoryGbp());
			//Combined query min coverage and percentage of kmers
			int minClusterSize = calculateMinimumClusterSize(queryIdx, clustersForward, clustersReverse, minHits);
			for(UngappedSearchHitsCluster cluster:clustersForward) processCluster(queryIdx, queryF, false, cluster, minClusterSize, relationships);
			for(UngappedSearchHitsCluster cluster:clustersReverse) processCluster(queryIdx, queryR, true, cluster, minClusterSize, relationships);
		} else {
			int i=0;
			int minClusterSize = minHits;
			while(i<subjectIdxsF.size() && i<subjectIdxsR.size() && i<10) {
				int subjectIdxF = subjectIdxsF.get(i);
				List<UngappedSearchHitsCluster> subjectClustersF = createClusters(queryIdx, queryLength, subjectIdxF, hitsForward.getHits(subjectIdxF));
				for(UngappedSearchHitsCluster clusterF:subjectClustersF) {
					if (processCluster(queryIdx, queryF, false, clusterF, minClusterSize, relationships)) return relationships;
					minClusterSize = Math.max(minClusterSize, clusterF.getNumDifferentKmers()/2);
				}
				int subjectIdxR = subjectIdxsR.get(i);
				List<UngappedSearchHitsCluster> subjectClustersR = createClusters(queryIdx, queryLength, subjectIdxR, hitsReverse.getHits(subjectIdxR));
				for(UngappedSearchHitsCluster clusterR:subjectClustersR) {
					if (processCluster(queryIdx, queryR, true, clusterR, minClusterSize, relationships)) return relationships;
					minClusterSize = Math.max(minClusterSize, clusterR.getNumDifferentKmers()/2);
				}
				i++;
			}
			if(i==10) return relationships;
			for(;i<subjectIdxsF.size();i++) {
				int subjectIdxF = subjectIdxsF.get(i);
				List<UngappedSearchHitsCluster> subjectClusters = createClusters(queryIdx, queryLength, subjectIdxF, hitsForward.getHits(subjectIdxF));
				for(UngappedSearchHitsCluster cluster:subjectClusters) if(processCluster(queryIdx, queryF, false, cluster, minClusterSize, relationships)) return relationships;
			}
			for(;i<subjectIdxsR.size();i++) {
				int subjectIdxR = subjectIdxsR.get(i);
				List<UngappedSearchHitsCluster> subjectClusters = createClusters(queryIdx, queryLength, subjectIdxR, hitsReverse.getHits(subjectIdxR));
				for(UngappedSearchHitsCluster cluster:subjectClusters) if(processCluster(queryIdx, queryR, true, cluster, minClusterSize, relationships)) return relationships;
			}
		}
		if (queryIdx == idxDebug) log.info("EdgesFinder. Query: "+queryIdx+" relationships: "+relationships.size()+" Memory: "+calculateMemoryGbp());
		return relationships;
	}

	public static long calculateMemoryGbp() {
		Runtime runtime = Runtime.getRuntime();
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		return usedMemory;
	}

	private List<Integer> filterAndSortSubjectIds(int queryIdx, KmerSearchResultsCompressedTable results, int minHits) {
		Map<Integer,Integer> hitCounts = results.getKmerHitCountsBySubjectId();
		List<Integer> subjectIdxs = new ArrayList<Integer>();
		Map<Integer,Integer> numQueryKmersSubject = new HashMap<Integer, Integer>();
		for(Map.Entry<Integer, Integer> entry:hitCounts.entrySet()) {
			int subjectIdx = entry.getKey();
			int subjectRawCount = entry.getValue(); 
			if(subjectIdx>= queryIdx) continue;
			int numQueryKmers = results.countDistinctKmerHits(subjectIdx);
			if (queryIdx == idxDebug && subjectRawCount>DEF_MIN_HITS) System.out.println("EdgesFinder. Query: "+queryIdx+" Subject sequence: "+subjectIdx+" hits: "+subjectRawCount+" distinct kmers: "+numQueryKmers+" min hits: "+minHits);
			if(subjectRawCount<minHits) continue;
			if(numQueryKmers<minHits) continue;
			numQueryKmersSubject.put(subjectIdx, numQueryKmers);
			subjectIdxs.add(subjectIdx);
		}
		Collections.sort(subjectIdxs, (s1,s2)->numQueryKmersSubject.get(s2)-numQueryKmersSubject.get(s1));
		return subjectIdxs;
	}
	
	private List<UngappedSearchHitsCluster> createClusters(int queryIdx, int queryLength, KmerSearchResultsCompressedTable results, List<Integer> subjectIdxs) {
		List<UngappedSearchHitsCluster> clusters = new ArrayList<UngappedSearchHitsCluster>(subjectIdxs.size());
		for(int subjectIdx:subjectIdxs) {
			List<UngappedSearchHit> hits = results.getHits(subjectIdx);
			List<UngappedSearchHitsCluster> subjectClusters = createClusters(queryIdx, queryLength, subjectIdx, hits);
			clusters.addAll(subjectClusters);
		}
		return clusters;
	}
	private List<UngappedSearchHitsCluster> createClusters (int queryIdx, int queryLength, int subjectIdx, List<UngappedSearchHit> hits) {
		QualifiedSequence subjectSequence = graph.getSequence(subjectIdx);
		int subjectLength = subjectSequence.getLength();
		List<UngappedSearchHitsCluster> subjectClusters = clustersBuilder.clusterRegionKmerAlns(queryLength, subjectIdx, subjectLength, hits, 0);
		List<UngappedSearchHitsCluster> answer = new ArrayList<UngappedSearchHitsCluster>(2);
		if(subjectClusters.size()==0) return answer;
		Collections.sort(subjectClusters, (o1,o2)-> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
		UngappedSearchHitsCluster bestCluster = subjectClusters.get(0);
		answer.add(bestCluster);
		int numKmers = bestCluster.getNumDifferentKmers();
		boolean queryBefore = bestCluster.getSubjectPredictedStart()<0;
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Query: "+queryIdx+" "+graph.getSequence(queryIdx).getName()+" Subject: "+subjectIdx+" "+subjectSequence.getName()+" hits: "+hits.size()+" subject clusters: "+subjectClusters.size()+" hits best subject cluster: "+numKmers+" subjectId in cluster: "+bestCluster.getSubjectIdx());
		if(subjectClusters.size()==1) return answer;
		UngappedSearchHitsCluster secondCluster = subjectClusters.get(1);
		if(secondCluster.getNumDifferentKmers()>0.8*numKmers && queryBefore != secondCluster.getSubjectPredictedStart()<0) {
			answer.add(secondCluster);
			if (queryIdx == idxDebug) System.out.println("EdgesFinder. Adding second cluster. Query: "+queryIdx+" "+graph.getSequence(queryIdx).getName()+" Subject: "+subjectIdx+" "+subjectSequence.getName()+" hits second subject cluster: "+secondCluster.getNumDifferentKmers());
		}
		return answer;
	}
	
	private int calculateMinimumClusterSize(int queryIdx, List<UngappedSearchHitsCluster> clustersForward, List<UngappedSearchHitsCluster> clustersReverse, int minHits) {
		List<UngappedSearchHitsCluster> allClusters = new ArrayList<UngappedSearchHitsCluster>(clustersForward.size()+clustersReverse.size());
		allClusters.addAll(clustersForward);
		allClusters.addAll(clustersReverse);
		int maxCount = 0;
		int passCount = 0;
		for(UngappedSearchHitsCluster cluster:allClusters) {
			int count = cluster.getNumDifferentKmers();
			maxCount = Math.max(maxCount, count);
			if(count >=minHits) passCount++;
		}
		if (queryIdx == idxDebug) System.out.println("EdgesFinder. Min hits: "+minHits+" number of passing clusters: "+passCount+" max count: "+maxCount);
		return Math.max(minHits, maxCount/5);
	}
	public void printSubjectHits(int subjectIdx, List<UngappedSearchHit> subjectHits) {
		for(UngappedSearchHit hit:subjectHits) {
			System.out.println(hit.getQueryStart()+" "+subjectIdx+":"+hit.getSubjectStart());
		}
		
	}

	private boolean processCluster(int querySequenceId, CharSequence query, boolean queryRC, UngappedSearchHitsCluster cluster, int minClusterSize, List<AssemblySequencesRelationship> relationships) {
		int queryLength = query.length();
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		if(querySequenceId==idxDebug) System.out.println("Processing cluster. Query: "+querySequenceId+" length: "+queryLength+ " subject: "+cluster.getSubjectIdx()+" length: "+subjectLength+" subject predicted "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd());
		if(!passFilters(querySequenceId, queryLength, minClusterSize, cluster)) {
			cluster.disposeHits();
			return false;
		}
		
		//Zero based limits
		int startSubject = cluster.getSubjectPredictedStart();
		int endSubject = cluster.getSubjectPredictedEnd();
		
		if(startSubject>=0 && endSubject<=subjectLength) {
			boolean answer = addEmbedded(querySequenceId, query, queryRC, cluster, relationships);
			cluster.disposeHits();
			return answer;
		} else if (startSubject>=0) {
			addQueryAfterSubjectEdge(querySequenceId, query, queryRC, cluster, relationships);
		} else if (endSubject<=subjectLength) {
			addQueryBeforeSubjectEdge(querySequenceId, query, queryRC, cluster, relationships);
		} else {
			// Similar sequences. Add possible embedded
			addEmbedded(querySequenceId, query, queryRC, cluster, relationships);
		}
		cluster.disposeHits();
		return false;
	}
	private boolean passFilters (int querySequenceId, int queryLength, int minClusterSize, UngappedSearchHitsCluster cluster) {
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		double overlap = cluster.getPredictedOverlap();
		int queryEvidenceLength = cluster.getQueryEvidenceEnd()-cluster.getQueryEvidenceStart();
		int subjectEvidenceLength = cluster.getSubjectEvidenceEnd() - cluster.getSubjectEvidenceStart();
		if(querySequenceId==idxDebug) System.out.println("EdgesFinder. Evaluating cluster. kmers: "+cluster.getNumDifferentKmers()+" min "+minClusterSize+ " subject: "+cluster.getSubjectIdx()+" length: "+subjectLength+" subject predicted "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd());
		if(cluster.getNumDifferentKmers()<minClusterSize) return false;
		if(querySequenceId==idxDebug) System.out.println("EdgesFinder. Evaluating cluster. qlen "+queryLength+" QPred: "+cluster.getQueryPredictedStart()+" - "+cluster.getQueryPredictedEnd()+" QEv: "+cluster.getQueryEvidenceStart()+" - "+cluster.getQueryEvidenceEnd()+" subject len: "+subjectLength+" Subject: "+cluster.getSubjectIdx()+" sPred: "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd()+" sEv: "+cluster.getSubjectEvidenceStart()+" - "+cluster.getSubjectEvidenceEnd()+" overlap1 "+overlap+" overlap2: "+cluster.getPredictedOverlap() +" plain count: "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount()+" min: "+(minProportionOverlap*queryLength));
		if(overlap < minProportionOverlap*queryLength) return false;
		if(overlap < minProportionOverlap*subjectLength) return false;
		if(queryEvidenceLength < minProportionEvidence*overlap) return false;
		if(subjectEvidenceLength < minProportionEvidence*overlap) return false;
		
		return true;
	}
	private boolean addEmbedded(int querySequenceId, CharSequence query, boolean queryRC, UngappedSearchHitsCluster cluster, List<AssemblySequencesRelationship> relationships) {
		int startSubject = cluster.getSubjectPredictedStart();
		int endSubject = cluster.getSubjectPredictedEnd();
		int subjectSeqIdx = cluster.getSubjectIdx();
		int queryLength = query.length();
		QualifiedSequence subjectSequence = graph.getSequence(subjectSeqIdx);
		int subjectLength = subjectSequence.getLength();
		double proportionEvidence = cluster.getSubjectEvidenceEnd()-cluster.getSubjectEvidenceStart();
		proportionEvidence/=queryLength;
		
		ReadAlignment aln = null;
		int [] simulatedAlnData = null;
		if(querySequenceId==idxDebug) System.out.println("Candidate embedded: "+subjectSeqIdx+" "+subjectSequence.getName()+" propEv "+proportionEvidence);
		if(completeAlignment) {
			if(querySequenceId==idxDebug) System.out.println("Performing complete alignment for embedded candidate: "+querySequenceId+" "+graph.getSequence(querySequenceId).getName()+" subject: "+subjectSeqIdx+" "+subjectSequence.getName());
			LongReadsUngappedSearchHitsClusterAligner aligner = new LongReadsUngappedSearchHitsClusterAligner(LongReadsUngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
			aln = aligner.buildAlignment(query, graph.getSequence(subjectSeqIdx).getCharacters(),  cluster);
			if(aln==null) {
				System.err.println("Alignment could not be performed for embedded candidate. Query: "+graph.getSequence(querySequenceId).getName()+" subject: "+subjectSequence.getName());
				return false;
			}
			if(querySequenceId==idxDebug) System.out.println("Alignment limits subject: "+aln.getFirst()+" - "+aln.getLast()+" query: "+aln.getAlignedReadPosition(aln.getFirst())+" - "+aln.getAlignedReadPosition(aln.getLast())+" CIGAR: "+aln.getCigarString());
			int skipStart = aln.getSoftClipStart();
			int skipEnd = aln.getSoftClipEnd();
			if(aln.getFirst()-skipStart<0) {
				if(querySequenceId==idxDebug) System.out.println("The alignment start minus the skipped base pairs goes beyond the subject sequence. Aln first: "+aln.getFirst()+" skip: "+skipStart);
				addQueryBeforeSubjectEdge(querySequenceId, query, queryRC, cluster ,relationships);
				return false;
			}
			if(aln.getLast()+skipEnd>subjectLength) {
				if(querySequenceId==idxDebug) System.out.println("The alignment end plus the skipped base pairs goes beyond the subject sequence. Aln end: "+aln.getLast()+" skip: "+skipEnd);
				addQueryAfterSubjectEdge(querySequenceId, query, queryRC, cluster, relationships);
				return false;
			}
		} else {
			simulatedAlnData =  PairwiseAlignerDynamicKmers.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, queryLength, cluster);
		}
		AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(querySequenceId, graph.getSequence(querySequenceId), queryRC, subjectSeqIdx, graph.getSequence(subjectSeqIdx), startSubject, endSubject);
		embeddedEvent.setNumSharedKmers(cluster.getNumDifferentKmers());
		embeddedEvent.setHostStartStandardDeviation((int) Math.round(cluster.getSubjectStartSD()));
		embeddedEvent.setRawKmerHits(cluster.getRawKmerHits());
		embeddedEvent.setRawKmerHitsSubjectStartSD((int)Math.round(cluster.getRawKmerHitsSubjectStartSD()));
		if(aln!=null) {
			embeddedEvent.setHostEvidenceStart(aln.getFirst()-1);
			embeddedEvent.setHostEvidenceEnd(aln.getLast());
			embeddedEvent.setSequenceEvidenceStart(aln.getAlignedReadPosition(aln.getFirst()));
			embeddedEvent.setSequenceEvidenceEnd(aln.getAlignedReadPosition(aln.getLast())+1);
			embeddedEvent.setNumMismatches(aln.getNumMismatches());
			embeddedEvent.setCoverageSharedKmers(aln.getCoverageSharedKmers());
			embeddedEvent.setWeightedCoverageSharedKmers(aln.getWeightedCoverageSharedKmers());
			embeddedEvent.setNumIndels(aln.getTotalLengthIndelCalls());
			embeddedEvent.setAliginment(aln);
		} else {
			embeddedEvent.setHostEvidenceStart(cluster.getSubjectEvidenceStart());
			embeddedEvent.setHostEvidenceEnd(cluster.getSubjectEvidenceEnd());
			embeddedEvent.setSequenceEvidenceStart(cluster.getQueryEvidenceStart());
			embeddedEvent.setSequenceEvidenceEnd(cluster.getQueryEvidenceEnd());
			embeddedEvent.setCoverageSharedKmers(simulatedAlnData[0]);
			embeddedEvent.setWeightedCoverageSharedKmers(simulatedAlnData[1]);
			embeddedEvent.setNumIndels(simulatedAlnData[2]);
			if(embeddedEvent.getEvidenceProportion()>0.9 && embeddedEvent.getCoverageSharedKmers()>0.5*queryLength && simulatedAlnData[2]>50 ) {
				QualifiedSequence subject = graph.getSequence(subjectSeqIdx);
				simulatedAlnData = redoSimulation(subjectSeqIdx, subject.getCharacters(),Math.max(0,startSubject),Math.min(subjectLength, endSubject),querySequenceId,query,0,query.length(),simulatedAlnData);
				if(querySequenceId==idxDebug) System.out.println("Repeated simulation for embedded: "+embeddedEvent+" new indels: "+simulatedAlnData[2]);
				embeddedEvent.setNumIndels(simulatedAlnData[2]);
			}
		}
		if(queryRC) {
			int reversedStart = queryLength - embeddedEvent.getSequenceEvidenceEnd();
			int reversedEnd = queryLength - embeddedEvent.getSequenceEvidenceStart();
			embeddedEvent.setSequenceEvidenceStart(reversedStart);
			embeddedEvent.setSequenceEvidenceEnd(reversedEnd);
		}
		
		/*synchronized (graph) {
			graph.addEmbedded(embeddedEvent);
		}*/
		relationships.add(embeddedEvent);
		
		if (querySequenceId==idxDebug) System.out.println("Query: "+querySequenceId+" embedded in "+subjectSeqIdx+" proportion evidence: "+proportionEvidence);
		return embeddedEvent.getEvidenceProportion()>0.99 && embeddedEvent.getIndelsPerKbp()<5 && embeddedEvent.getWeightedCoverageSharedKmers()>0.7*queryLength;
	}
	private void addQueryAfterSubjectEdge(int querySequenceId, CharSequence query, boolean queryRC, UngappedSearchHitsCluster cluster, List<AssemblySequencesRelationship> relationships) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		AssemblyVertex vertexSubject = graph.getVertex(subjectSeqIdx, false);
		AssemblyVertex vertexQuery = graph.getVertex(querySequenceId, !queryRC);
		int overlap = (int) ((double)cluster.getPredictedOverlap());
		double proportionEvidence = cluster.getSubjectEvidenceEnd()-cluster.getSubjectEvidenceStart();
		proportionEvidence/=(overlap+1);
		
		ReadAlignment aln = null;
		int [] simulatedAlnData = null;
		if(querySequenceId==idxDebug) System.out.println("Candidate edge: "+subjectSeqIdx+" "+graph.getSequence(subjectSeqIdx).getName()+" propEv "+proportionEvidence+" overlap: "+overlap+" qlen: "+queryLength+" prop: "+overlap/queryLength);
		if(completeAlignment) {
			LongReadsUngappedSearchHitsClusterAligner aligner = new LongReadsUngappedSearchHitsClusterAligner(LongReadsUngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
			aln = aligner.buildAlignment(query,graph.getSequence(subjectSeqIdx).getCharacters(), cluster);
			if(aln==null) {
				System.err.println("Alignment could not be performed for edge candidate. Query: "+graph.getSequence(querySequenceId).getName()+" subject: "+graph.getSequence(subjectSeqIdx).getName());
				return;
			}	
			if(querySequenceId==idxDebug) System.out.println("Alignment limits subject: "+aln.getFirst()+" - "+aln.getLast()+" query: "+aln.getAlignedReadPosition(aln.getFirst())+" - "+aln.getAlignedReadPosition(aln.getLast())+" CIGAR: "+aln.getCigarString());
			if(aln.getFirst()-aln.getSoftClipStart()>0 && aln.getLast()+aln.getSoftClipEnd()<subjectLength) {
				if(querySequenceId==idxDebug) System.out.println("Sequence looks embedded after alignment. New predicted limits: "+(aln.getFirst()-aln.getSoftClipStart())+" - "+(aln.getLast()+aln.getSoftClipEnd()));
				addEmbedded(querySequenceId, query, queryRC, cluster, relationships);
			}
			overlap = subjectLength-aln.getFirst();
		} else {
			simulatedAlnData =  PairwiseAlignerDynamicKmers.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, queryLength, cluster);
		}
		
		AssemblyEdge edge = new AssemblyEdge(vertexSubject, vertexQuery, overlap);
		edge.setAverageOverlap(cluster.getAveragePredictedOverlap());
		edge.setMedianOverlap(cluster.getMedianPredictedOverlap());
		edge.setFromLimitsOverlap(cluster.getFromLimitsPredictedOverlap());
		edge.setNumSharedKmers(cluster.getNumDifferentKmers());
		edge.setOverlapStandardDeviation((int) Math.round(cluster.getPredictedOverlapSD()));
		edge.setRawKmerHits(cluster.getRawKmerHits());
		edge.setRawKmerHitsSubjectStartSD((int)Math.round(cluster.getRawKmerHitsSubjectStartSD()));
		
		if(aln!=null) {
			edge.setVertex1EvidenceStart(aln.getFirst()-1);
			edge.setVertex1EvidenceEnd(aln.getLast());
			edge.setVertex2EvidenceStart(aln.getAlignedReadPosition(aln.getFirst()));
			edge.setVertex2EvidenceEnd(aln.getAlignedReadPosition(aln.getLast())+1);
			edge.setNumMismatches(aln.getNumMismatches());
			edge.setCoverageSharedKmers(aln.getCoverageSharedKmers());
			edge.setWeightedCoverageSharedKmers(aln.getWeightedCoverageSharedKmers());
			edge.setNumIndels(aln.getTotalLengthIndelCalls());
		} else {
			edge.setVertex1EvidenceStart(cluster.getSubjectEvidenceStart());
			edge.setVertex1EvidenceEnd(cluster.getSubjectEvidenceEnd());
			edge.setVertex2EvidenceStart(cluster.getQueryEvidenceStart());
			edge.setVertex2EvidenceEnd(cluster.getQueryEvidenceEnd());
			
			edge.setCoverageSharedKmers(simulatedAlnData[0]);
			edge.setWeightedCoverageSharedKmers(simulatedAlnData[1]);
			edge.setNumIndels(simulatedAlnData[2]);
			if(edge.getEvidenceProportion()>0.9 && edge.getCoverageSharedKmers()>0.5*overlap && simulatedAlnData[2]>50  ) {
				QualifiedSequence subject = graph.getSequence(subjectSeqIdx);
				simulatedAlnData = redoSimulation(subjectSeqIdx, subject.getCharacters(),Math.max(0, subjectLength-overlap),subject.getLength(),querySequenceId,query,0,overlap,simulatedAlnData);
				if(querySequenceId==idxDebug) System.out.println("Repeated simulation for edge: "+edge+" new indels: "+simulatedAlnData[2]);
				edge.setNumIndels(simulatedAlnData[2]);
			}
		}
		if(queryRC) {
			int reversedStart = queryLength - edge.getVertex2EvidenceEnd();
			int reversedEnd = queryLength - edge.getVertex2EvidenceStart();
			edge.setVertex2EvidenceStart(reversedStart);
			edge.setVertex2EvidenceEnd(reversedEnd);
		}
		
		/*synchronized (graph) {
			graph.addEdge(edge);
		}*/
		relationships.add(edge);
		if(querySequenceId==idxDebug) System.out.println("New edge: "+edge);
	}

	private void addQueryBeforeSubjectEdge(int querySequenceId, CharSequence query, boolean queryRC, UngappedSearchHitsCluster cluster, List<AssemblySequencesRelationship> relationships) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int subjectSeqIdx = cluster.getSubjectIdx();
		int subjectLength = graph.getSequenceLength(subjectSeqIdx);
		AssemblyVertex vertexSubject = graph.getVertex(subjectSeqIdx, true);
		AssemblyVertex vertexQuery = graph.getVertex(querySequenceId, queryRC);
		int overlap = cluster.getPredictedOverlap();
		double proportionEvidence = cluster.getSubjectEvidenceEnd()-cluster.getSubjectEvidenceStart();
		proportionEvidence/=(overlap+1);

		ReadAlignment aln = null;
		int [] simulatedAlnData = null;
		if(querySequenceId==idxDebug) System.out.println("Candidate edge: "+subjectSeqIdx+" "+graph.getSequence(subjectSeqIdx).getName()+" propEv "+proportionEvidence+" overlap: "+overlap+" qlen: "+queryLength+" prop: "+overlap/queryLength);
		if(completeAlignment) {
			LongReadsUngappedSearchHitsClusterAligner aligner = new LongReadsUngappedSearchHitsClusterAligner(LongReadsUngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
			aln = aligner.buildAlignment(query, graph.getSequence(subjectSeqIdx).getCharacters(), cluster);
			if(aln==null) {
				System.err.println("Alignment could not be performed for edge candidate. Query: "+graph.getSequence(querySequenceId).getName()+" subject: "+graph.getSequence(subjectSeqIdx).getName());
				return;
			}
			if(querySequenceId==idxDebug) System.out.println("Alignment limits subject: "+aln.getFirst()+" - "+aln.getLast()+" query: "+aln.getAlignedReadPosition(aln.getFirst())+" - "+aln.getAlignedReadPosition(aln.getLast())+" CIGAR: "+aln.getCigarString());
			if(aln.getFirst()-aln.getSoftClipStart()>0 && aln.getLast()+aln.getSoftClipEnd()<subjectLength) {
				if(querySequenceId==idxDebug) System.out.println("Sequence looks embedded after alignment. New predicted limits: "+(aln.getFirst()-aln.getSoftClipStart())+" - "+(aln.getLast()+aln.getSoftClipEnd()));
				addEmbedded(querySequenceId, query, queryRC, cluster, relationships);
			}
			overlap = aln.getLast();
		} else {
			simulatedAlnData =  PairwiseAlignerDynamicKmers.simulateAlignment(subjectSeqIdx, subjectLength, querySequenceId, queryLength, cluster);
		}
		AssemblyEdge edge = new AssemblyEdge(vertexQuery, vertexSubject, overlap);
		edge.setAverageOverlap(cluster.getAveragePredictedOverlap());
		edge.setMedianOverlap(cluster.getMedianPredictedOverlap());
		edge.setFromLimitsOverlap(cluster.getFromLimitsPredictedOverlap());
		
		edge.setNumSharedKmers(cluster.getNumDifferentKmers());
		edge.setOverlapStandardDeviation((int) Math.round(cluster.getPredictedOverlapSD()));
		edge.setRawKmerHits(cluster.getRawKmerHits());
		edge.setRawKmerHitsSubjectStartSD((int)Math.round(cluster.getRawKmerHitsSubjectStartSD()));
		
		if(aln!=null) {
			edge.setVertex1EvidenceStart(aln.getAlignedReadPosition(aln.getFirst()));
			edge.setVertex1EvidenceEnd(aln.getAlignedReadPosition(aln.getLast()));
			edge.setVertex2EvidenceStart(aln.getFirst());
			edge.setVertex2EvidenceEnd(aln.getLast());
			edge.setNumMismatches(aln.getNumMismatches());
			edge.setCoverageSharedKmers(aln.getCoverageSharedKmers());
			edge.setWeightedCoverageSharedKmers(aln.getWeightedCoverageSharedKmers());
			edge.setNumIndels(aln.getTotalLengthIndelCalls());
		} else {
			edge.setVertex1EvidenceStart(cluster.getQueryEvidenceStart());
			edge.setVertex1EvidenceEnd(cluster.getQueryEvidenceEnd());
			edge.setVertex2EvidenceStart(cluster.getSubjectEvidenceStart());
			edge.setVertex2EvidenceEnd(cluster.getSubjectEvidenceEnd());
			edge.setCoverageSharedKmers(simulatedAlnData[0]);
			edge.setWeightedCoverageSharedKmers(simulatedAlnData[1]);
			edge.setNumIndels(simulatedAlnData[2]);
			if(edge.getEvidenceProportion()>0.9 && edge.getCoverageSharedKmers()>0.5*overlap && simulatedAlnData[2]>50  ) {
				QualifiedSequence subject = graph.getSequence(subjectSeqIdx);
				simulatedAlnData = redoSimulation(subjectSeqIdx, subject.getCharacters(),0,overlap,querySequenceId,query,Math.max(0, queryLength-overlap),queryLength,simulatedAlnData);
				if(querySequenceId==idxDebug) System.out.println("Repeated simulation for edge: "+edge+" new indels: "+simulatedAlnData[2]);
				edge.setNumIndels(simulatedAlnData[2]);
			}
		}
		if(queryRC) {
			int reversedStart = queryLength - edge.getVertex1EvidenceEnd();
			int reversedEnd = queryLength - edge.getVertex1EvidenceStart();
			edge.setVertex1EvidenceStart(reversedStart);
			edge.setVertex1EvidenceEnd(reversedEnd);
		}
		/*synchronized (graph) {
			graph.addEdge(edge);
		}*/
		relationships.add(edge);
		if(querySequenceId==idxDebug) System.out.println("New edge: "+edge);
	}
	private int[] redoSimulation(int subjectSeqIdx, CharSequence subjectSequence, int startSubject, int endSubject, int querySequenceId, CharSequence query, int queryStart, int queryEnd, int [] initialSimulation) {
		UngappedSearchHitsCluster cluster = PairwiseAlignerDynamicKmers.findBestKmersCluster(subjectSequence, startSubject, endSubject, query, queryStart, queryEnd,15);
		if(cluster ==null) return initialSimulation;
		int [] newSim = PairwiseAlignerDynamicKmers.simulateAlignment(subjectSeqIdx, subjectSequence.length(), querySequenceId, query.length(), cluster);
		if(newSim[2]>=initialSimulation[2]) return initialSimulation;
		return newSim;
	}
}
