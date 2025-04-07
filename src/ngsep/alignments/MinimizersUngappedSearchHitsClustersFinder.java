/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.alignments;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerSearchResultsCompressedTable;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.ShortKmerCodesSampler;
import ngsep.sequences.ShortKmerCodesTable;

/**
 * @author Jorge Duitama
 */
public class MinimizersUngappedSearchHitsClustersFinder implements UngappedSearchHitsClustersFinder {

	private Logger log = Logger.getLogger(MinimizersUngappedSearchHitsClustersFinder.class.getName());
	
	
	private int minRawHits = 10;
	private double minProportionReadLength = 0.001;
	
	private ReferenceGenome genome;
	private ShortKmerCodesTable kmerCodesTable;
	private int tableKmerLength;
	
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	
	public int getMinRawHits() {
		return minRawHits;
	}
	public void setMinRawHits(int minRawHits) {
		this.minRawHits = minRawHits;
	}
	public double getMinProportionReadLength() {
		return minProportionReadLength;
	}
	public void setMinProportionReadLength(double minProportionReadLength) {
		this.minProportionReadLength = minProportionReadLength;
	}
	public void loadGenome(ReferenceGenome genome, int kmerLength, int windowLength, int numThreads) {
		loadGenome(genome, kmerLength, windowLength, numThreads,false);
	}
	public void loadGenome(ReferenceGenome genome, int kmerLength, int windowLength, int numThreads, boolean buildKmersTable) {
		this.genome = genome;
		int n = genome.getNumSequences();
		int longestSequenceLengthMbp = 1+genome.getLongestSequenceLength()/1000000;
		//KmersMapAnalyzer analyzer = new KmersMapAnalyzer(extractor.getKmersMap(), true);
		log.info("Creating kmer codes table for genome with "+n+" sequences loaded from file: "+genome.getFilename());
		ShortKmerCodesSampler sampler = new ShortKmerCodesSampler();
		sampler.setKmerLength(kmerLength);
		sampler.setWindowLength(windowLength);
		kmerCodesTable = new ShortKmerCodesTable(sampler,10*n*longestSequenceLengthMbp,true);
		tableKmerLength = kmerLength;
		if(buildKmersTable) {
			//log.info("Calculating kmers distribution");
			KmersExtractor extractor = new KmersExtractor();
			extractor.setLog(log);
			extractor.setNumThreads(numThreads);
			log.info("Calculating kmer counts from reference sequence");
			extractor.processQualifiedSequences(genome.getSequencesList());
			log.info("Calculated kmer counts from reference sequence");
			//kmerCodesTable.setKmersMap(extractor.getKmersMap());
		}
		kmerCodesTable.setLog(log);
		//TODO: Fix introduced bias
		kmerCodesTable.setMaxHitsKmerCode(1000);
		log.info("Filling kmer codes");
		
		int nt = Math.max(1, numThreads-1);
		
		ThreadPoolManager poolTable = new ThreadPoolManager(nt, nt);
		//ThreadPoolManager poolTable = new ThreadPoolManager(1, n);
		poolTable.setSecondsPerTask(longestSequenceLengthMbp);
		for (int i=0;i<n;i++) {
			final int seqId = i;
			try {
				poolTable.queueTask(()->kmerCodesTable.addSequence(seqId, genome.getSequenceCharacters(seqId).toString()));
			} catch (InterruptedException e) {
				e.printStackTrace();
				throw new RuntimeException("Concurrence error creating minimizers table",e);
			}
		}
		try {
			poolTable.terminatePool();
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException("Concurrence error creating codes table",e);
		}
		kmerCodesTable.endAddingSequences();
		//minimizersTable.calculateDistributionHits().printDistribution(System.err);
		log.info("Calculated kmer codes. Total: "+kmerCodesTable.size());
	}
	
	public ShortKmerCodesTable getKmerCodesTable() {
		return kmerCodesTable;
	}
	public void setKmerCodesTable(ReferenceGenome genome, ShortKmerCodesTable kmerCodesTable) {
		this.genome = genome;
		this.kmerCodesTable = kmerCodesTable;
	}
	
	
	@Override
	public List<UngappedSearchHitsCluster> findHitClusters(CharSequence query) {
		return buildHitClusters(query, false);
	}
	public List<UngappedSearchHitsCluster> buildHitClusters (CharSequence query, boolean extensiveKmersSearch) {
		int queryLength = query.length();
		//System.out.println("Mapping sequence with length: "+queryLength);
		UngappedSearchHitsClusterBuilder clustersBuilder = new UngappedSearchHitsClusterBuilder();
		Map<Integer,List<UngappedSearchHit>> hitsByReference;
		if(extensiveKmersSearch) {
			//By now this only affects TEs
			Map<Integer,Long> codes = KmersExtractor.extractDNAKmerCodesAsMap(query.toString(), tableKmerLength , 0, query.length());
			KmerSearchResultsCompressedTable results = kmerCodesTable.matchCompressed(-1, query.length(), codes, -1);
			hitsByReference = results.getAllHits();
			//System.out.println("Searched "+codes.size()+" kmers of length: "+tableKmerLength+". Hit sequences: "+hitsByReference.size()+" first seq hits: "+hitsByReference.get(0).size());
		}
		else hitsByReference = kmerCodesTable.match(-1,query);
		//System.out.println("Mapping sequence with length: "+queryLength+" references hits: "+hitsByReference.size());
		List<UngappedSearchHitsCluster> clusters = new ArrayList<UngappedSearchHitsCluster>();
		double minRawHitsSize = Math.max(minRawHits, minProportionReadLength*query.length());
		for (int sequenceIdx:hitsByReference.keySet()) {
			int sequenceLength = genome.getSequenceByIndex(sequenceIdx).getLength();
			List<UngappedSearchHit> totalHitsSubject = hitsByReference.get(sequenceIdx);
			//System.out.println("Reference id: "+sequenceIdx+" total raw hits subject: "+totalHitsSubject.size());
			//for(UngappedSearchHit hit: totalHitsSubject) if(hit.getQueryStart()>20000 && hit.getQueryStart()<27000) System.out.println("Next hit. "+hit.getQueryIdx()+" "+hit.getQuery()+" "+hit.getStart()+" "+hit.getWeight()+" "+(hit.getStart()-hit.getQueryIdx()));
			//for(UngappedSearchHit hit: totalHitsSubject) if(hit.getSubjectStart()>50000 && hit.getSubjectStart()<70000) System.out.println("Next hit. "+hit.getQueryStart()+" "+hit.getSubjectStart()+" "+hit.getWeight()+" "+(hit.getSubjectStart()-hit.getQueryStart()));
			Collections.sort(totalHitsSubject, (h1,h2)->h1.getSubjectStart()-h2.getSubjectStart());
			List<UngappedSearchHit> rawClusterKmers = new ArrayList<UngappedSearchHit>();
			UngappedSearchHitsCluster cluster = null;
			for(UngappedSearchHit hit:totalHitsSubject) {
				//if(hit.getWeight()<0.01) {
				//	System.out.println("Hit with low weight. Pos: "+hit.getQueryStart()+" weight: "+hit.getWeight());
				//	continue;
				//}
				if(cluster==null) {
					cluster = new UngappedSearchHitsCluster(queryLength, sequenceIdx, sequenceLength, hit);
					//System.out.println("Created first cluster. Evidence limits: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hit length: "+hit.getHitLength()+" weight: "+hit.getWeight());
				} else if (!cluster.addKmerHit(hit, 0)) {
					if (rawClusterKmers.size()>=minRawHitsSize) {
						List<UngappedSearchHitsCluster> regionClusters = clustersBuilder.clusterRegionKmerAlns(queryLength, sequenceIdx, sequenceLength, rawClusterKmers);
						//System.out.println("Qlen: "+query.length()+" next raw cluster inside "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+rawClusterKmers.size()+" subclusters "+regionClusters.size());
						clusters.addAll(regionClusters);
					} 
					//else System.out.println("Qlen: "+query.length()+" next raw small cluster discarded "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+rawClusterKmers.size()+" limit: "+minRawHitsSize);
					cluster = new UngappedSearchHitsCluster(queryLength, sequenceIdx, sequenceLength, hit);
					rawClusterKmers.clear();
				}
				rawClusterKmers.add(hit);
			}
			if(cluster!=null && rawClusterKmers.size()>=minRawHitsSize) {
				List<UngappedSearchHitsCluster> regionClusters = clustersBuilder.clusterRegionKmerAlns(queryLength, sequenceIdx, sequenceLength, rawClusterKmers);
				//System.out.println("Qlen: "+query.length()+" next raw cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+cluster.getCountKmerHitsCluster()+" subclusters "+regionClusters.size());
				clusters.addAll(regionClusters);
			}
			//else System.out.println("Qlen: "+query.length()+" next raw small cluster discarded "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+rawClusterKmers.size()+" limit: "+minRawHitsSize);
			//System.out.println("Reference id: "+sequenceIdx+" total clusters: "+clusters.size());
		}
		//System.out.println("Final number of clusters: "+clusters.size());
		return clusters;
	}
	/*public static UngappedSearchHitsCluster findBestKmersCluster (CharSequence subjectSequence, int subjectFirst, int subjectLast, CharSequence querySequence, int queryFirst, int queryLast, int kmerLength) {
		Map<Integer,Long> codesSubject = KmersExtractor.extractDNAKmerCodesAsMap(subjectSequence, kmerLength, subjectFirst, subjectLast);
		Map<Integer,Long> codesQuery = KmersExtractor.extractDNAKmerCodesAsMap(querySequence, kmerLength, queryFirst, queryLast);
		List<UngappedSearchHit> initialKmerHits = alignKmerCodes(codesSubject, codesQuery, kmerLength);
		//System.out.println("Number of kmer hits: "+initialKmerHits.size());
		if(initialKmerHits.size()==0) return null;
		List<UngappedSearchHit> filteredKmerHits = new ArrayList<>();
		for(UngappedSearchHit hit: initialKmerHits) {
			int distanceStartSubject = hit.getSubjectStart()-subjectFirst;
			int distanceStartQuery = hit.getQueryStart()-queryFirst;
			int distanceEndSubject = subjectLast-hit.getSubjectStart();
			int distanceEndQuery = queryLast-hit.getQueryStart();
			int diff1 = Math.abs(distanceStartSubject-distanceStartQuery);
			int diff2 = Math.abs(distanceEndSubject-distanceEndQuery);
			if(diff1 < 200 || diff2<200) filteredKmerHits.add(hit);
		}
		
		List<UngappedSearchHitsCluster> clusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(querySequence.length(), 0, subjectSequence.length(), filteredKmerHits, 0);
		//System.out.println("Number of clusters: "+clusters.size());
		//printClusters(clusters);
		if(clusters.size()>1) Collections.sort(clusters, (o1,o2)->o2.getNumDifferentKmers()-o1.getNumDifferentKmers());	
		else if (clusters.size()==0) return null;
		return clusters.get(0);
	}*/
	
	public static UngappedSearchHitsCluster findBestKmersCluster (int subjectLength, Map<Integer, Long> codesSubject, int queryLength, Map<Integer, Long> codesQuery, int kmerLength) {
		List<UngappedSearchHit> initialKmerHits = alignKmerCodes(codesSubject, codesQuery, kmerLength);
		//System.out.println("Number of kmer hits: "+initialKmerHits.size());
		if(initialKmerHits.size()==0) return null;
		
		List<UngappedSearchHitsCluster> clusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(queryLength, 0, subjectLength, initialKmerHits);
		if(clusters.size()>1) Collections.sort(clusters, (o1,o2)->o2.getCountKmerHitsCluster()-o1.getCountKmerHitsCluster());	
		else if (clusters.size()==0) return null;
		return clusters.get(0);
	}
	public static void printClusters(List<UngappedSearchHitsCluster> clusters) {
		System.out.println("Clusters: "+clusters.size());
		for(UngappedSearchHitsCluster cluster:clusters) {
			System.out.println("kmers: "+cluster.getCountKmerHitsCluster()+" predicted limits: "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd()+" query limits "+cluster.getQueryEvidenceStart()+"-"+cluster.getQueryEvidenceEnd());
		}
		
	}
	private static List<UngappedSearchHit> alignKmerCodes(Map<Integer, Long> codesSubject, Map<Integer, Long> codesQuery, int kmerLength) {
		Map<Long,List<Integer>> reverseSubjectMap = new HashMap<Long, List<Integer>>();
		for(Map.Entry<Integer, Long> entry:codesSubject.entrySet()) {
			List<Integer> starts = reverseSubjectMap.computeIfAbsent(entry.getValue(), v->new ArrayList<Integer>());
			starts.add(entry.getKey());
		}
		List<UngappedSearchHit> initialKmerHits = new ArrayList<UngappedSearchHit>();
		for(int i:codesQuery.keySet()) {
			Long codeRead = codesQuery.get(i);
			List<Integer> subjectPosList = reverseSubjectMap.get(codeRead);
			if(subjectPosList==null) continue;
			for(int subjectPos:subjectPosList) {
				UngappedSearchHit hit = new UngappedSearchHit(0 , subjectPos);
				hit.setQueryStart(i);
				hit.setHitLength((short)kmerLength);
				initialKmerHits.add(hit);
				//if(i<30) System.out.println("Next hit for kmer at "+i+" pos "+subjectPos); 
			}
		}
		return initialKmerHits;
	}
}

