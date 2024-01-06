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
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerSearchResultsCompressedTable;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.ShortKmerCodesTable;

/**
 * @author Jorge Duitama
 */
public class MinimizersUngappedSearchHitsClustersFinder implements UngappedSearchHitsClustersFinder {

	private Logger log = Logger.getLogger(MinimizersUngappedSearchHitsClustersFinder.class.getName());
	
	
	private int minRawHits = 10;
	private double minProportionReadLength = 0.01;
	
	private ReferenceGenome genome;
	private ShortKmerCodesTable kmerCodesTable;
	
	
	
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
		kmerCodesTable = new ShortKmerCodesTable(kmerLength, windowLength);
		
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
		kmerCodesTable.setLimitHitsPerSequence(1000);
		log.info("Filling kmer codes");
		
		//TODO: Improve parallelization
		//ThreadPoolManager poolTable = new ThreadPoolManager(numThreads, n);
		ThreadPoolManager poolTable = new ThreadPoolManager(1, n);
		poolTable.setSecondsPerTask(longestSequenceLengthMbp);
		for (int i=0;i<n;i++) {
			final int seqId = i;
			try {
				poolTable.queueTask(()->kmerCodesTable.addSequence(seqId, genome.getSequenceCharacters(seqId)));
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
		return buildHitClusters(query, false, false);
	}
	public List<UngappedSearchHitsCluster> buildHitClusters (CharSequence query, boolean extensiveKmersSearch, boolean filterClusters) {
		int queryLength = query.length();
		Map<Integer,List<UngappedSearchHit>> hitsByReference;
		if(extensiveKmersSearch) {
			Map<Integer, Long> codes = KmersExtractor.extractDNAKmerCodes(query.toString(), kmerCodesTable.getKmerLength() , 0, query.length());
			KmerSearchResultsCompressedTable results = kmerCodesTable.matchCompressed(-1, query.length(), codes, -1);
			hitsByReference = results.getAllHits();	
		}
		else hitsByReference = kmerCodesTable.match(-1,query);
		List<UngappedSearchHitsCluster> clusters = new ArrayList<UngappedSearchHitsCluster>();
		double minRawHitsSize = Math.max(minRawHits, minProportionReadLength*query.length());
		for (int sequenceIdx:hitsByReference.keySet()) {
			int sequenceLength = genome.getSequenceByIndex(sequenceIdx).getLength();
			List<UngappedSearchHit> totalHitsSubject = hitsByReference.get(sequenceIdx);
			//System.out.println("Reference id: "+sequenceIdx+" total raw hits: "+totalHitsSubject.size());
			//for(UngappedSearchHit hit: totalHitsSubject) if(hit.getQueryIdx()>21000 && hit.getQueryIdx()<27000) System.out.println("Next hit. "+hit.getQueryIdx()+" "+hit.getQuery()+" "+hit.getStart()+" "+hit.getWeight()+" "+(hit.getStart()-hit.getQueryIdx()));
			//for(UngappedSearchHit hit: totalHitsSubject) if(query.length()==11805 && hit.getStart()>0 && hit.getStart()<1000000) System.out.println("Next hit. "+hit.getQueryIdx()+" "+hit.getQuery()+" "+hit.getStart()+" "+hit.getWeight()+" "+(hit.getStart()-hit.getQueryIdx()));
			Collections.sort(totalHitsSubject, (h1,h2)->h1.getSubjectStart()-h2.getSubjectStart());
			List<UngappedSearchHit> rawClusterKmers = new ArrayList<UngappedSearchHit>();
			UngappedSearchHitsCluster cluster = null;
			for(UngappedSearchHit hit:totalHitsSubject) {
				if(hit.getWeight()<0.01) {
					//System.out.println("Hit with low weight. Pos: "+hit.getQueryStart()+" weight: "+hit.getWeight());
					continue;
				}
				if(cluster==null) {
					cluster = new UngappedSearchHitsCluster(queryLength, sequenceIdx, sequenceLength, hit);
					//System.out.println("Created first cluster. Evidence limits: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hit length: "+hit.getHitLength()+" weight: "+hit.getWeight());
				} else if (!cluster.addKmerHit(hit, 0)) {
					if (rawClusterKmers.size()>=minRawHitsSize) {
						List<UngappedSearchHitsCluster> regionClusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(queryLength, sequenceIdx, sequenceLength, rawClusterKmers, 0);
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
				List<UngappedSearchHitsCluster> regionClusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(queryLength, sequenceIdx, sequenceLength, rawClusterKmers, 0);
				//System.out.println("Qlen: "+query.length()+" next raw cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+cluster.getNumDifferentKmers()+" subclusters "+regionClusters.size()+" cluster0 end: "+regionClusters.get(0).getSubjectPredictedEnd());
				clusters.addAll(regionClusters);
			}
			//else System.out.println("Qlen: "+query.length()+" next raw small cluster discarded "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+rawClusterKmers.size()+" limit: "+minRawHitsSize);
			//System.out.println("Reference id: "+sequenceIdx+" total clusters: "+clusters.size());
		}
		//System.out.println("Clusters before filtering: "+clusters.size()+" filter: "+filterClusters);
		if(clusters.size()>3 && filterClusters) clusters = filterClusters(clusters); 
		//System.out.println("Final number of clusters: "+clusters.size());
		return clusters;
	}
	
	

	private List<UngappedSearchHitsCluster> filterClusters(List<UngappedSearchHitsCluster> clusters) {
		List<UngappedSearchHitsCluster> filteredClusters = new ArrayList<>();
		Collections.sort(clusters, (c1,c2)-> c2.getNumDifferentKmers()-c1.getNumDifferentKmers());
		int max = clusters.get(0).getNumDifferentKmers();
		int limit = max*6/10;
		if(max==3) limit++;
		//System.out.println("Filtering clusters. Max: "+max+" limit: "+limit);
		for(UngappedSearchHitsCluster cluster:clusters) {
			//System.out.println("Filtering clusters. Next cluster count: "+cluster.getNumDifferentKmers()+" limit: "+limit);
			if(cluster.getNumDifferentKmers()<limit) break;
			filteredClusters.add(cluster);
		}
		return filteredClusters;
	}
}

