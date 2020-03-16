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
package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.FMIndexUngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;

/**
 * @author Jorge Duitama
 * @author Juan Camilo Bojaca
 * @author David Guevara
 */
public class GraphBuilderFMIndex implements GraphBuilder {
	private Logger log = Logger.getLogger(GraphBuilderFMIndex.class.getName());
	private final static int TALLY_DISTANCE = 100;
	private final static int SUFFIX_FRACTION = 20;

	private int kmerLength;
	private int kmerOffset;
	private int minKmerPercentage;
	
	private static int idxDebug = -1;
	
	

	public GraphBuilderFMIndex(int kmerLength, int kmerOffset, int minKmerPercentage) {
		this.kmerLength = kmerLength;
		this.kmerOffset = kmerOffset;
		this.minKmerPercentage = minKmerPercentage;
	}
	
	/**
	 * @return the log
	 */
	public Logger getLog() {
		return log;
	}
	
	/**
	 * @param log the log to set
	 */
	public void setLog(Logger log) {
		this.log = log;
	}

	@Override
	public AssemblyGraph buildAssemblyGraph(List<CharSequence> sequences) {
		
		AssemblyGraph graph = new AssemblyGraph(sequences);
		log.info("Created graph vertices. Edges: "+graph.getEdges().size());
		// Create FM-Index
		FMIndex fmIndex = new FMIndex();
		fmIndex.loadUnnamedSequences(sequences, TALLY_DISTANCE, SUFFIX_FRACTION);
		log.info("Created FM-Index");
	
		for (int seqId = 0; seqId < sequences.size(); seqId++) {
			String seq = sequences.get(seqId).toString();
			updateGraph(graph, seqId, seq, false, fmIndex);
			String complement = DNAMaskedSequence.getReverseComplement(seq).toString();
			updateGraph(graph, seqId, complement, true, fmIndex);
			graph.filterEdgesAndEmbedded (seqId);
			if(seqId == idxDebug) System.out.println("Edges start: "+graph.getEdges(graph.getVertex(seqId, true)).size()+" edges end: "+graph.getEdges(graph.getVertex(seqId, false)).size()+" Embedded: "+graph.getEmbeddedBySequenceId(seqId));
			if ((seqId+1)%100==0) log.info("Processed "+(seqId+1) +" sequences. Number of edges: "+graph.getEdges().size()+ " Embedded: "+graph.getEmbeddedCount());
		}
		log.info("Built graph. Edges: "+graph.getEdges().size()+" Embedded: "+graph.getEmbeddedCount()+" Prunning embedded sequences");
		graph.pruneEmbeddedSequences();
		log.info("Prunned graph. Edges: "+graph.getEdges().size());
		
		return graph;
	}


	private void updateGraph(AssemblyGraph graph, int querySequenceId, String query, boolean queryRC, FMIndex fmIndex) {
		Map<Integer,CharSequence> kmersMap = KmersExtractor.extractKmersAsMap(query, kmerLength, kmerOffset, true, true, true);
		
		
		//Search kmers using the FM index
		if(kmersMap.size()==0) return;
		int kmersCount=0;
		double averageHits = 0;
		List<FMIndexUngappedSearchHit> initialKmerHits = new ArrayList<>();
		for (int start:kmersMap.keySet()) {
			String kmer = kmersMap.get(start).toString();
			//List<FMIndexUngappedSearchHit> kmerHits=fmIndex.exactSearch(kmer,0,querySequenceId-1);
			List<FMIndexUngappedSearchHit> kmerHits=fmIndex.exactSearch(kmer);
			kmersCount++;
			averageHits+=kmerHits.size();
			//if(querySequenceId==52) System.out.println("Query: "+querySequenceId+" complement: "+queryRC+" Found "+kmerHits.size()+" hits for kmer: "+kmer);
			for(FMIndexUngappedSearchHit hit:kmerHits) {
				//if(querySequenceId==52) System.out.println("Kmer start: "+hit.getStart()+" Next alignment: "+aln.getSequenceIndex()+": "+aln.getFirst()+"-"+aln.getLast()+" rc: "+aln.isNegativeStrand());
				if(hit.getSequenceIdx()>=querySequenceId) continue;
				hit.setQueryIdx(start);
				initialKmerHits.add(hit);
			}
		}
		if(kmersCount==0) return;
		averageHits/=kmersCount;
		
		
		//Filter repetitive kmer hits
		if(querySequenceId==idxDebug) System.out.println("Query: "+querySequenceId+" complement: "+queryRC+" Average hits "+averageHits);
		List<FMIndexUngappedSearchHit> filteredHits = new ArrayList<FMIndexUngappedSearchHit>();
		Set<Integer> kmerStarts = new HashSet<Integer>();
		for(FMIndexUngappedSearchHit hit:initialKmerHits) {
			if(hit.getTotalHitsQuery()>5*averageHits) continue;
			filteredHits.add(hit);
			kmerStarts.add(hit.getQueryIdx());
		}
		
		// Cluster hits by target region
		kmersCount = kmerStarts.size();
		if(kmersCount==0) return;
		int minKmers = (int) (0.5*minKmerPercentage*kmersCount/100);
		List<KmerHitsCluster> clusteredKmerAlns = clusterKmerHits(querySequenceId, query, filteredHits, Math.max(10, minKmers));
		if(querySequenceId==idxDebug) System.out.println("Query id: "+querySequenceId+" RC: "+queryRC+" kmers: "+kmersCount+" Clusters: "+clusteredKmerAlns.size());
		
		//Process clusters
		Collections.sort(clusteredKmerAlns, (o1,o2)-> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
		double firstScoreCluster = 0;
		double firstQueryCoverage = 0;
		for (int i=0;i<clusteredKmerAlns.size() && i<10;i++) {
			KmerHitsCluster cluster = clusteredKmerAlns.get(i);
			cluster.summarize(averageHits, kmersCount);
			double pct = 100.0*cluster.getProportionKmers();
			if(querySequenceId==idxDebug) System.out.println("Processing cluster. Subject: "+cluster.getSequenceIdx()+" first: "+cluster.getFirst()+" last: "+cluster.getLast()+" plain count: "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount()+" pct: "+pct+" coverage: "+cluster.getQueryCoverage());
			if(pct<minKmerPercentage) break;
			if(i==0) {
				firstScoreCluster = pct;
				firstQueryCoverage = cluster.getQueryCoverage();
			}
			else if (3.0*pct<firstScoreCluster || 2.0*cluster.getQueryCoverage()<firstQueryCoverage) break;
			processAlignment(graph, querySequenceId, queryRC, cluster);
		}
	}

	private List<KmerHitsCluster> clusterKmerHits(int querySequenceId, CharSequence query, List<FMIndexUngappedSearchHit> kmerHits, int minKmers) {
		List<KmerHitsCluster> clusters = new ArrayList<>();
		Map<Integer,List<FMIndexUngappedSearchHit>> hitsByTargetSequence = new HashMap<>();
		for(FMIndexUngappedSearchHit kmerHit: kmerHits) {
			int targetSequenceId = kmerHit.getSequenceIdx();
			List<FMIndexUngappedSearchHit> targetHits = hitsByTargetSequence.computeIfAbsent(targetSequenceId, k-> new ArrayList<FMIndexUngappedSearchHit>());
			targetHits.add(kmerHit);
		}
		
		for(int targetIdx:hitsByTargetSequence.keySet()) {
			List<FMIndexUngappedSearchHit> targetHits = hitsByTargetSequence.get(targetIdx);
			if(targetHits.size()>=minKmers) clusters.addAll(clusterSequenceKmerAlns(querySequenceId, query, targetHits));
		}
		return clusters;
	}
	public static List<KmerHitsCluster> clusterSequenceKmerAlns(int querySequenceId, CharSequence query, List<FMIndexUngappedSearchHit> sequenceKmerHits) {
		List<KmerHitsCluster> answer = new ArrayList<>();
		
		KmerHitsCluster uniqueCluster = new KmerHitsCluster(query, sequenceKmerHits);
		if(querySequenceId==idxDebug) System.out.println("Hits to cluster: "+sequenceKmerHits.size()+" target: "+uniqueCluster.getSequenceIdx()+" first: "+uniqueCluster.getFirst()+" last: "+uniqueCluster.getLast()+" kmers: "+uniqueCluster.getNumDifferentKmers());
		answer.add(uniqueCluster);
		if(uniqueCluster.getNumDifferentKmers()>0.8*sequenceKmerHits.size()) return answer;
		//Cluster remaining hits
		List<FMIndexUngappedSearchHit> remainingHits = new ArrayList<FMIndexUngappedSearchHit>();
		for(FMIndexUngappedSearchHit hit:sequenceKmerHits) {
			if (hit!=uniqueCluster.getKmerHit(hit.getQueryIdx())) remainingHits.add(hit);
		}
		answer.add(new KmerHitsCluster(query, remainingHits));
		return answer;
	}
	
	public void printTargetHits(List<FMIndexUngappedSearchHit> targetHits) {
		for(FMIndexUngappedSearchHit hit:targetHits) {
			System.out.println(hit.getQueryIdx()+" "+hit.getSequenceIdx()+":"+hit.getStart());
		}
		
	}

	private void processAlignment(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerHitsCluster cluster) {
		int queryLength = graph.getSequenceLength(querySequenceId);
		int targetSeqIdx = cluster.getSequenceIdx();
		int targetLength = graph.getSequenceLength(targetSeqIdx);
		double pct = 100.0*cluster.getProportionKmers();
		double queryCoverage = cluster.getQueryCoverage();
		int queryRegionLength = cluster.getQueryEnd()-cluster.getQueryStart();
		//Zero based limits
		int firstTarget = cluster.getFirst()-1;
		int lastTarget = cluster.getLast()-1;
		//System.out.println("Processing cluster. Query: "+querySequenceId+" length: "+queryLength+ " target: "+cluster.getSequenceIdx()+" length: "+targetLength);
		if(firstTarget>=0 && lastTarget<targetLength) {
			//TODO: improve rules for embedded sequences
			if(pct>=2*minKmerPercentage && queryCoverage>=0.5) {
				AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(querySequenceId, graph.getSequence(querySequenceId), queryRC, targetSeqIdx, firstTarget);
				embeddedEvent.setEvidence(cluster);
				graph.addEmbedded(embeddedEvent);
				if (querySequenceId==53) System.out.println("Query: "+querySequenceId+" embedded in "+targetSeqIdx);
			}
		} else if (firstTarget>=0) {
			//Query after target
			AssemblyVertex vertexTarget = graph.getVertex(targetSeqIdx, false);
			AssemblyVertex vertexQuery;
			if(queryRC) {
				vertexQuery = graph.getVertex(querySequenceId, false); 
			} else {
				vertexQuery = graph.getVertex(querySequenceId, true);
			}
			int overlap = targetLength-firstTarget;
			if (queryRegionLength > 0.5*overlap) {
				// TODO: design well cost function
				//double overlapCost = weightedProportion*overlap;
				//int cost = targetLength + queryLength - (int)Math.round(overlapCost);
				int cost = targetLength + queryLength - overlap;
				AssemblyEdge edge = new AssemblyEdge(vertexTarget, vertexQuery, cost, overlap);
				edge.setEvidence(cluster);
				graph.addEdge(edge);
			}
			
			//System.out.println("Edge between target: "+targetSeqIdx+" and query "+querySequenceId+" overlap: "+overlap+" weight: "+weight);
		} else if (lastTarget<targetLength) {
			//Query before target
			AssemblyVertex vertexTarget = graph.getVertex(targetSeqIdx, true);
			AssemblyVertex vertexQuery;
			if(queryRC) {
				vertexQuery = graph.getVertex(querySequenceId, true); 
			} else {
				vertexQuery = graph.getVertex(querySequenceId, false);
			}
			int overlap = lastTarget+1;
			if (cluster.getQueryCoverage()*queryLength > 0.5*overlap) {
				int cost = targetLength + queryLength -overlap;
				AssemblyEdge edge = new AssemblyEdge(vertexQuery, vertexTarget, cost, overlap);
				edge.setEvidence(cluster);
				graph.addEdge(edge);
			}
			//System.out.println("Edge between query: "+querySequenceId+" and target "+targetSeqIdx+" overlap: "+overlap+" weight: "+weight);
		} else if(pct>=2*minKmerPercentage && cluster.getQueryCoverage()>=0.5) {
			log.warning("Possible reverse embedded. Query id: "+querySequenceId+" length: "+queryLength+" target id: "+targetSeqIdx+" length: "+targetLength+" first target: "+firstTarget+" last target: "+lastTarget);
		}
	}
}