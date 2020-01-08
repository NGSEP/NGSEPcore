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
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.FMIndexUngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersCounter;

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
			boolean isEmbedded = updateGraph(graph, seqId, seq, false, fmIndex);
			if(!isEmbedded) {
				String complement = DNAMaskedSequence.getReverseComplement(seq).toString();
				updateGraph(graph, seqId, complement, true, fmIndex);
			}
			if ((seqId+1)%100==0) log.info("Processed "+(seqId+1) +" sequences. Number of edges: "+graph.getEdges().size());
		}
		log.info("Built graph. Edges: "+graph.getEdges().size()+" Prunning embedded sequences");
		graph.pruneEmbeddedSequences();
		log.info("Prunned graph. Edges: "+graph.getEdges().size());
		
		return graph;
	}


	private boolean updateGraph(AssemblyGraph graph, int querySequenceId, String query, boolean queryRC, FMIndex fmIndex) {
		boolean isEmbedded = false;
		Map<Integer,CharSequence> kmersMap = KmersCounter.extractKmersAsMap(query, kmerLength, kmerOffset, true, true, true);
		
		int kmersCount=kmersMap.size();
		if(kmersCount==0) return isEmbedded;
		List<FMIndexUngappedSearchHit> initialKmerHits = searchKmers (querySequenceId, kmersMap, fmIndex);
		
		List<KmerHitsCluster> clusteredKmerAlns = clusterKmerHits(query, initialKmerHits, kmersCount);
		//System.out.println("Query id: "+querySequenceId+" length "+sequence.length()+" RC: "+queryRC+" kmers: "+kmers.size()+" Kmer alignments: "+initialKmerAlns.size()+" Clusters: "+clusteredKmerAlns.size());
		Collections.sort(clusteredKmerAlns, (o1,o2)-> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());

		int kmersMaxCluster = 0;
		for (int i=0;i<clusteredKmerAlns.size() && i<10;i++) {
			KmerHitsCluster cluster = clusteredKmerAlns.get(i);
			int numKmersCluster = cluster.getNumDifferentKmers();
			double pct = 100.0*(double)numKmersCluster/kmersCount;
			if(pct<minKmerPercentage) break;
			if(i==0) kmersMaxCluster = cluster.getNumDifferentKmers();
			else if (2*numKmersCluster<kmersMaxCluster) break;
			isEmbedded = processAlignment(graph, querySequenceId, queryRC, cluster);
			if(isEmbedded) break;
		}
		return isEmbedded;
	}

	/**
	 * Searches the given kmers in the fmIndex 
	 * @param kmers to search
	 * @return List of alignments of each kmer. The read number of each alignment contains the kmer number.
	 */
	private List<FMIndexUngappedSearchHit> searchKmers(int querySequenceId, Map<Integer,CharSequence> kmersMap, FMIndex fmIndex) {
		List<FMIndexUngappedSearchHit> answer = new ArrayList<>();
		for (int start:kmersMap.keySet()) {
			String kmer = kmersMap.get(start).toString();
			List<FMIndexUngappedSearchHit> kmerHits=fmIndex.exactSearch(kmer,0,querySequenceId-1);
			//if(querySequenceId==1) System.out.println("Query: "+querySequenceId+" Found "+kmerAlns.size()+" alignments for kmer: "+kmer.getKmer().toString());
			for(FMIndexUngappedSearchHit hit:kmerHits) {
				//if(querySequenceId==1) System.out.println("Kmer start: "+kmer.getStart()+" Next alignment: "+aln.getSequenceIndex()+": "+aln.getFirst()+"-"+aln.getLast()+" rc: "+aln.isNegativeStrand());
				hit.setQueryIdx(start);
				answer.add(hit);
			}
		}
		return answer;
	}
	private List<KmerHitsCluster> clusterKmerHits(CharSequence query, List<FMIndexUngappedSearchHit> initialKmerHits, int numKmers) {
		List<KmerHitsCluster> clusters = new ArrayList<>();
		Map<Integer,List<FMIndexUngappedSearchHit>> hitsByTargetSequence = new HashMap<>();
		for(FMIndexUngappedSearchHit kmerHit: initialKmerHits) {
			int targetSequenceId = kmerHit.getSequenceIdx();
			List<FMIndexUngappedSearchHit> targetHits = hitsByTargetSequence.computeIfAbsent(targetSequenceId, k-> new ArrayList<FMIndexUngappedSearchHit>());
			targetHits.add(kmerHit);
		}
		
		for(int targetIdx:hitsByTargetSequence.keySet()) {
			List<FMIndexUngappedSearchHit> targetHits = hitsByTargetSequence.get(targetIdx);
			
			double pct = 100.0*(double)targetHits.size()/numKmers;
			if(pct<minKmerPercentage) continue;
			Collections.sort(targetHits,new Comparator<FMIndexUngappedSearchHit>() {

				@Override
				public int compare(FMIndexUngappedSearchHit hit0, FMIndexUngappedSearchHit hit1) {
					return hit0.getQueryIdx()-hit1.getQueryIdx();
				}
			});
			//if(targetIdx==0) printTargetHits(targetHits);
			clusters.addAll(clusterSequenceKmerAlns(query, targetHits));
		}
		return clusters;
	}
	public static List<KmerHitsCluster> clusterSequenceKmerAlns(CharSequence query, List<FMIndexUngappedSearchHit> sequenceKmerHits) {
		List<KmerHitsCluster> answer = new ArrayList<>();
		//System.out.println("Alns to cluster: "+sequenceAlns.size());
		for(FMIndexUngappedSearchHit kmerHit: sequenceKmerHits) {
			boolean clustered = false;
			for(KmerHitsCluster cluster:answer) {
				if(cluster.addKmerHit(kmerHit, 50)) {
					clustered=true;
					break;
				}
			}
			if(!clustered) {
				answer.add(new KmerHitsCluster(query, kmerHit));
			}
		}
		return answer;
	}
	public void printTargetHits(List<FMIndexUngappedSearchHit> targetHits) {
		for(FMIndexUngappedSearchHit hit:targetHits) {
			System.out.println(hit.getQueryIdx()+" "+hit.getSequenceIdx()+":"+hit.getStart());
		}
		
	}

	private boolean processAlignment(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerHitsCluster cluster) {
		boolean isEmbedded = false;
		int queryLength = graph.getSequenceLength(querySequenceId);
		int targetSeqIdx = cluster.getSequenceIdx();
		int targetLength = graph.getSequenceLength(targetSeqIdx);
		//Zero based limits
		int firstTarget = cluster.getFirst()-1;
		int lastTarget = cluster.getLast()-1;
		//System.out.println("Processing cluster. Query: "+querySequenceId+" length: "+queryLength+ " target: "+cluster.getSequenceIdx()+" length: "+targetLength);
		if(firstTarget>=0 && lastTarget<targetLength) {
			//Embedded
			AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(querySequenceId, graph.getSequence(querySequenceId), firstTarget, queryRC);
			graph.addEmbedded(targetSeqIdx, embeddedEvent);
			isEmbedded = true;
			//System.out.println("Query: "+querySequenceId+" embedded in "+targetSeqIdx);
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
			int cost = targetLength + queryLength - overlap;
			graph.addEdge(vertexTarget, vertexQuery, cost, overlap);
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
			int cost = targetLength + queryLength -overlap;
			graph.addEdge(vertexQuery, vertexTarget, cost, overlap);
			//System.out.println("Edge between query: "+querySequenceId+" and target "+targetSeqIdx+" overlap: "+overlap+" weight: "+weight);
		} else {
			log.warning("Possible reverse embedded. Query id: "+querySequenceId+" length: "+queryLength+" target id: "+targetSeqIdx+" length: "+targetLength+" first target: "+firstTarget+" last target: "+lastTarget);
		}
		return isEmbedded;
	}
}