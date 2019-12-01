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
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.KmerAlignmentCluster;
import ngsep.alignments.KmerWithStart;
import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;

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


	private boolean updateGraph(AssemblyGraph graph, int querySequenceId, String sequence, boolean queryRC, FMIndex fmIndex) {
		boolean isEmbedded = false;
		List<KmerWithStart> kmers = KmerWithStart.selectKmers(sequence, kmerLength, kmerOffset);
		if(kmers==null) return isEmbedded;
		
		int kmersCount=kmers.size();
		List<ReadAlignment> initialKmerAlns = searchKmers (querySequenceId, kmers, fmIndex);
		
		List<KmerAlignmentCluster> clusteredKmerAlns = clusterKmerAlignments(sequence, initialKmerAlns, kmersCount);
		//System.out.println("Query id: "+querySequenceId+" length "+sequence.length()+" RC: "+queryRC+" kmers: "+kmers.size()+" Kmer alignments: "+initialKmerAlns.size()+" Clusters: "+clusteredKmerAlns.size());
		Collections.sort(clusteredKmerAlns, new Comparator<KmerAlignmentCluster>() {
			@Override
			public int compare(KmerAlignmentCluster o1, KmerAlignmentCluster o2) {
				return o2.getNumDifferentKmers()-o1.getNumDifferentKmers();
			}
		});

		int kmersMaxCluster = 0;
		for (int i=0;i<clusteredKmerAlns.size() && i<10;i++) {
			KmerAlignmentCluster cluster = clusteredKmerAlns.get(i);
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
	private List<ReadAlignment> searchKmers(int querySequenceId, List<KmerWithStart> kmers, FMIndex fmIndex) {
		List<ReadAlignment> answer = new ArrayList<>();
		for (KmerWithStart kmer:kmers) {
			List<ReadAlignment> kmerAlns=fmIndex.search(kmer.getKmer().toString(),0,querySequenceId-1);
			//if(querySequenceId==1) System.out.println("Query: "+querySequenceId+" Found "+kmerAlns.size()+" alignments for kmer: "+kmer.getKmer().toString());
			for(ReadAlignment aln:kmerAlns) {
				//if(querySequenceId==1) System.out.println("Kmer start: "+kmer.getStart()+" Next alignment: "+aln.getSequenceIndex()+": "+aln.getFirst()+"-"+aln.getLast()+" rc: "+aln.isNegativeStrand());
				aln.setReadNumber(kmer.getStart());
				answer.add(aln);
			}
		}
		return answer;
	}
	private List<KmerAlignmentCluster> clusterKmerAlignments(CharSequence query, List<ReadAlignment> initialKmerAlns, int numKmers) {
		List<KmerAlignmentCluster> clusters = new ArrayList<>();
		Map<Integer,List<ReadAlignment>> alnsByTargetSequence = new HashMap<>();
		for(ReadAlignment aln: initialKmerAlns) {
			int targetSequenceId = aln.getSequenceIndex();
			List<ReadAlignment> targetHits = alnsByTargetSequence.get(targetSequenceId);
			if(targetHits==null) {
				targetHits = new ArrayList<>();
				alnsByTargetSequence.put(targetSequenceId, targetHits);
			}
			targetHits.add(aln);
		}
		
		for(int targetIdx:alnsByTargetSequence.keySet()) {
			List<ReadAlignment> targetHits = alnsByTargetSequence.get(targetIdx);
			
			double pct = 100.0*(double)targetHits.size()/numKmers;
			if(pct<minKmerPercentage) continue;
			Collections.sort(targetHits,GenomicRegionPositionComparator.getInstance());
			//if(targetIdx==0) printTargetHits(targetHits);
			clusters.addAll(clusterSequenceKmerAlns(query, targetHits));
		}
		return clusters;
	}
	private Collection<KmerAlignmentCluster> clusterSequenceKmerAlns(CharSequence query, List<ReadAlignment> sequenceAlns) {
		Collection<KmerAlignmentCluster> answer = new ArrayList<>();
		//System.out.println("Alns to cluster: "+sequenceAlns.size());
		for(ReadAlignment aln:sequenceAlns) {
			boolean clustered = false;
			for(KmerAlignmentCluster cluster:answer) {
				if(cluster.addAlignment(aln,50)) {
					clustered=true;
					break;
				}
			}
			if(!clustered) {
				answer.add(new KmerAlignmentCluster(query, aln));
			}
		}
		return answer;
	}
	public void printTargetHits(List<ReadAlignment> targetHits) {
		for(ReadAlignment aln:targetHits) {
			System.out.println(aln.getReadNumber()+" "+aln.getSequenceIndex()+":"+aln.getFirst()+"-"+aln.getLast());
		}
		
	}

	private boolean processAlignment(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerAlignmentCluster cluster) {
		boolean isEmbedded = false;
		int queryLength = graph.getSequenceLength(querySequenceId);
		int targetSeqIdx = cluster.getSequenceIdx();
		int targetLength = graph.getSequenceLength(targetSeqIdx);
		int firstTarget = cluster.getFirst();
		int lastTarget = cluster.getLast();
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