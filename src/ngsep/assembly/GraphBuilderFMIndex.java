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
		log.info("Created graph vertices");
		Collections.sort(sequences, (l1, l2) -> l2.length() - l1.length());
		
		log.info("Sorted sequences");
		// Create FM-Index
		FMIndex fmIndex = new FMIndex();
		fmIndex.loadUnnamedSequences(sequences, TALLY_DISTANCE, SUFFIX_FRACTION);
		log.info("Created FM-Index");
	
		for (int seqId = 0; seqId < sequences.size(); seqId++) {
			if (graph.isEmbedded(seqId)) {
				continue;
			}
			String seq = sequences.get(seqId).toString();
			updateGraph(graph, seqId, seq, false, fmIndex);
			String complement = DNAMaskedSequence.getReverseComplement(seq).toString();
			updateGraph(graph, seqId, complement, true, fmIndex);
			if ((seqId+1)%100==0) log.info("Processed "+(seqId+1) +" sequences");
		}
		log.info("Built graph");
		//TODO: Check if needed
		//assemblyGraph.removeAllEmbeddedsIntoGraph();
		//assemblyGraph.ExtrapolateAligns();
		return graph;
	}


	private void updateGraph(AssemblyGraph graph, int querySequenceId, String sequence, boolean queryRC, FMIndex fmIndex) {
		List<KmerWithStart> kmers = KmerWithStart.selectKmers(sequence, kmerLength, kmerOffset);
		if(kmers==null) return;
		//System.out.println("Query: "+query.toString()+" kmers: "+kmers.size());
		int kmersCount=kmers.size();
		List<ReadAlignment> initialKmerAlns = searchKmers (querySequenceId, kmers, fmIndex);
		List<KmerAlignmentCluster> clusteredKmerAlns = clusterKmerAlignments(sequence, initialKmerAlns, kmersCount); 
		//System.out.println("Clusters: "+clusteredKmerAlns.size());
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
			processAlignment(graph, querySequenceId, queryRC, cluster);
			
		}
	}

	/**
	 * Searches the given kmers in the fmIndex 
	 * @param kmers to search
	 * @return List of alignments of each kmer. The read number of each alignment contains the kmer number.
	 */
	private List<ReadAlignment> searchKmers(int querySequenceId, List<KmerWithStart> kmers, FMIndex fmIndex) {
		List<ReadAlignment> answer = new ArrayList<>();
		for (KmerWithStart kmer:kmers) {
			List<ReadAlignment> kmerAlns=fmIndex.search(kmer.getKmer().toString());
			for(ReadAlignment aln:kmerAlns) {
				if(aln.getSequenceIndex()>querySequenceId) {
					aln.setReadNumber(kmer.getStart());
					answer.add(aln);
				}
				
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
		
		for(List<ReadAlignment> targetHits:alnsByTargetSequence.values()) {
			double pct = 100.0*(double)targetHits.size()/numKmers;
			if(pct<minKmerPercentage) continue;
			Collections.sort(targetHits,GenomicRegionPositionComparator.getInstance());
			clusters.addAll(KmerAlignmentCluster.clusterSequenceKmerAlns(query, targetHits));
		}
		return clusters;
	}
	private void processAlignment(AssemblyGraph graph, int querySequenceId, boolean queryRC, KmerAlignmentCluster cluster) {
		int targetSeqIdx = cluster.getSequenceIdx();
		int sequenceLength = graph.getSequenceLength(targetSeqIdx);
		int firstTarget = cluster.getFirst();
		int lastTarget = cluster.getLast();
		if(firstTarget>=0 && lastTarget<sequenceLength) {
			//Embedded
			AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(querySequenceId, graph.getSequence(querySequenceId), firstTarget-1, queryRC);
			graph.addEmbedded(targetSeqIdx, embeddedEvent);
		} else if (firstTarget>=0) {
			//Query after the target
			AssemblyVertex vertexTarget = graph.getVertex(targetSeqIdx, false);
			AssemblyVertex vertexQuery;
			if(queryRC) {
				vertexQuery = graph.getVertex(querySequenceId, false); 
			} else {
				vertexQuery = graph.getVertex(querySequenceId, true);
			}
			graph.addEdge(vertexTarget, vertexQuery, sequenceLength-firstTarget-1);
		} else {
			//Query before target
			AssemblyVertex vertexTarget = graph.getVertex(targetSeqIdx, true);
			AssemblyVertex vertexQuery;
			if(queryRC) {
				vertexQuery = graph.getVertex(querySequenceId, true); 
			} else {
				vertexQuery = graph.getVertex(querySequenceId, false);
			}
			graph.addEdge(vertexQuery, vertexTarget, lastTarget-1);
		}
	}
}