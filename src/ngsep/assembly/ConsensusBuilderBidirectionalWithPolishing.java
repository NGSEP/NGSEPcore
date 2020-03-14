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
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.LongReadsAligner;
import ngsep.alignments.ReadAlignment;
import ngsep.discovery.AlignmentsPileupGenerator;
import ngsep.discovery.IndelRealignerPileupListener;
import ngsep.discovery.VariantPileupListener;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.CalledGenomicVariant;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class ConsensusBuilderBidirectionalWithPolishing implements ConsensusBuilder {
	
	private Logger log = Logger.getLogger(ConsensusBuilderBidirectionalWithPolishing.class.getName());
	private static final String MOCK_REFERENCE_NAME = "Consensus";
	
	private LongReadsAligner aligner = new LongReadsAligner();
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	@Override
	public List<CharSequence> makeConsensus(AssemblyGraph graph) 
	{
		//List of final contigs
		List<CharSequence> consensusList = new ArrayList<CharSequence>();
		List<List<AssemblyEdge>> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			List<AssemblyEdge> path = paths.get(i);
			CharSequence consensusPath = makeConsensus (graph, path);
			consensusList.add(consensusPath);
		}
		
		return consensusList;
	}
	
	private CharSequence makeConsensus(AssemblyGraph graph, List<AssemblyEdge> path) {
		StringBuilder rawConsensus = new StringBuilder();
		AssemblyVertex lastVertex = null;
		List<ReadAlignment> alignments = new ArrayList<ReadAlignment>();
		int totalReads = 0;
		int unalignedReads = 0;
		String pathS = "";
		if(path.size()==1) {
			rawConsensus.append(path.get(0).getVertex1().getRead());
			return rawConsensus;
		}
		for(int j = 0; j < path.size(); j++) {
			//Needed to find which is the origin vertex
			AssemblyEdge edge = path.get(j);
			AssemblyVertex vertexPreviousEdge;
			AssemblyVertex vertexNextEdge;
			//If the first edge is being checked, compare to the second edge to find the origin vertex
			if(j == 0) {
				AssemblyEdge nextEdge = path.get(j + 1);
				vertexNextEdge = edge.getSharedVertex(nextEdge);
				if(vertexNextEdge== null) throw new RuntimeException("Inconsistency found in first edge of path");
				vertexPreviousEdge = edge.getVertex1();
				if(vertexPreviousEdge == vertexNextEdge) vertexPreviousEdge = edge.getVertex2();
			}
			else if (lastVertex == edge.getVertex1()) {
				vertexPreviousEdge = edge.getVertex1();
				vertexNextEdge = edge.getVertex2();
			}
			else if (lastVertex == edge.getVertex2()) {
				vertexPreviousEdge = edge.getVertex2();
				vertexNextEdge = edge.getVertex1();
			}
			else {
				throw new RuntimeException("Inconsistency found in path");
			}
			if(j == 0) {
				pathS = pathS.concat(vertexPreviousEdge.getSequenceIndex() + ",");
				String seq = vertexPreviousEdge.getRead().toString();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) seq = DNASequence.getReverseComplement(seq);
				rawConsensus.append(seq.toUpperCase());
			} 
			else if(vertexPreviousEdge.getRead()!=vertexNextEdge.getRead()) {
				//If the second string is not start, then the reverse complement is added to the consensus
				String seq = vertexNextEdge.getRead().toString();
				boolean reverse = !vertexNextEdge.isStart();
				if(reverse) seq = DNASequence.getReverseComplement(seq);
				int overlap = edge.getOverlap();
				if(seq.length() - overlap > 0) {
					pathS = pathS.concat(vertexNextEdge.getSequenceIndex() + ",");
					//String overlapSegment = nextSequence.substring(0, edge.getOverlap());
					String remainingSegment = seq.substring(edge.getOverlap());
					rawConsensus.append(remainingSegment.toUpperCase());
				} else {
					log.warning("Non embedded edge has overlap: "+edge.getOverlap()+ " and length: "+seq.length());
				}
				
			}
			if(vertexPreviousEdge.getRead()==vertexNextEdge.getRead()) {
				CharSequence read = vertexPreviousEdge.getRead();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) read = DNASequence.getReverseComplement(read.toString());
				int startConsensus = Math.max(0, rawConsensus.length() - read.length()-10);
				Map<CharSequence, Integer> uniqueKmersSubject = aligner.extractUniqueKmers(rawConsensus,startConsensus,rawConsensus.length());
				totalReads++;
				ReadAlignment alnRead = aligner.alignRead(rawConsensus, read, uniqueKmersSubject, MOCK_REFERENCE_NAME);
				if (alnRead!=null) alignments.add(alnRead);
				else unalignedReads++;
				
				List<AssemblyEmbedded> embeddedList = graph.getEmbedded(vertexPreviousEdge.getSequenceIndex());
				for(AssemblyEmbedded embedded:embeddedList) {
					CharSequence embeddedRead = embedded.getRead();
					boolean reverseE = (reverse!=embedded.isReverse());
					if(reverseE) embeddedRead = DNASequence.getReverseComplement(embeddedRead.toString());
					totalReads++;
					ReadAlignment alnEmbedded = aligner.alignRead(rawConsensus, embeddedRead, uniqueKmersSubject, MOCK_REFERENCE_NAME);
					if(alnEmbedded!=null) alignments.add(alnEmbedded);
					else unalignedReads++;
				}
			}
			lastVertex = vertexNextEdge;
		}
		log.info("Total reads: "+totalReads+" alignments: "+alignments.size()+" unaligned: "+unalignedReads);
		log.info("Path: "+pathS);
		String consensus = rawConsensus.toString();
		List<CalledGenomicVariant> variants = callVariants(consensus,alignments);
		log.info("Identified "+variants.size()+" total variants from read alignments");
		return applyVariants(consensus, variants);
	}


	private List<CalledGenomicVariant> callVariants(String consensus, List<ReadAlignment> alignments) {
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		ReferenceGenome genome = new ReferenceGenome(new QualifiedSequence(MOCK_REFERENCE_NAME, consensus));
		generator.setSequencesMetadata(genome.getSequencesMetadata());
		IndelRealignerPileupListener realignerListener = new IndelRealignerPileupListener();
		realignerListener.setGenome(genome);
		generator.addListener(realignerListener);
		VariantPileupListener varListener = new VariantPileupListener();
		varListener.setGenome(genome);
		generator.addListener(varListener);
		Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
		for(ReadAlignment aln:alignments) generator.processAlignment(aln);
		return varListener.getCalledVariants();
	}

	private CharSequence applyVariants(String consensus, List<CalledGenomicVariant> variants) {
		StringBuilder polishedConsensus = new StringBuilder();
		int l = consensus.length();
		int nextPos = 1;
		int appliedVariants = 0;
		for(CalledGenomicVariant call:variants) {
			if(call.isUndecided()|| call.isHeterozygous()) continue;
			appliedVariants++;
			String [] calledAlleles = call.getCalledAlleles();
			if(nextPos<call.getFirst()) {
				//Fill haplotypes with non variant segment
				String segment = consensus.substring(nextPos-1, call.getFirst()-1);
				polishedConsensus.append(segment);
			}
			polishedConsensus.append(calledAlleles[0]);
			nextPos = call.getLast()+1;
		}
		if(nextPos<=l) {
			//End of a chromosome
			CharSequence nonVarLast = consensus.substring(nextPos-1);
			polishedConsensus.append(nonVarLast);
		}
		log.info("Applied "+appliedVariants+" homozygous alternative variants");
		return new DNASequence(polishedConsensus.toString());
	}
}
