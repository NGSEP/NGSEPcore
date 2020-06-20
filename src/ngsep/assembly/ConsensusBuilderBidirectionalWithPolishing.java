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

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.discovery.AlignmentsPileupGenerator;
import ngsep.discovery.IndelRealignerPileupListener;
import ngsep.discovery.SingleSampleVariantPileupListener;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmerHitsCluster;
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
	
	private MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm();
	
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
				pathS = pathS.concat(vertexPreviousEdge.getUniqueNumber() + ",");
				CharSequence seq = vertexPreviousEdge.getRead().getCharacters();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
				rawConsensus.append(seq);
			} 
			else if(vertexPreviousEdge.getRead()!=vertexNextEdge.getRead()) {
				// Augment consensus with the next path read
				CharSequence nextPathSequence = vertexNextEdge.getRead().getCharacters();
				boolean reverse = !vertexNextEdge.isStart();
				if(reverse) nextPathSequence = DNAMaskedSequence.getReverseComplement(nextPathSequence);
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) printAllOverlappingSeqs(graph,path,j,vertexPreviousEdge);
				
				int overlap = edge.getOverlap();
				int startSuffix = overlap;
				KmerHitsCluster cluster = edge.getEvidence();
				if(cluster!=null) {
					if (vertexPreviousEdge.getSequenceIndex()== cluster.getSequenceIdx()) {
						if (!vertexPreviousEdge.isStart()) startSuffix = cluster.getQueryPredictedEnd();
						else startSuffix = nextPathSequence.length() - cluster.getQueryPredictedStart();
					} else if (vertexNextEdge.getSequenceIndex() == cluster.getSequenceIdx()) {
						if (vertexNextEdge.isStart()) startSuffix = cluster.getSubjectPredictedEnd();
						else startSuffix = nextPathSequence.length() - cluster.getSubjectPredictedStart();
					}
				}
				
				if(startSuffix<nextPathSequence.length()) {
					pathS = pathS.concat(vertexNextEdge.getUniqueNumber() + ",");
					String remainingSegment = nextPathSequence.subSequence(startSuffix, nextPathSequence.length()).toString();
					//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) System.out.println("Consensus length: "+rawConsensus.length()+" Vertex: "+vertexNextEdge.getUniqueNumber()+" read length: "+nextPathSequence.length()+" overlap: "+edge.getOverlap()+" remaining: "+remainingSegment.length());
					rawConsensus.append(remainingSegment.toUpperCase());
				} else {
					log.warning("Non embedded edge btween vertices"+vertexPreviousEdge.getUniqueNumber()+" and "+vertexNextEdge.getUniqueNumber()+" has overlap: "+overlap+ " start suffix "+startSuffix+" and length: "+nextPathSequence.length());
				}
				
			}
			if(vertexPreviousEdge.getRead()==vertexNextEdge.getRead()) {
				//Align to consensus next path read and its embedded sequences
				CharSequence read = vertexPreviousEdge.getRead().getCharacters();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) read = DNAMaskedSequence.getReverseComplement(read);
				Map<CharSequence, Integer> uniqueKmersSubject = aligner.extractUniqueKmers(rawConsensus,Math.max(0, rawConsensus.length()-read.length()),rawConsensus.length());
				totalReads++;
				ReadAlignment alnRead = aligner.alignRead(vertexPreviousEdge.getSequenceIndex(), rawConsensus, read, uniqueKmersSubject, 0.5);
				if (alnRead!=null) {
					alnRead.setSequenceName(MOCK_REFERENCE_NAME);
					//alnRead.setQualityScores(RawRead.generateFixedQSString('5', read.length()));
					alignments.add(alnRead);
				}
				else unalignedReads++;
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) System.out.println("Consensus length: "+rawConsensus.length()+" Vertex: "+vertexNextEdge.getUniqueNumber()+" sequence: "+read.length()+" alignment: "+alnRead);
				if (totalReads%100==0) log.info("Aligning. Processed reads: "+totalReads+" alignments: "+alignments.size()+" unaligned: "+unalignedReads);
				
				List<AssemblyEmbedded> embeddedList = graph.getAllEmbedded(vertexPreviousEdge.getSequenceIndex());
				//List<AssemblyEmbedded> embeddedList = graph.getEmbeddedByHostId(vertexPreviousEdge.getSequenceIndex());
				for(AssemblyEmbedded embedded:embeddedList) {
					CharSequence embeddedRead = embedded.getRead();
					boolean reverseE = (reverse!=embedded.isReverse());
					if(reverseE) embeddedRead = DNAMaskedSequence.getReverseComplement(embeddedRead);
					totalReads++;
					ReadAlignment alnEmbedded = aligner.alignRead(embedded.getSequenceId(), rawConsensus, embeddedRead, uniqueKmersSubject, 0.5);
					if(alnEmbedded!=null) {
						alnEmbedded.setSequenceName(MOCK_REFERENCE_NAME);
						//alnEmbedded.setQualityScores(RawRead.generateFixedQSString('5', read.length()));
						alignments.add(alnEmbedded);
					}
					else unalignedReads++;
					if (totalReads%100==0) log.info("Aligning. Processed reads: "+totalReads+" alignments: "+alignments.size()+" unaligned: "+unalignedReads);
				}
			}
			lastVertex = vertexNextEdge;
		}
		log.info("Total reads: "+totalReads+" alignments: "+alignments.size()+" unaligned: "+unalignedReads);
		log.info("Path: "+pathS);
		String consensus = rawConsensus.toString();
		//return consensus;
		List<CalledGenomicVariant> variants = callVariants(consensus,alignments);
		log.info("Identified "+variants.size()+" total variants from read alignments");
		return applyVariants(consensus, variants);
	}


	public void printAllOverlappingSeqs(AssemblyGraph graph, List<AssemblyEdge> path, int pathPos, AssemblyVertex vertexPreviousEdge) {
		for(int j = pathPos; j < path.size(); j++) {
			AssemblyEdge edge = path.get(j);
			if(edge.isSameSequenceEdge()) continue;
			AssemblyEdge alt1 = graph.getEdge(vertexPreviousEdge, edge.getVertex1());
			AssemblyEdge alt2 = graph.getEdge(vertexPreviousEdge, edge.getVertex2());
			if(alt1==null && alt2==null) continue;
			if(alt1!=null) System.out.println("Edge between: "+alt1.getVertex1().getUniqueNumber()+" and "+alt1.getVertex2().getUniqueNumber()+" overlap: "+alt1.getOverlap());
			if(alt2!=null) System.out.println("Edge between: "+alt2.getVertex1().getUniqueNumber()+" and "+alt2.getVertex2().getUniqueNumber()+" overlap: "+alt2.getOverlap());	
		}
		
	}

	private List<CalledGenomicVariant> callVariants(String consensus, List<ReadAlignment> alignments) {
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		generator.setLog(log);
		ReferenceGenome genome = new ReferenceGenome(new QualifiedSequence(MOCK_REFERENCE_NAME, consensus));
		generator.setSequencesMetadata(genome.getSequencesMetadata());
		generator.setMaxAlnsPerStartPos(100);
		generator.setMinMQ(80);
		IndelRealignerPileupListener realignerListener = new IndelRealignerPileupListener();
		realignerListener.setGenome(genome);
		generator.addListener(realignerListener);
		SingleSampleVariantPileupListener varListener = new SingleSampleVariantPileupListener();
		varListener.setMinQuality((short) 30);
		varListener.setGenome(genome);
		generator.addListener(varListener);
		Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
		int count = 0;
		for(ReadAlignment aln:alignments) {
			generator.processAlignment(aln);
			count++;
			if(count%100==0) log.info("Processed "+count+" alignments. Called variants: "+varListener.getCalledVariants().size()); 
		}
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
		return new DNAMaskedSequence(polishedConsensus.toString());
	}
}
