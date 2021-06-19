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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;

public class AlignmentBasedIndelErrorsCorrector {
	
	private Logger log = Logger.getAnonymousLogger();
	private LayoutBuilderKruskalPath pathsFinder;
	private AssemblySequencesRelationshipFilter filter;
	private int numThreads=1;
	
	
	
	public AlignmentBasedIndelErrorsCorrector() {
		super();
		pathsFinder = new LayoutBuilderKruskalPath();
		pathsFinder.setMinPathLength(1);
		pathsFinder.setRunImprovementAlgorithms(false);
		filter = new AssemblySequencesRelationshipFilter();	
	}

	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
		pathsFinder.setLog(log);
	}

	public int getNumThreads() {
		return numThreads;
	}

	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}

	public void correctErrors(AssemblyGraph graph) {
		AssemblyPathReadsAligner aligner = new AssemblyPathReadsAligner();
		aligner.setLog(log);
		aligner.setNumThreads(numThreads);
		aligner.setAlignEmbedded(true);
		AssemblyGraph copyGraph = graph.buildSubgraph(null);
		copyGraph.removeVerticesChimericReads();
		copyGraph.updateScores(0.5);
		filter.filterEdgesAndEmbedded(copyGraph, 0.5);
		pathsFinder.findPaths(copyGraph);
		int n =graph.getNumSequences();
		int [] sequencePaths = new int [n];
		Arrays.fill(sequencePaths, 0);
		List<AssemblyPath> paths = copyGraph.getPaths();
		for(int i = 0; i < paths.size(); i++)
		{
			AssemblyPath path = paths.get(i);
			path.setPathId(i+1);
			correctErrors(graph, copyGraph, path, sequencePaths, aligner);
		}
		
	}

	private void correctErrors(AssemblyGraph graph, AssemblyGraph copyGraph, AssemblyPath path, int[] sequencePaths, AssemblyPathReadsAligner aligner) {
		String sequenceName = "Consensus_"+path.getPathId();
		aligner.alignPathReads(copyGraph, path);
		StringBuilder rawConsensus = aligner.getConsensus();
		List<ReadAlignment> alignments = aligner.getAlignedReads();
		for(ReadAlignment aln:alignments) aln.setSequenceName(sequenceName);
		Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
		//TODO: Define better ploidy
		List<CalledGenomicVariant> calledVars = ConsensusBuilderBidirectionalWithPolishing.callVariants(sequenceName, rawConsensus, alignments, 2);
		List<CalledGenomicVariant> filteredVars = new ArrayList<CalledGenomicVariant>();
		for(CalledGenomicVariant call:calledVars) {
			if(!call.isUndecided() && !call.isHomozygousReference()) filteredVars.add(call);
		}
		int indexNextActive = 0;
		int correctedErrors = 0;
		int correctedReads = 0;
		
		for(ReadAlignment aln:alignments) {
			int readId = aln.getReadNumber();
			if(sequencePaths[readId]>0) {
				log.warning("Read "+readId+" "+aln.getReadName()+" already correctd by path: "+sequencePaths[readId]+" current path: "+path.getPathId());
				continue;
			}
			while (indexNextActive<filteredVars.size()) {
				GenomicRegion activeSegment = filteredVars.get(indexNextActive);
				if(activeSegment.getLast()<aln.getFirst()) {
					indexNextActive++;
				} else break;
			}
			String alignedRead = aln.getReadCharacters().toString();
			StringBuilder correctedRead = new StringBuilder(alignedRead.length());
			Map<Integer,GenomicVariant> calls = aln.getIndelCallsByAlignedReadPos();
			if(calls == null) continue;
			int nextPos = 0;
			int lastRef = 0;
			for(Map.Entry<Integer, GenomicVariant> entry:calls.entrySet()) {
				int posRead = entry.getKey();
				GenomicVariant indel = entry.getValue();
				if(indel.length()>2) continue;
				if(posRead>nextPos) {
					//Apply homozygous indels not called in this alignment
					for(int j=indexNextActive;j<filteredVars.size();j++) {
						CalledGenomicVariant calledVar = filteredVars.get(j);
						if(calledVar.isHeterozygous()) continue;
						if(calledVar.getFirst() > lastRef && calledVar.getLast()<indel.getFirst()) {
							int readPosStartVar = aln.getAlignedReadPosition(calledVar.getFirst());
							int readPosEndVar = aln.getAlignedReadPosition(calledVar.getLast());
							if(readPosStartVar> nextPos && readPosEndVar<posRead && readPosStartVar<readPosEndVar) {
								if(readPosEndVar-readPosStartVar>10) System.err.println("WARN: Correcting indel spanning "+readPosStartVar+" "+readPosEndVar+" segment: "+alignedRead.substring(readPosStartVar,readPosEndVar+1)+" alleles var: "+calledVar.getAlleles()[0]+" "+calledVar.getAlleles()[1]);
								correctedRead.append(alignedRead.substring(nextPos, readPosStartVar));
								correctedRead.append(calledVar.getCalledAlleles()[0]);
								nextPos = readPosEndVar+1;
							}
						} else {
							break;
						}
					}
					correctedRead.append(alignedRead.substring(nextPos, posRead+1));
				}
				nextPos = posRead+1;
				//Correct if indel not called in this region
				boolean correctIndel = true;
				for(int j=indexNextActive;j<filteredVars.size();j++) {
					GenomicRegion activeSegment = filteredVars.get(j);
					if(activeSegment.getLast()>=indel.getFirst()) {
						correctIndel = (activeSegment.getFirst()>indel.getLast());
						break;
					}
				}
				if(correctIndel) {
					if(indel.getLast()==indel.getFirst()+1) {
						//Insertion
						nextPos +=indel.length();
					} else {
						//Deletion
						correctedRead.append(rawConsensus.substring(indel.getFirst(), indel.getLast()-1));
					}
					correctedErrors++;
				}
				lastRef = indel.getLast();
			}
			if(nextPos<rawConsensus.length()) correctedRead.append(alignedRead.substring(nextPos));
			aln.setReadCharacters(null);
			CharSequence correctedReadS;
			if(aln.isNegativeStrand()) {
				correctedReadS = DNAMaskedSequence.getReverseComplement(correctedRead);
			} else {
				correctedReadS = new DNAMaskedSequence(correctedRead);
			}
			QualifiedSequence qseq = graph.getSequence(readId);
			qseq.setCharacters(correctedReadS);
			sequencePaths[readId] = path.getPathId();
			correctedReads++;
		}
		log.info("IndelErrorsCorrector. Path: "+path.getPathId()+" Corrected reads: "+correctedReads+" correctedErrors: "+correctedErrors);
	}
}
