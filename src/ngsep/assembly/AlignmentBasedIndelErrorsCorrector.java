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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
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
		filter = new AssemblySequencesRelationshipFilter();	
	}

	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
		pathsFinder.setLog(log);
	}

	/**
	 * @return Number of threads
	 */
	public int getNumThreads() {
		return numThreads;
	}

	/**
	 * Changes the number of threads
	 * @param numThreads
	 */
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	/**
	 * Corrects errors in reads
	 * @param graph Input graph with reads
	 */
	public void correctErrors(AssemblyGraph graph) {
		Runtime runtime = Runtime.getRuntime();
		long start = System.currentTimeMillis(); 
		AssemblyPathReadsAligner aligner = new AssemblyPathReadsAligner();
		aligner.setLog(log);
		aligner.setNumThreads(numThreads);
		aligner.setAlignEmbedded(true);
		AssemblyGraph copyGraph = graph.buildSubgraph(null);
		copyGraph.removeVerticesChimericReads();
		copyGraph.updateScores(0);
		filter.filterEdgesAndEmbedded(copyGraph, 0.5);
		pathsFinder.findPaths(copyGraph);
		int n =graph.getNumSequences();
		List<AssemblyPath> paths = copyGraph.getPaths();
		Map<Integer,AssemblyPath> selectedPathsMap = new HashMap<Integer, AssemblyPath>(paths.size());
		Map<Integer,QualifiedSequence> selectedPathsConsensus = new HashMap<Integer, QualifiedSequence>(paths.size());
		Map<Integer,List<ReadAlignment>> selectedPathsAlns = new HashMap<Integer, List<ReadAlignment>>(paths.size());
		Set<Integer> alignedReadIds = new HashSet<Integer>();
		
		for(int i = 0; i < paths.size(); i++)
		{
			AssemblyPath path = paths.get(i);
			int pathId = i+1;
			path.setPathId(pathId);
			String sequenceName = ""+pathId;
			if(path.getPathLength()>10) {
				aligner.alignPathReads(copyGraph, path);
				String pathConsensus = aligner.getConsensus().toString();
				List<ReadAlignment> alignments = aligner.getAlignedReads();
				for(ReadAlignment aln:alignments) {
					aln.setSequenceName(sequenceName);
					alignedReadIds.add(aln.getReadNumber());
				}
				selectedPathsMap.put(pathId, path);
				selectedPathsConsensus.put(pathId, new QualifiedSequence(sequenceName,pathConsensus));
				selectedPathsAlns.put(pathId, alignments);
			}
		}
		long endAln = System.currentTimeMillis();
		long timeAln = (endAln-start)/1000; 
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		log.info("Aligned reads within paths: "+alignedReadIds.size()+". Time Alignment (s): "+timeAln+" Memory: "+usedMemory+" Aligning remaining reads.");
		QualifiedSequenceList sequences = new QualifiedSequenceList(selectedPathsConsensus.values());
		MinimizersTableReadAlignmentAlgorithm aligner2 = new MinimizersTableReadAlignmentAlgorithm(MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		ReferenceGenome genome = new ReferenceGenome(sequences);
		try {
			aligner2.loadGenome(genome, 15, 30, numThreads, false);
		} catch (InterruptedException e) {
			log.severe("Error building index for assembled genome "+e.getMessage());
			e.printStackTrace();
			return;
		}
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		log.info("Built minimizers for consensus sequences of: "+genome.getNumSequences()+" paths. Memory: "+usedMemory);
		int numRemaininigReads = 0;
		ThreadPoolManager poolAlign = new ThreadPoolManager(numThreads, 1000);
		for(int i=0;i<n;i++) {
			if(alignedReadIds.contains(i)) continue;
			numRemaininigReads++;
			QualifiedSequence seq = graph.getSequence(i);
			if(seq==null) continue;
			final int id = i;
			try {
				poolAlign.queueTask(()->alignReadProcess(aligner2, id, seq, selectedPathsAlns));
			} catch (InterruptedException e) {
				// TODO: Better handling
				e.printStackTrace();
			}
		}
		try {
			poolAlign.terminatePool();
		} catch (InterruptedException e) {
			// TODO Better handling
			e.printStackTrace();
		}
		
		int numAligned = 0;
		int [] sequencePaths = new int [n];
		Arrays.fill(sequencePaths, 0);
		for(Map.Entry<Integer,List<ReadAlignment>> entry:selectedPathsAlns.entrySet()) {
			int pathId = entry.getKey();
			List<ReadAlignment> alns = entry.getValue();
			numAligned+=alns.size();
			for(ReadAlignment aln:alns) sequencePaths[aln.getReadNumber()] = pathId;
		}
		long endAln2 = System.currentTimeMillis();
		long timeAln2 = (endAln2-endAln)/1000;
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		log.info("Aligned "+numRemaininigReads+" remaining reads. Final number unaligned: "+(graph.getNumSequences()-numAligned)+". Time process (s): "+timeAln2+" Memory: "+usedMemory);
		ThreadPoolManager poolCorrection = new ThreadPoolManager(numThreads, 2*numThreads);
		long startCorrection = System.currentTimeMillis();
		for(Map.Entry<Integer, AssemblyPath>entry:selectedPathsMap.entrySet())
		{
			int pathId = entry.getKey();
			AssemblyPath path = entry.getValue();
			log.info("Correcting errors for reads aligned to path: "+pathId);
			QualifiedSequence consensus = selectedPathsConsensus.get(pathId);
			try {
				poolCorrection.queueTask(()->correctErrors(graph, selectedPathsAlns.get(path.getPathId()), path, consensus.getName(), consensus.getCharacters().toString(), sequencePaths));
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		try {
			poolCorrection.terminatePool();
		} catch (InterruptedException e) {
			// TODO Better handling
			e.printStackTrace();
		}
		long endCorr = System.currentTimeMillis();
		long timeCorr = (endCorr-startCorrection)/1000; 
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		log.info("Corrected errors within paths: "+alignedReadIds.size()+". Time error correction (s): "+timeCorr+" Memory: "+usedMemory);
		
	}

	private void alignReadProcess(MinimizersTableReadAlignmentAlgorithm aligner, int id, QualifiedSequence seq, Map<Integer,List<ReadAlignment>> selectedPathsAlns) {
		List<ReadAlignment> alns = aligner.alignRead(new RawRead(seq.getName(), seq.getCharacters(), null));
		if(alns.size()==0) {
			System.out.println("Unaligned read "+seq.getName()+" to consensus");
			return;
		}
		ReadAlignment aln = alns.get(0);
		aln.setReadNumber(id);
		int pathId = Integer.parseInt(aln.getSequenceName());
		List<ReadAlignment> pathAlns = selectedPathsAlns.get(pathId);
		synchronized (pathAlns) {
			pathAlns.add(aln);
		}
	}

	private void correctErrors(AssemblyGraph graph, List<ReadAlignment> alignments, AssemblyPath path, String sequenceName, String rawConsensus, int [] sequencePaths) {
		int debugIdx = 0;
		Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
		//TODO: Define better ploidy
		List<CalledGenomicVariant> calledVars = ConsensusBuilderBidirectionalWithPolishing.callVariants(sequenceName, rawConsensus, alignments, 2);
		List<CalledGenomicVariant> filteredVars = new ArrayList<CalledGenomicVariant>();
		for(CalledGenomicVariant call:calledVars) {
			if(!call.isUndecided() && !call.isHomozygousReference()) filteredVars.add(call);
		}
		int pathId = path.getPathId();
		int indexNextActive = 0;
		int correctedErrors = 0;
		int correctedReads = 0;
		
		for(ReadAlignment aln:alignments) {
			int readId = aln.getReadNumber();
			if(sequencePaths[readId]!=pathId) {
				log.warning("Read "+readId+" "+aln.getReadName()+" corrected by path: "+sequencePaths[readId]+" current path: "+pathId);
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
			if(readId == debugIdx) System.out.println("AlignmentBasedErrorCorrection. Correcting read: "+readId+" indel calls: "+calls.size()+" aln: "+aln);
			int nextPos = 0;
			int lastRef = 0;
			for(Map.Entry<Integer, GenomicVariant> entry:calls.entrySet()) {
				int posRead = entry.getKey();
				GenomicVariant indel = entry.getValue();
				if(indel.length()>3) continue;
				if(posRead>nextPos) {
					//Apply homozygous indels not called in this alignment
					for(int j=indexNextActive;j<filteredVars.size();j++) {
						CalledGenomicVariant calledVar = filteredVars.get(j);
						if(readId == debugIdx) System.out.println("AlignmentBasedErrorCorrection. ReadId: "+readId+" Next active var: "+calledVar.getFirst()+" "+calledVar.getLast()+" ref limits: "+lastRef+" "+indel.getFirst()+" read limits: "+nextPos+" "+posRead);
						if(calledVar.getFirst() > lastRef && calledVar.getLast()<indel.getFirst()) {
							if(calledVar.isHeterozygous()) continue;
							int diffLength = calledVar.getAlleles()[0].length()-calledVar.getAlleles()[1].length();
							if(Math.abs(diffLength)>3) continue;
							int readPosStartVar = aln.getAlignedReadPosition(calledVar.getFirst());
							int readPosEndVar = aln.getAlignedReadPosition(calledVar.getLast()); 
							if(readPosStartVar> nextPos && readPosEndVar<posRead && readPosStartVar<readPosEndVar) {
								String currentSegment = alignedRead.substring(readPosStartVar,readPosEndVar+1);
								String correctedAllele = calledVar.getCalledAlleles()[0];
								if(Math.abs(currentSegment.length()-correctedAllele.length())>3) continue;
								if(readId==debugIdx || currentSegment.length()>30) System.out.println("AlignmentBasedErrorCorrection. Adding called indel for read: "+readId+" spanning "+readPosStartVar+" "+readPosEndVar+" segment: "+currentSegment+" alleles var: "+calledVar.getAlleles()[0]+" "+calledVar.getAlleles()[1]+" called allele: "+correctedAllele);
								correctedRead.append(alignedRead.substring(nextPos, readPosStartVar));
								correctedRead.append(correctedAllele);
								nextPos = readPosEndVar+1;
							}
						} else if (calledVar.getLast()>=indel.getFirst()) {
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
					if(readId == debugIdx) System.out.println("AlignmentBasedErrorCorrection. Read id: "+readId+" Correcting indel at: "+indel.getFirst()+" "+indel.getLast()+" length "+indel.length()+" outside active regions");
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
			if(nextPos<aln.getReadLength()) {
				//Apply homozygous indels not called in this alignment
				for(int j=indexNextActive;j<filteredVars.size();j++) {
					CalledGenomicVariant calledVar = filteredVars.get(j);
					if(readId == debugIdx) System.out.println("AlignmentBasedErrorCorrection. ReadId: "+readId+" Next active var: "+calledVar.getFirst()+" "+calledVar.getLast()+" alleles: "+calledVar.getAlleles()[0]+" "+calledVar.getAlleles()[1]+" heterozygous: "+calledVar.isHeterozygous()+" ref limits: "+lastRef+" "+aln.getLast()+" read limits: "+nextPos+" "+aln.getReadLength());
					if(calledVar.getFirst() > lastRef && calledVar.getLast()<aln.getLast()) {
						if(calledVar.isHeterozygous()) continue;
						int diffLength = calledVar.getAlleles()[0].length()-calledVar.getAlleles()[1].length();
						if(Math.abs(diffLength)>3) continue;
						int readPosStartVar = aln.getAlignedReadPosition(calledVar.getFirst());
						int readPosEndVar = aln.getAlignedReadPosition(calledVar.getLast()); 
						if(readPosStartVar> nextPos && readPosEndVar<aln.getReadLength() && readPosStartVar<readPosEndVar) {
							String currentSegment = alignedRead.substring(readPosStartVar,readPosEndVar+1);
							String correctedAllele = calledVar.getCalledAlleles()[0];
							if(Math.abs(currentSegment.length()-correctedAllele.length())>3) continue;
							if(readId==debugIdx || currentSegment.length()>30) System.out.println("AlignmentBasedErrorCorrection. Adding called indel for read: "+readId+" spanning "+readPosStartVar+" "+readPosEndVar+" segment: "+currentSegment+" alleles var: "+calledVar.getAlleles()[0]+" "+calledVar.getAlleles()[1]+" called allele: "+correctedAllele);
							correctedRead.append(alignedRead.substring(nextPos, readPosStartVar));
							correctedRead.append(correctedAllele);
							nextPos = readPosEndVar+1;
						} else {
							if(readId==debugIdx) System.out.println("AlignmentBasedErrorCorrection. Discarded called indel for read: "+readId+" spanning "+readPosStartVar+" "+readPosEndVar+" read limits: "+nextPos+" "+aln.getReadLength());
						}
					} else if (calledVar.getLast()>=aln.getLast()) {
						break;
					}
				}
				correctedRead.append(alignedRead.substring(nextPos));
			}
			//Removing soft clips within the sequence
			 
			/*if(aln.getLast()+aln.getSoftClipEnd()<consensusLength) {
				if(readId == debugIdx || aln.getSoftClipEnd()>1000) System.out.println("Removing "+aln.getSoftClipEnd()+" bp from end of read "+readId+" "+aln);
				correctedRead.delete(correctedRead.length()-aln.getSoftClipEnd(),correctedRead.length());
			}
			if(aln.getFirst()-aln.getSoftClipStart()>0) {
				if(readId == debugIdx || aln.getSoftClipStart()>1000) System.out.println("Removing "+aln.getSoftClipStart()+" bp from start of read "+readId+" "+aln);
				correctedRead.delete(0, aln.getSoftClipStart());
			}*/
			
			//if(readId == debugIdx) System.out.println("AlignmentBasedErrorCorrection. Correcting read: "+readId+" initial read: "+alignedRead+" corrected: "+correctedRead);
			aln.setReadCharacters(null);
			CharSequence correctedReadS;
			if(aln.isNegativeStrand()) {
				correctedReadS = DNAMaskedSequence.getReverseComplement(correctedRead);
			} else {
				correctedReadS = new DNAMaskedSequence(correctedRead);
			}
			QualifiedSequence qseq = graph.getSequence(readId);
			qseq.setCharacters(correctedReadS);
			correctedReads++;
			
		}
		Runtime runtime = Runtime.getRuntime();
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		log.info("IndelErrorsCorrector. Path: "+path.getPathId()+" Corrected reads: "+correctedReads+" correctedErrors: "+correctedErrors+" Memory: "+usedMemory);
	}
}
