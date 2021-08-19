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
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.discovery.AlignmentsPileupGenerator;
import ngsep.discovery.PileupListener;
import ngsep.discovery.PileupRecord;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class ConsensusBuilderBidirectionalWithPolishing implements ConsensusBuilder {
	
	private Logger log = Logger.getLogger(ConsensusBuilderBidirectionalWithPolishing.class.getName());
	public static final int DEF_NUM_THREADS = 1;
	
	private String sequenceNamePrefix = "Contig";
	
	private short normalPloidy = 1;
	
	private int numThreads = DEF_NUM_THREADS;
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	
	public String getSequenceNamePrefix() {
		return sequenceNamePrefix;
	}

	public void setSequenceNamePrefix(String sequenceNamePrefix) {
		this.sequenceNamePrefix = sequenceNamePrefix;
	}
	
	public int getNumThreads() {
		return numThreads;
	}

	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}

	@Override
	public List<QualifiedSequence> makeConsensus(AssemblyGraph graph) 
	{
		AssemblyPathReadsAligner aligner = new AssemblyPathReadsAligner();
		aligner.setLog(log);
		List<QualifiedSequence> consensusList = new ArrayList<QualifiedSequence>();
		List<AssemblyPath> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			AssemblyPath path = paths.get(i);
			String sequenceName = ""+sequenceNamePrefix+"_"+(i+1);
			path.setPathId(i+1);
			path.setSequenceName(sequenceName);
			CharSequence consensusSequence = makeConsensus (graph, path, aligner);
			consensusList.add(new QualifiedSequence(sequenceName,consensusSequence));	
		}
		return consensusList;
	}
	
	
	private CharSequence makeConsensus(AssemblyGraph graph, AssemblyPath path, AssemblyPathReadsAligner aligner) {
		aligner.calculateConsensus(path);
		List<ReadAlignment> alignments = aligner.alignPathReads(graph, path, numThreads);
		int pathIdx = path.getPathId();
		String sequenceName = path.getSequenceName();
		StringBuilder rawConsensus = new StringBuilder(path.getConsensus());
		for(ReadAlignment aln:alignments) aln.setSequenceName(sequenceName);
		Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
		
		List<CalledGenomicVariant> variants = aligner.callIndels(path, alignments, normalPloidy);
		log.info("Path "+pathIdx+" "+sequenceName+" Identified "+variants.size()+" total variants from read alignments");
		//Identify and correct SNV errors first
		correctSNVErrors(sequenceName, rawConsensus, alignments, variants);
		return applyVariants(rawConsensus, variants);
	}

	/*private boolean containsLargeIndels(ReadAlignment alnRead) {
		Map<Integer,GenomicVariant> indelCalls = alnRead.getIndelCalls();
		if(indelCalls==null) return false;
		for(GenomicVariant call:indelCalls.values()) if(call.length()>10) return true;
		return false;
	}*/

	public void printAllOverlappingSeqs(AssemblyGraph graph, List<AssemblyEdge> path, int pathPos, AssemblyVertex vertexPreviousEdge) {
		System.out.println("Vertex to check: "+vertexPreviousEdge);
		for(int j = pathPos; j < path.size(); j++) {
			AssemblyEdge edge = path.get(j);
			if(edge.isSameSequenceEdge()) continue;
			AssemblyEdge alt1 = graph.getEdge(vertexPreviousEdge, edge.getVertex1());
			AssemblyEdge alt2 = graph.getEdge(vertexPreviousEdge, edge.getVertex2());
			System.out.println("Edge 1: "+alt1);
			System.out.println("Edge 2: "+alt2);
		}		
	}

	

	private void correctSNVErrors(String sequenceName, StringBuilder consensus, List<ReadAlignment> alignments, List<CalledGenomicVariant> variants) {
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		generator.setLog(log);
		QualifiedSequenceList metadata = new QualifiedSequenceList();
		metadata.add(new QualifiedSequence(sequenceName,consensus.length()));
		generator.setSequencesMetadata(metadata);
		generator.setMaxAlnsPerStartPos(0);
		SimpleSNVErrorCorrectorPileupListener snvsCorrectorListener = new SimpleSNVErrorCorrectorPileupListener(consensus, variants);
		generator.addListener(snvsCorrectorListener);
		
		int count = 0;
		for(ReadAlignment aln:alignments) {
			generator.processAlignment(aln);
			count++;
			if(count%1000==0) log.info("Sequence: "+sequenceName+". Corrected SNVs from "+count+" alignments"); 
		}
		generator.notifyEndOfAlignments();
	}

	private CharSequence applyVariants(StringBuilder consensus, List<CalledGenomicVariant> variants) {
		StringBuilder polishedConsensus = new StringBuilder(consensus.length());
		int l = consensus.length();
		int nextPos = 1;
		int appliedVariants = 0;
		for(CalledGenomicVariant call:variants) {
			if(call.isUndecided()) continue;
			if(call.isHomozygousReference()) continue;
			String [] alleles = call.getAlleles();
			if(alleles.length==1) continue;
			appliedVariants++;
			if(nextPos<call.getFirst()) {
				//Fill haplotypes with non variant segment
				String segment = consensus.substring(nextPos-1, call.getFirst()-1);
				polishedConsensus.append(segment);
			}
			//The reconstructed consensus allele is the first alternative allele in the call
			polishedConsensus.append(alleles[1]);
			nextPos = call.getLast()+1;
		}
		if(nextPos<=l) {
			//Consensus end
			CharSequence nonVarLast = consensus.substring(nextPos-1);
			polishedConsensus.append(nonVarLast);
		}
		log.info("Applied "+appliedVariants+" variants");
		return new DNAMaskedSequence(polishedConsensus.toString());
	}
}

class SimpleSNVErrorCorrectorPileupListener implements PileupListener {

	private StringBuilder consensus;
	private List<CalledGenomicVariant> indelRegions;
	private int nextIndelPos = 0;

	public SimpleSNVErrorCorrectorPileupListener(StringBuilder consensus, List<CalledGenomicVariant> indelRegions) {
		super();
		this.consensus = consensus;
		this.indelRegions = indelRegions;
	}
	
	@Override
	public void onPileup(PileupRecord pileup) {
		int pos = pileup.getPosition();
		//Check if pileup is located within an indel region
		while(nextIndelPos<indelRegions.size()) {
			CalledGenomicVariant region = indelRegions.get(nextIndelPos);
			if(region.getFirst()<=pos && pos<=region.getLast()) return;
			else if (pos<region.getFirst()) break;
			nextIndelPos++;
		}
		List<ReadAlignment> alns = pileup.getAlignments();
		//Index alignments per nucleotide call
		int n = DNASequence.BASES_STRING.length();
		Map<Character,List<ReadAlignment>> alnsPerNucleotide = new HashMap<Character, List<ReadAlignment>>(n);
		for(int i=0;i<n;i++) {
			alnsPerNucleotide.put(DNASequence.BASES_STRING.charAt(i), new ArrayList<ReadAlignment>(alns.size()));
		}
		for(ReadAlignment aln:alns) {
			CharSequence call = aln.getAlleleCall(pos);
			if(call == null) continue;
			char c = call.charAt(0);
			List<ReadAlignment> alnsAllele = alnsPerNucleotide.get(c);
			if(alnsAllele==null) continue;
			alnsAllele.add(aln);
		}
		//Extract counts from map of allele calls
		int [] acgtCounts = new int [n];
		for(int i=0;i<n;i++) {
			char c = DNASequence.BASES_STRING.charAt(i);
			acgtCounts[i] = alnsPerNucleotide.get(c).size();
		}
		int maxIdx = NumberArrays.getIndexMaximum(acgtCounts);
		int maxCount = acgtCounts[maxIdx];
		char maxBP = DNASequence.BASES_STRING.charAt(maxIdx);
		int consensusPos = pileup.getPosition()-1;
		char refBase = consensus.charAt(consensusPos);
		int refIdx = DNASequence.BASES_STRING.indexOf(refBase);
		int refCount = (refIdx>=0?acgtCounts[refIdx]:0);
		if(maxIdx!=refIdx && maxCount>refCount) {
			consensus.setCharAt(consensusPos, maxBP);
		} 
		if (maxCount<10) return;
	}

	@Override
	public void onSequenceStart(QualifiedSequence sequence) {
	}

	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {
	}
	
}
