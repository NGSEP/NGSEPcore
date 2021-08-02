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
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.math.CountsRankHelper;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.DeBruijnGraphExplorationMiniAssembler;
import ngsep.sequences.DefaultKmersMapImpl;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.KmersMap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;

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
		aligner.setNumThreads(numThreads);
		aligner.setAlignEmbedded(true);
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
		aligner.alignPathReads(graph, path);
		int pathIdx = path.getPathId();
		String sequenceName = path.getSequenceName();
		StringBuilder rawConsensus = aligner.getConsensus();
		List<ReadAlignment> alignments = aligner.getAlignedReads();
		for(ReadAlignment aln:alignments) aln.setSequenceName(sequenceName);
		Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
		
		List<CalledGenomicVariant> variants = callVariants(sequenceName, rawConsensus, alignments, normalPloidy);
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

	public static List<CalledGenomicVariant> callVariants(String sequenceName, CharSequence consensus, List<ReadAlignment> alignments, int normalPloidy) {
		List<GenomicRegion> activeSegments = calculateActiveSegments(sequenceName, alignments);
		List<CalledGenomicVariant> answer=new ArrayList<CalledGenomicVariant>(activeSegments.size());
		System.out.println("Number of active segments "+activeSegments.size());
		int firstIdxAln = 0;
		for(GenomicRegion region:activeSegments) {
			int first = Math.max(1, region.getFirst());
			int last = Math.min(consensus.length(), region.getLast());
			while(firstIdxAln<alignments.size()) {
				ReadAlignment aln = alignments.get(firstIdxAln);
				if(aln.getLast()>=first) break;
				firstIdxAln++;
			}
			String currentConsensus = consensus.subSequence(first-1,last).toString();
			String localConsensus = calculateLocalConsensus(first, last, alignments, firstIdxAln, null);
			String altConsensus = null;
			if(normalPloidy>1 && localConsensus!=null) altConsensus = calculateLocalConsensus(first, last, alignments, firstIdxAln, localConsensus);
			CalledGenomicVariant call = buildCall(sequenceName, first, currentConsensus, localConsensus, altConsensus);
			answer.add(call);
		}
		return answer;
	}
	public static List<GenomicRegion> calculateActiveSegments(String sequenceName, List<ReadAlignment> alignments) {
		//Extract indel calls adding one bp on the sides for insertions
		List<GenomicRegion> rawRegions = new ArrayList<GenomicRegion>();
		for(ReadAlignment aln:alignments) {
			Map<Integer,GenomicVariant> indelCalls = aln.getIndelCalls();
			if(indelCalls==null) continue;
			for(GenomicVariant indelCall:indelCalls.values()) {
				if(indelCall.length()>10) {
					if(indelCall.length()>100) System.out.println("WARN: Long indel from alignment: "+aln + "coordinates: "+indelCall.getFirst()+"-"+indelCall.getLast()+" Ignoring.");
					continue;
				}
				if(indelCall.getLast()-indelCall.getFirst()>1) rawRegions.add(new GenomicRegionImpl(sequenceName, indelCall.getFirst(), indelCall.getLast()));
				else rawRegions.add(new GenomicRegionImpl(sequenceName, indelCall.getFirst()-1, indelCall.getLast()+1));
			}
			
		}
		if(rawRegions.size()<2) return rawRegions;
		//Merge overlapping regions
		Collections.sort(rawRegions,GenomicRegionPositionComparator.getInstance());
		List<GenomicRegion> mergedRegions = new ArrayList<GenomicRegion>();
		GenomicRegionImpl nextRegion = (GenomicRegionImpl)rawRegions.get(0);
		int countSupport = 1;
		for(GenomicRegion rawRegion: rawRegions) {
			if(GenomicRegionSpanComparator.getInstance().span(nextRegion, rawRegion) ) {
				nextRegion.setLast(Math.max(nextRegion.getLast(), rawRegion.getLast()));
				countSupport++;
			} else {
				if(nextRegion.length()>20) System.out.println("Adding long region "+nextRegion.getSequenceName()+":"+nextRegion.getFirst()+"-"+nextRegion.getLast()+" support: "+countSupport);
				if(countSupport>=5) mergedRegions.add(nextRegion);
				nextRegion = (GenomicRegionImpl)rawRegion;
				countSupport=1;
			}
		}
		if(nextRegion.length()>20) System.out.println("Adding long region "+nextRegion.getSequenceName()+":"+nextRegion.getFirst()+"-"+nextRegion.getLast()+" support: "+countSupport);
		if(countSupport>=5) mergedRegions.add(nextRegion);
		return mergedRegions;
	}

	private static String calculateLocalConsensus(int first, int last, List<ReadAlignment> alignments, int firstIdxAln, String consensusAllele) {
		Map<Integer,List<String>> alleleCallsByLength = new HashMap<Integer, List<String>>();
		List<String> allCalls = new ArrayList<String>();
		int count = 0;
		for(int i=firstIdxAln;i<alignments.size();i++) {
			ReadAlignment aln = alignments.get(i);
			if(aln.getFirst()>last) break;
			CharSequence call = aln.getAlleleCall(first, last);
			if(call==null) continue;
			String callStr = call.toString();
			if(consensusAllele!=null && callStr.equals(consensusAllele)) continue;
			count++;
			List<String> lengthCalls = alleleCallsByLength.computeIfAbsent(call.length(), (v)->new ArrayList<String>());
			lengthCalls.add(callStr);
			allCalls.add(callStr);
		}
		//if(first < 10000) System.out.println("Active site: "+first +" "+last+" Alignments: "+count);
		if(count<5) return null;
		List<String> maxLength = null;
		for(List<String> nextList:alleleCallsByLength.values()) {
			if(maxLength==null || maxLength.size()<nextList.size()) {
				maxLength = nextList;
			}
		}
		if(maxLength==null) return null;
		//Double check that the majority length actually has at least half of the total reads
		if(count <10 && maxLength.size()<0.8*count) return null;
		if(2*maxLength.size()<count) {
			if(last-first+1>=8) {
				boolean debug = first ==-1 || first == -2; 
				if(debug) System.out.println("DeBruijn consensus for active site: "+first +" "+last+" calls: "+allCalls);
				String assembly = makeDeBruijnConsensus(last-first+1, allCalls);
				return assembly;
			}
			return null;
		}
		String consensus = HammingSequenceDistanceMeasure.makeHammingConsensus(maxLength);
		//if(first < 10000) System.out.println("Active site: "+first +" "+last+" Majority alleles: "+maxLength+" consensus "+consensus);
		return consensus;
	}

	private static String makeDeBruijnConsensus(int currentLength, List<String> allCalls) {
		KmersMap kmersMap = new DefaultKmersMapImpl();
		CountsRankHelper<String> firstKmerCounts = new CountsRankHelper<String>();
		CountsRankHelper<String> lastKmerCounts = new CountsRankHelper<String>();
		int kmerLength = Math.max(6, currentLength/4);
		kmerLength = Math.min(kmerLength, 15);
		int minCallLength = allCalls.get(0).length();
		int maxCallLength = 0;
		List<Integer> lengths = new ArrayList<Integer>(allCalls.size());
		for(String call:allCalls) {
			minCallLength = Math.min(minCallLength, call.length());
			maxCallLength = Math.max(maxCallLength, call.length());
			lengths.add(call.length());
			int last = call.length()-kmerLength;
			for(int i=0;i<=last;i++) {
				String kmer = call.substring(i,i+kmerLength);
				kmersMap.addOcurrance(kmer);
				if(i==0) firstKmerCounts.add(kmer);
				else if (i==last) lastKmerCounts.add(kmer);
			}
		}
		if(firstKmerCounts.getNumDifferent()==0 || lastKmerCounts.getNumDifferent()==0) return null;
		Collections.sort(lengths);
		int medianCallLength = lengths.get(allCalls.size()/2);
		
		String bestKmerStart = firstKmerCounts.selectBest(1).keySet().iterator().next();
		String bestKmerEnd = lastKmerCounts.selectBest(1).keySet().iterator().next();
		//System.out.println("First kmer: "+bestKmerStart+" lastKmer: "+bestKmerEnd+" total calls: "+allCalls.size()+" median length: "+medianCallLength+" maxlength: "+maxCallLength+" minlength: "+minCallLength);
		DeBruijnGraphExplorationMiniAssembler miniAssembler = new DeBruijnGraphExplorationMiniAssembler(kmersMap, allCalls.size()/3);
		String assembly = miniAssembler.assemble(bestKmerStart, bestKmerEnd, medianCallLength-1, medianCallLength, maxCallLength);
		//System.out.println("Assembly: "+assembly);
		return assembly;
	}

	private static CalledGenomicVariant buildCall(String sequenceName, int first, String currentConsensus, String localConsensus, String altConsensus) {
		List<String> alleles = new ArrayList<String>(2);
		alleles.add(currentConsensus);
		boolean hetero = false;
		if(localConsensus!=null && !localConsensus.equals(currentConsensus)) alleles.add(localConsensus);
		if(altConsensus!=null && !altConsensus.equals(currentConsensus) && !altConsensus.equals(localConsensus)) {
			alleles.add(altConsensus);
		}
		hetero = altConsensus!=null && !altConsensus.equals(localConsensus);
		GenomicVariantImpl variant = new GenomicVariantImpl(sequenceName, first, alleles);
		CalledGenomicVariantImpl call;
		if(alleles.size()==1 || currentConsensus.equals(localConsensus)) call = new CalledGenomicVariantImpl(variant, CalledGenomicVariant.GENOTYPE_HOMOREF); 
		else if (hetero) call = new CalledGenomicVariantImpl(variant, CalledGenomicVariant.GENOTYPE_HETERO);
		else call = new CalledGenomicVariantImpl(variant, CalledGenomicVariant.GENOTYPE_HOMOALT);
		return call;
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
	
	

	private String calculateVariantSegment(ReadAlignment alignment, GenomicVariant indelReadCall, CalledGenomicVariant calledVariant) {
		//Check that indel call is contained
		if(calledVariant.getFirst()>indelReadCall.getFirst()) return null;
		if(calledVariant.getLast()<indelReadCall.getLast()) return null;
		//String readName = alignment.getReadName();
		CharSequence extendedReadCall = alignment.getAlleleCall(calledVariant.getFirst(), calledVariant.getLast());
		if(extendedReadCall==null) return null;
		String [] alleles = calledVariant.getAlleles();
		if(alleles.length<2) return null;
		String majorAllele = calledVariant.isHomozygousReference()?alleles[0]:alleles[1];
		String secondAllele = calledVariant.isHomozygousReference()?alleles[1]:alleles[0];
		if(alleles.length>2) secondAllele = alleles[2];
		
		
		String extendedCallStr = extendedReadCall.toString();
		//if(readName.equals("ref1M_977918_0")) System.out.println("CorrectRead. Read: "+readName+" indel coords: "+indelReadCall.getFirst()+"-"+indelReadCall.getLast()+" variant coords "+calledVariant.getFirst()+"-"+calledVariant.getLast()+" major allele "+majorAllele);
		//if(readName.equals("ref1M_977918_0")) System.out.println("CorrectRead. Read: "+readName+" indel coords: "+indelReadCall.getFirst()+"-"+indelReadCall.getLast()+" variant coords "+calledVariant.getFirst()+"-"+calledVariant.getLast()+" secnd allele "+secondAllele);
		//if(readName.equals("ref1M_977918_0")) System.out.println("CorrectRead. Read: "+readName+" indel coords: "+indelReadCall.getFirst()+"-"+indelReadCall.getLast()+" variant coords "+calledVariant.getFirst()+"-"+calledVariant.getLast()+" extendedcall "+extendedCallStr);
		if(majorAllele.equals(extendedCallStr)) return null;
		if(normalPloidy>1 && secondAllele.equals(extendedCallStr)) return null;
		//TODO: Take into account heterozygosity
		String bestAllele = majorAllele;
		int refOffsetLeft = indelReadCall.getFirst()-calledVariant.getFirst()+1;
		if(refOffsetLeft>=10) return null;
		int bestAlleleStart = refOffsetLeft;
		//int bestAlleleEnd = bestAllele.length()-refOffsetRight;
		if(bestAlleleStart>=bestAllele.length()) return null;
		String answer = bestAllele.substring(bestAlleleStart);
		if(answer.length()>=10) return null;
		//System.out.println("CorrectRead. Read: "+readName+" indel coords: "+indelReadCall.getFirst()+"-"+indelReadCall.getLast()+" variant coords "+calledVariant.getFirst()+"-"+calledVariant.getLast()+" major allele "+majorAllele+" refOffsetLeft: "+refOffsetLeft+" best allele start: "+bestAlleleStart+" extended call "+extendedCallStr);
		//if(bestAlleleStart<0) throw new RuntimeException("CorrectRead. Read: "+readName+" indel coords: "+indelReadCall.getFirst()+"-"+indelReadCall.getLast()+" variant coords "+calledVariant.getFirst()+"-"+calledVariant.getLast()+" major allele "+majorAllele+" refOffsetLeft: "+refOffsetLeft+" best allele start: "+bestAlleleStart+" extended call "+extendedCallStr);
		//if(readName.equals("ref1M_977918_0")) System.out.println("CorrectRead. Read: "+readName+" indel coords: "+indelReadCall.getFirst()+"-"+indelReadCall.getLast()+" consensus segment "+answer);
		return answer;
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
