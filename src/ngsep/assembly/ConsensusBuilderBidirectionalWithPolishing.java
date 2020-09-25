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

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.discovery.AlignmentsPileupGenerator;
import ngsep.discovery.CountsHelper;
import ngsep.discovery.PileupAlleleCall;
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
import ngsep.sequences.KmersExtractor;
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
	private static final int TIMEOUT_SECONDS = 30;
	
	private String correctedReadsFile = null;
	private int kmerLength = KmersExtractor.DEF_KMER_LENGTH;
	
	private PrintStream outCorrectedReads = null;
	
	private int numThreads = DEF_NUM_THREADS;
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	
	public int getNumThreads() {
		return numThreads;
	}

	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	
	public String getCorrectedReadsFile() {
		return correctedReadsFile;
	}

	public void setCorrectedReadsFile(String correctedReadsFile) {
		this.correctedReadsFile = correctedReadsFile;
	}

	@Override
	public List<QualifiedSequence> makeConsensus(AssemblyGraph graph) 
	{
		if(correctedReadsFile==null) return makeConsensus2(graph);
		try (OutputStream os = new GZIPOutputStream(new FileOutputStream(correctedReadsFile));
			 PrintStream out = new PrintStream(os)) {
			this.outCorrectedReads = out;
			return makeConsensus2(graph);
		} catch (IOException e) {
			e.printStackTrace();
			log.warning("Error opening file for corrected reads "+e.getMessage());
			return null;
		}
	}
	
	private List<QualifiedSequence> makeConsensus2(AssemblyGraph graph) {
		ThreadPoolExecutor poolConsensus = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		List<QualifiedSequence> consensusList = new ArrayList<QualifiedSequence>();
		Map<String, BuildConsensusTask> tasksMap = new LinkedHashMap<String, BuildConsensusTask>();
		List<List<AssemblyEdge>> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			List<AssemblyEdge> path = paths.get(i);
			String sequenceName = "Contig_"+(i+1);
			if(numThreads==1) {
				CharSequence consensusSequence = makeConsensus (graph, path, i, sequenceName);
				consensusList.add(new QualifiedSequence(sequenceName,consensusSequence));
			} else {
				BuildConsensusTask task = new BuildConsensusTask(this, graph, path, i, sequenceName);
				poolConsensus.execute(task);
				tasksMap.put(sequenceName, task);
			}	
		}
		int finishTime = 10*graph.getNumSequences();
		poolConsensus.shutdown();
		try {
			poolConsensus.awaitTermination(finishTime, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
    	if(!poolConsensus.isShutdown()) {
			throw new RuntimeException("The ThreadPoolExecutor was not shutdown after an await Termination call");
		}
    	for(String sequenceName:tasksMap.keySet()) {
    		BuildConsensusTask task = tasksMap.get(sequenceName);
    		CharSequence seq = task.getConsensusSequence();
    		if(seq==null) {
    			log.warning("Null consensus sequence for contig: "+sequenceName);
    			continue;
    		}
    		consensusList.add(new QualifiedSequence(sequenceName,seq));
    	}
		return consensusList;
	}
	
	
	CharSequence makeConsensus(AssemblyGraph graph, List<AssemblyEdge> path, int sequenceIdx, String sequenceName ) {
		StringBuilder rawConsensus = new StringBuilder();
		AssemblyVertex lastVertex = null;
		MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm();
		List<ReadAlignment> alignments = new ArrayList<ReadAlignment>();
		int totalReads = 0;
		int unalignedReads = 0;
		String pathS = "";
		if(path.size()==1) {
			rawConsensus.append(path.get(0).getVertex1().getRead());
			return rawConsensus;
		}
		ReadAlignment lastPartialAln = null;
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
				pathS = pathS.concat(vertexNextEdge.getUniqueNumber() + ",");
				// Augment consensus with the next path read
				CharSequence nextPathSequence = vertexNextEdge.getRead().getCharacters();
				boolean reverse = !vertexNextEdge.isStart();
				if(reverse) nextPathSequence = DNAMaskedSequence.getReverseComplement(nextPathSequence);
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) printAllOverlappingSeqs(graph,path,j,vertexPreviousEdge);
				
				//int startSuffix = edge.getOverlap();
				Map<Long, Integer> uniqueKmersSubject = KmersExtractor.extractLocallyUniqueKmerCodes(rawConsensus,kmerLength, Math.max(0, rawConsensus.length()-nextPathSequence.length()),rawConsensus.length());
				ReadAlignment alnRead = ConsensusBuilderBidirectionalSimple.alignRead(aligner, sequenceIdx, rawConsensus, nextPathSequence, uniqueKmersSubject, 0.5);
				int startSuffix;
				if(alnRead!=null) {
					alnRead.setSequenceName(sequenceName);
					alnRead.setReadName(vertexNextEdge.getRead().getName());
					int posAlnRead = nextPathSequence.length()-1-alnRead.getSoftClipEnd();
					int lastPosSubject = alnRead.getReferencePositionAlignedRead(posAlnRead);
					//Just in case cycle but if the read aligns this should not enter
					while(posAlnRead>0 && lastPosSubject<0) {
						posAlnRead--;
						lastPosSubject = alnRead.getReferencePositionAlignedRead(posAlnRead);
					}
					if(lastPosSubject>=0) {
						startSuffix = posAlnRead + (rawConsensus.length()-lastPosSubject+1);
					} else {
						startSuffix = edge.getOverlap();
					}
					//System.out.println("Calculated overlap from alignment: "+startSuffix+" alignment: "+alnRead+" edge: "+edge );
				} else {
					startSuffix = edge.getOverlap();
				}
				if(startSuffix<nextPathSequence.length()) {
					String remainingSegment = nextPathSequence.subSequence(startSuffix, nextPathSequence.length()).toString();
					//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) System.out.println("Consensus length: "+rawConsensus.length()+" Vertex: "+vertexNextEdge.getUniqueNumber()+" read length: "+nextPathSequence.length()+" overlap: "+edge.getOverlap()+" remaining: "+remainingSegment.length());
					rawConsensus.append(remainingSegment.toUpperCase());
				}
				lastPartialAln = alnRead;
			}
			if(vertexPreviousEdge.getRead()==vertexNextEdge.getRead()) {
				//Align to consensus next path read and its embedded sequences
				QualifiedSequence read = vertexPreviousEdge.getRead(); 
				CharSequence seq = read.getCharacters();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
				Map<Long, Integer> uniqueKmersSubject = KmersExtractor.extractLocallyUniqueKmerCodes(rawConsensus,kmerLength,Math.max(0, rawConsensus.length()-seq.length()),rawConsensus.length());
				totalReads++;
				ReadAlignment alnRead = ConsensusBuilderBidirectionalSimple.alignRead(aligner, sequenceIdx, rawConsensus, seq, uniqueKmersSubject, 0.5);
				if (alnRead!=null) {
					alnRead.setSequenceName(sequenceName);
					alnRead.setReadName(read.getName());
					alignments.add(alnRead);
					if (lastPartialAln == null) {
						log.warning("Consensus backbone read did not align to last consensus" );
					}
					if(alnRead.getSoftClipEnd()>0) {
						log.warning("Weird alignment of consensus backbone read. Partial alignment: "+lastPartialAln+" soft clipped sequence: "+lastPartialAln.getReadCharacters().subSequence(lastPartialAln.getReadLength()-lastPartialAln.getSoftClipEnd()-5, lastPartialAln.getReadLength()) );
						log.warning("Alignment to enlarged consensus: "+alnRead+" consensus end: "+rawConsensus.substring(rawConsensus.length()-lastPartialAln.getSoftClipEnd()-5));
					}
				}
				else unalignedReads++;
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) System.out.println("Consensus length: "+rawConsensus.length()+" Vertex: "+vertexNextEdge.getUniqueNumber()+" sequence: "+read.length()+" alignment: "+alnRead);
				if (totalReads%100==0) log.info("Sequence "+sequenceName+". Aligning. Processed reads: "+totalReads+" alignments: "+alignments.size()+" unaligned: "+unalignedReads);
				
				List<AssemblyEmbedded> embeddedList = graph.getAllEmbedded(vertexPreviousEdge.getSequenceIndex());
				//List<AssemblyEmbedded> embeddedList = graph.getEmbeddedByHostId(vertexPreviousEdge.getSequenceIndex());
				for(AssemblyEmbedded embedded:embeddedList) {
					QualifiedSequence embeddedRead = embedded.getRead(); 
					CharSequence embeddedSeq = embeddedRead.getCharacters();
					boolean reverseE = (reverse!=embedded.isReverse());
					if(reverseE) embeddedSeq = DNAMaskedSequence.getReverseComplement(embeddedSeq);
					totalReads++;
					ReadAlignment alnEmbedded = ConsensusBuilderBidirectionalSimple.alignRead(aligner,sequenceIdx, rawConsensus, embeddedSeq, uniqueKmersSubject, 0.5);
					if(alnEmbedded!=null) {
						alnEmbedded.setSequenceName(sequenceName);
						alnEmbedded.setReadName(embeddedRead.getName());
						alignments.add(alnEmbedded);
					}
					else unalignedReads++;
					if (totalReads%100==0) log.info("Sequence "+sequenceName+". Aligning. Processed reads: "+totalReads+" alignments: "+alignments.size()+" unaligned: "+unalignedReads);
				}
			}
			lastVertex = vertexNextEdge;
		}
		log.info("Path "+sequenceIdx+". Total reads: "+totalReads+" alignments: "+alignments.size()+" unaligned: "+unalignedReads);
		log.info("Path "+sequenceIdx+". Path: "+pathS);
		
		Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
		//Identify and correct SNV errors first
		correctSNVErrors(sequenceName, rawConsensus, alignments);
		String consensus = rawConsensus.toString();
		List<CalledGenomicVariant> variants = callVariants(sequenceName, consensus, alignments);
		log.info("Path "+sequenceIdx+". Identified "+variants.size()+" total variants from read alignments");
		if(outCorrectedReads!=null) {
			for(ReadAlignment aln:alignments) correctRead(consensus, aln, variants);
		}
		return applyVariants(consensus, variants);
	}

	/*private boolean containsLargeIndels(ReadAlignment alnRead) {
		Map<Integer,GenomicVariant> indelCalls = alnRead.getIndelCalls();
		if(indelCalls==null) return false;
		for(GenomicVariant call:indelCalls.values()) if(call.length()>10) return true;
		return false;
	}*/

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

	private List<CalledGenomicVariant> callVariants(String sequenceName, String consensus, List<ReadAlignment> alignments) {
		List<CalledGenomicVariant> answer=new ArrayList<CalledGenomicVariant>();
		List<GenomicRegion> activeSegments = calculateActiveSegments(sequenceName, alignments);
		log.info("Number of active segments "+activeSegments.size());
		int firstIdxAln = 0;
		for(GenomicRegion region:activeSegments) {
			int first = Math.max(1, region.getFirst());
			int last = Math.min(consensus.length(), region.getLast());
			while(firstIdxAln<alignments.size()) {
				ReadAlignment aln = alignments.get(firstIdxAln);
				if(aln.getLast()>=first) break;
				firstIdxAln++;
			}
			String currentConsensus = consensus.substring(first-1,last);
			String localConsensus = calculateLocalConsensus(first, last, alignments, firstIdxAln);
			if(localConsensus!=null && !localConsensus.equals(currentConsensus)) {
				CalledGenomicVariant call = buildCall(sequenceName, first, currentConsensus, localConsensus);
				/*if(Math.abs(currentConsensus.length()-localConsensus.length())>=5) {
					System.out.println("Next call: "+sequenceName+" "+call.getFirst()+" "+call.getLast()+" length: "+call.length());
					System.out.println(currentConsensus);
					System.out.println(localConsensus);
				}*/
				answer.add(call);
			}
		}
		return answer;
	}
	private List<GenomicRegion> calculateActiveSegments(String sequenceName, List<ReadAlignment> alignments) {
		List<GenomicRegion> rawRegions = new ArrayList<GenomicRegion>();
		for(ReadAlignment aln:alignments) {
			Map<Integer,GenomicVariant> indelCalls = aln.getIndelCalls();
			if(indelCalls==null) continue;
			for(GenomicVariant indelCall:indelCalls.values()) {
				if(indelCall.getLast()-indelCall.getFirst()>1) rawRegions.add(new GenomicRegionImpl(sequenceName, indelCall.getFirst(), indelCall.getLast()));
				else rawRegions.add(new GenomicRegionImpl(sequenceName, indelCall.getFirst()-1, indelCall.getLast()+1));
			}
			
		}
		if(rawRegions.size()<2) return rawRegions;
		Collections.sort(rawRegions,GenomicRegionPositionComparator.getInstance());
		List<GenomicRegion> mergedRegions = new ArrayList<GenomicRegion>();
		GenomicRegionImpl nextRegion = (GenomicRegionImpl)rawRegions.get(0);
		mergedRegions.add(nextRegion);
		for(GenomicRegion rawRegion: rawRegions) {
			if(GenomicRegionSpanComparator.getInstance().span(nextRegion, rawRegion) ) {
				nextRegion.setLast(Math.max(nextRegion.getLast(), rawRegion.getLast()));
			} else {
				nextRegion = (GenomicRegionImpl)rawRegion;
				mergedRegions.add(nextRegion);
			}
		}
		return mergedRegions;
	}

	private String calculateLocalConsensus(int first, int last, List<ReadAlignment> alignments, int firstIdxAln) {
		Map<Integer,List<String>> alleleCallsByLength = new HashMap<Integer, List<String>>();
		List<String> allCalls = new ArrayList<String>();
		int count = 0;
		for(int i=firstIdxAln;i<alignments.size();i++) {
			ReadAlignment aln = alignments.get(i);
			if(aln.getFirst()>last) break;
			CharSequence call = aln.getAlleleCall(first, last);
			if(call==null) continue;
			count++;
			List<String> lengthCalls = alleleCallsByLength.computeIfAbsent(call.length(), (v)->new ArrayList<String>());
			String callStr = call.toString();
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

	private String makeDeBruijnConsensus(int currentLength, List<String> allCalls) {
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

	private CalledGenomicVariant buildCall(String sequenceName, int first, String currentConsensus, String localConsensus) {
		List<String> alleles = new ArrayList<String>(2);
		alleles.add(currentConsensus);
		alleles.add(localConsensus);
		GenomicVariantImpl variant = new GenomicVariantImpl(sequenceName, first, alleles);
		CalledGenomicVariantImpl call = new CalledGenomicVariantImpl(variant, CalledGenomicVariant.GENOTYPE_HOMOALT);
		return call;
	}

	private void correctSNVErrors(String sequenceName, StringBuilder consensus, List<ReadAlignment> alignments) {
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		generator.setLog(log);
		QualifiedSequenceList metadata = new QualifiedSequenceList();
		metadata.add(new QualifiedSequence(sequenceName,consensus.length()));
		generator.setSequencesMetadata(metadata);
		generator.setMaxAlnsPerStartPos(0);
		SimpleSNVErrorCorrectorPileupListener snvsCorrectorListener = new SimpleSNVErrorCorrectorPileupListener(consensus);
		generator.addListener(snvsCorrectorListener);
		
		int count = 0;
		for(ReadAlignment aln:alignments) {
			generator.processAlignment(aln);
			count++;
			if(count%100==0) log.info("Sequence: "+sequenceName+". Corrected SNVs from "+count+" alignments"); 
		}
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
			//Consensus end
			CharSequence nonVarLast = consensus.substring(nextPos-1);
			polishedConsensus.append(nonVarLast);
		}
		log.info("Applied "+appliedVariants+" homozygous alternative variants");
		return new DNAMaskedSequence(polishedConsensus.toString());
	}
	
	private void correctRead(String consensus, ReadAlignment alignment, List<CalledGenomicVariant> variants) {
		StringBuilder correctedRead = new StringBuilder();
		String readName = alignment.getReadName();
		CharSequence sequence = alignment.getReadCharacters();
		Map<Integer,GenomicVariant> indelCalls = alignment.getIndelCalls();
		if(indelCalls==null || indelCalls.size()==0) {
			printRead (readName, sequence);
			return;
		}
		int nextVarIdx = 0;
		int nextReadPos = 0;
		for(int refpos:indelCalls.keySet()) {
			GenomicVariant indelReadCall = indelCalls.get(refpos);
			int indelReadPos = alignment.getAlignedReadPosition(refpos);
			if(indelReadPos<nextReadPos) {
				log.warning("Inconsistency correcting errors of "+readName+". Next indel pos: "+indelReadPos+" next pos: "+nextReadPos);
				printRead (readName, sequence);
				return;
			}
			correctedRead.append(sequence.subSequence(nextReadPos, indelReadPos+1));
			nextReadPos = indelReadPos+1;
			boolean correctIndel = true;
			//Check if called indel
			while(nextVarIdx<variants.size()) {
				CalledGenomicVariant calledVariant = variants.get(nextVarIdx);
				if(GenomicRegionSpanComparator.getInstance().span(calledVariant, indelReadCall.getFirst(), indelReadCall.getLast())) {
					correctIndel=false;
					break;
				} else if (indelReadCall.getLast() < calledVariant.getFirst()) {
					break;
				}
				nextVarIdx++;
			}
			//Process indel
			if(!correctIndel) continue;
			if(indelReadCall.getLast()-refpos==1) {
				//Ignore false insertion
				nextReadPos+=indelReadCall.length();
			} else {
				//add reference sequence for false deletion. 
				//Genomic coordinates in this case correspond to string indexes because the deletion coordinates include flanking basepairs
				correctedRead.append(consensus.subSequence(refpos, indelReadCall.getLast()-1));
				
			}
		}
		correctedRead.append(sequence.subSequence(nextReadPos,sequence.length()));
		printRead (readName, correctedRead.toString());
	}

	private void printRead(String readName, CharSequence sequence) {
		synchronized (outCorrectedReads) {
			outCorrectedReads.println(">"+readName);
			outCorrectedReads.println(sequence);
		}
	}
}
class BuildConsensusTask implements Runnable {
	private ConsensusBuilderBidirectionalWithPolishing parent;
	private AssemblyGraph graph;
	private List<AssemblyEdge> path;
	private int sequenceIdx;
	private String sequenceName;
	private CharSequence consensusSequence;
	public BuildConsensusTask(ConsensusBuilderBidirectionalWithPolishing parent, AssemblyGraph graph, List<AssemblyEdge> path, int sequenceIdx, String sequenceName) {
		super();
		this.parent = parent;
		this.graph = graph;
		this.path = path;
		this.sequenceIdx = sequenceIdx;
		this.sequenceName = sequenceName;
	}
	public CharSequence getConsensusSequence() {
		return consensusSequence;
	}
	@Override
	public void run() {
		consensusSequence = parent.makeConsensus(graph,path, sequenceIdx, sequenceName);
	}
}

class SimpleSNVErrorCorrectorPileupListener implements PileupListener {

	private StringBuilder consensus;

	public SimpleSNVErrorCorrectorPileupListener(StringBuilder consensus) {
		super();
		this.consensus = consensus;
	}
	
	@Override
	public void onPileup(PileupRecord pileup) {
		List<PileupAlleleCall> calls = pileup.getAlleleCalls(1);
		CountsHelper helper = CountsHelper.calculateCountsSNV(calls, (byte)30, 0.0);
		int [] acgtCounts = helper.getCounts();
		int maxIdx = NumberArrays.getIndexMaximum(acgtCounts);
		int consensusPos = pileup.getPosition()-1;
		char refBase = consensus.charAt(consensusPos);
		int refIdx = DNASequence.BASES_STRING.indexOf(refBase);
		if(maxIdx!=refIdx) {
			consensus.setCharAt(consensusPos, DNASequence.BASES_STRING.charAt(maxIdx));
		}
	}

	@Override
	public void onSequenceStart(QualifiedSequence sequence) {
	}

	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {
	}
	
}
