package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.PairwiseAlignerDynamicKmers;
import ngsep.alignments.UngappedSearchHitsCluster;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.alignments.LongReadsUngappedSearchHitsClusterAligner;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.ReadAlignmentPositionComparator;
import ngsep.alignments.ReadsAligner;
import ngsep.main.ThreadPoolManager;
import ngsep.math.CountsRankHelper;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.DeBruijnGraphExplorationMiniAssembler;
import ngsep.sequences.DefaultKmersMapImpl;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.ShortKmerCodesTable;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;

public class AssemblyPathReadsAligner {
	private Logger log = Logger.getLogger(AssemblyPathReadsAligner.class.getName());
	private static Runtime runtime = Runtime.getRuntime();
	public static final int KMER_LENGTH_LOCAL_ALN = 25;
	private boolean haploid = true;
	private boolean buildUnalignedReadRecords = false;
	private ShortKmerCodesTable kmerCodesTable = new ShortKmerCodesTable(KMER_LENGTH_LOCAL_ALN, 40,0);
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
		//factory.setLog(log);
	}
	
	
	
	public boolean isHaploid() {
		return haploid;
	}
	public void setHaploid(boolean haploid) {
		this.haploid = haploid;
	}
	
	
	public boolean isBuildUnalignedReadRecords() {
		return buildUnalignedReadRecords;
	}
	public void setBuildUnalignedReadRecords(boolean buildUnalignedReadRecords) {
		this.buildUnalignedReadRecords = buildUnalignedReadRecords;
	}
	public void calculateConsensus(AssemblyPath path) {
		int debugIdx = -1;
		int n = path.getPathLength();
		int pathIdx = path.getPathId();
		long usedMemory = (runtime.totalMemory()-runtime.freeMemory())/1000000;
		log.info("Building consensus for path "+pathIdx+" "+path.getSequenceName()+" with length: "+n+" Memory (Mbp): "+usedMemory);
		List<AssemblyEdge> edges = path.getEdges();
		Map<Integer,Integer> pathVerticesEnds = new HashMap<Integer, Integer>();
		StringBuilder rawConsensus = new StringBuilder();
		AssemblyVertex lastVertex = path.getVertexLeft();
		//Build consensus first
		LongReadsUngappedSearchHitsClusterAligner aligner = new LongReadsUngappedSearchHitsClusterAligner(LongReadsUngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		for(int j = 0; j < n; j++) {
			AssemblyEdge edge = edges.get(j);
			AssemblyVertex nextVertex = edge.getConnectingVertex(lastVertex);
			if(nextVertex == null) throw new RuntimeException("Inconsistency found in path. Previous vertex: "+lastVertex+" edge: "+edge);
			if(j == 0) {
				CharSequence seq = lastVertex.getRead().getCharacters();
				boolean reverse = !lastVertex.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
				pathVerticesEnds.put(lastVertex.getSequenceIndex(), seq.length());
				rawConsensus.append(seq);
				if(pathIdx == debugIdx) System.err.println("Added sequence: "+lastVertex.getRead().getName()+" reverse: "+reverse+" length: "+seq.length());
			} else if(!edge.isSameSequenceEdge()) {
				// Augment consensus with the next path read
				QualifiedSequence nextRead = nextVertex.getRead();
				CharSequence nextPathSequence = nextRead.getCharacters();
				boolean reverse = !nextVertex.isStart();
				if(reverse) nextPathSequence = DNAMaskedSequence.getReverseComplement(nextPathSequence);
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) printAllOverlappingSeqs(graph,path,j,vertexPreviousEdge);
				if(pathIdx == debugIdx) System.err.println("Aligning next path read "+nextRead.getName()+". length1: "+nextRead.getLength()+" length2: "+nextPathSequence.length()+" Reverse "+reverse+ " edge: "+edge);
				int startSuffixConsensus = Math.max(0, rawConsensus.length()-500);
				Map<Integer, Long> kmersSubject = KmersExtractor.extractDNAKmerCodes(rawConsensus, KMER_LENGTH_LOCAL_ALN, startSuffixConsensus,rawConsensus.length());
				int endSegmentQuery = Math.min(nextPathSequence.length(), edge.getOverlap()+10);
				int startSegmentQuery = Math.max(0, endSegmentQuery-520);
				String segmentQuery = nextPathSequence.subSequence(startSegmentQuery, endSegmentQuery).toString();
				ReadAlignment alnRead = alignRead(aligner, pathIdx, rawConsensus, segmentQuery, kmersSubject);
				int startRemove = -1;
				int startSuffixQuery;
				if(alnRead!=null) {
					alnRead.setReadName(nextVertex.getRead().getName());
					int posAlnRead = segmentQuery.length()-1-alnRead.getSoftClipEnd();
					int lastPosSubject = alnRead.getReferencePositionAlignedRead(posAlnRead);
					int tailSubject = rawConsensus.length()-lastPosSubject-1;
					if(pathIdx == debugIdx) System.err.println("Sequence length: "+nextPathSequence.length()+" subject length: "+rawConsensus.length()+" Soft clip end: "+alnRead.getSoftClipEnd()+" pos aln: "+posAlnRead+" pos subject: "+lastPosSubject+" tail: "+tailSubject+" aln: "+alnRead);
					//if(alnRead.getSoftClipEnd()>0 && lastPosSubject>=0 && tailSubject>50) System.err.println("Large subject tail not aligned. Tail length "+tailSubject+" new sequence suffix: "+(lastPosSubject+1)+" Next path alignment: "+alnRead+" Tail of subject: "+rawConsensus.substring(lastPosSubject+1)+" end read: "+nextPathSequence.subSequence(posAlnRead+1, nextPathSequence.length()));
					//Just in case cycle but if the read aligns this should not enter
					while(posAlnRead>0 && lastPosSubject<0) {
						if(pathIdx == debugIdx && j<10) System.err.println("Negative pos subject: "+lastPosSubject+" for read position: "+posAlnRead+" read length: "+nextPathSequence.length());
						posAlnRead--;
						lastPosSubject = alnRead.getReferencePositionAlignedRead(posAlnRead);
						tailSubject = rawConsensus.length()-lastPosSubject-1;
					}
					if(lastPosSubject>=0) {
						startSuffixQuery = startSegmentQuery + posAlnRead + 1;
						int suffixLength = nextPathSequence.length()-startSuffixQuery;
						if(tailSubject>0 && 2*tailSubject<suffixLength) {
							//Replace tail with new read end
							startRemove = lastPosSubject+1;
						} else if (tailSubject<suffixLength) {
							startSuffixQuery+=tailSubject;
						} else {
							startSuffixQuery = nextPathSequence.length();
						}
					} else {
						startSuffixQuery = edge.getOverlap();
					}
					if(pathIdx == debugIdx && j<10) System.err.println("Calculated start new suffix from alignment: "+startSuffixQuery );
				} else {
					if(pathIdx == debugIdx && j<10) System.err.println("WARN: Consensus backbone read "+nextVertex.getRead().getName()+" did not align to last consensus. Using overlap: "+edge.getOverlap());
					startSuffixQuery = edge.getOverlap();
					
				}
				if(pathIdx == debugIdx && j<10) System.err.println("Start suffix: "+startSuffixQuery+" next seq len: "+nextPathSequence.length()+" start remove previous: "+startRemove);
				if(startSuffixQuery<nextPathSequence.length()) {
					String remainingSegment = nextPathSequence.subSequence(startSuffixQuery, nextPathSequence.length()).toString();
					if(pathIdx == debugIdx) System.err.println("Enlarging consensus. Current length: "+rawConsensus.length()+" remove from: "+startRemove+" SeqName: "+nextRead.getName()+" length: "+nextPathSequence.length()+" Start segment: "+startSuffixQuery+" length segment: "+remainingSegment.length());
					if(startRemove>0) rawConsensus.delete(startRemove, rawConsensus.length());
					rawConsensus.append(remainingSegment.toUpperCase());
				}
				pathVerticesEnds.put(nextVertex.getSequenceIndex(), rawConsensus.length());
			}
			lastVertex = nextVertex;
		}
		path.setConsensus(rawConsensus.toString());
		path.setPathVerticesConsensusEnds(pathVerticesEnds);
		log.info("Consensus built for path "+pathIdx+". Length: "+rawConsensus.length()+" Memory (Mbp): "+usedMemory);
		
	}
	
	public List<ReadAlignment> alignPathReads(AssemblyGraph graph, AssemblyPath path, int numThreads) {
		int debugIdx = -1;
		int n = path.getPathLength();
		int pathIdx = path.getPathId();
		List<ReadAlignment> alignedReads = new ArrayList<ReadAlignment>();
		String consensus = path.getConsensus();
		Map<Integer,Integer> pathVerticesEnds = path.getPathVerticesConsensusEnds();
		
		List<AssemblyEdge> edges = path.getEdges();
		int totalReads = 0;
		long usedMemory = (runtime.totalMemory()-runtime.freeMemory())/1000000;
		
		
		ThreadPoolManager poolAlign = new ThreadPoolManager(numThreads, Math.max(numThreads, 1000));
		AssemblyVertex lastVertex = path.getVertexLeft();
		for(int j = 0; j < n; j++) {
			AssemblyEdge edge = edges.get(j);
			AssemblyVertex nextVertex = edge.getConnectingVertex(lastVertex);
			if ((j+1)%500==0) {
				usedMemory = (runtime.totalMemory()-runtime.freeMemory())/1000000;
				log.info("Path "+pathIdx+". Aligning. Processed path edges: "+(j+1)+" of "+n+" alignments: "+alignedReads.size()+" Memory (Mbp): "+usedMemory);
			}
			if(!edge.isSameSequenceEdge()) {
				lastVertex = nextVertex;
				continue;		
			}
			//Align to consensus next path read and its embedded sequences
			int readIndex = lastVertex.getSequenceIndex();
			QualifiedSequence read = lastVertex.getRead(); 
			CharSequence seq = read.getCharacters();
			boolean reverse = !lastVertex.isStart();
			if(reverse) seq = DNASequence.getReverseComplement(seq);
			//String seqStr = seq.toString();
			int endConsensusPathVertex = Math.min(consensus.length(), pathVerticesEnds.get(readIndex));
			int startConsensusPathVertex = Math.max(0, endConsensusPathVertex-seq.length()-100);
			Map<Integer, Long> kmersSubject = kmerCodesTable.computeSequenceCodesAsMap(consensus, startConsensusPathVertex, endConsensusPathVertex);
			totalReads++;
			//Synchronic call to calculate actual backbone read ends
			ReadAlignment alnRead = alignReadProcess(pathIdx, consensus, kmersSubject, readIndex, read.getName(), seq, reverse, startConsensusPathVertex,endConsensusPathVertex, alignedReads);
			
			if(pathIdx == debugIdx) System.out.println("Consensus length: "+consensus.length()+" Limits consensus: "+startConsensusPathVertex+" "+endConsensusPathVertex+" Next path read: "+read.getName()+" sequence: "+seq.length()+" alignment: "+alnRead);
			
			if(alnRead!=null) {
				startConsensusPathVertex = alnRead.getFirst()-alnRead.getSoftClipStart();
				endConsensusPathVertex = alnRead.getLast()+alnRead.getSoftClipEnd();
				List<AssemblyEmbedded> embeddedList = graph.getAllEmbedded(readIndex);
				//List<AssemblyEmbedded> embeddedList = graph.getEmbeddedByHostId(vertexPreviousEdge.getSequenceIndex());
				for(AssemblyEmbedded embedded:embeddedList) {
					QualifiedSequence embeddedRead = embedded.getRead(); 
					CharSequence embeddedSeq = embeddedRead.getCharacters();
					boolean reverseE = (reverse!=embedded.isReverse());
					if(reverseE) embeddedSeq = DNAMaskedSequence.getReverseComplement(embeddedSeq);
					int startConsensusEmbedded = startConsensusPathVertex+embedded.getHostStart();
					int endConsensusEmbedded = startConsensusPathVertex+embedded.getHostEnd();
					if(reverse) {
						startConsensusEmbedded = startConsensusPathVertex+(seq.length()-embedded.getHostEnd());
						endConsensusEmbedded = startConsensusPathVertex+(seq.length()-embedded.getHostStart());
					}
					//if(embedded.getSequenceId()==1940) System.out.println("Consensus length: "+rawConsensus.length()+" limits: "+startConsensus+" "+endConsensus+" reverseEmb: "+reverseE+" host: "+readIndex+" "+read.getName()+" Reverse host: "+reverse+" rel: "+embedded);
					totalReads++;
					try {
						final int s = startConsensusEmbedded;
						final int e = endConsensusEmbedded;
						CharSequence q = embeddedSeq;
						poolAlign.queueTask(()->alignReadProcess(pathIdx, consensus, kmersSubject, embedded.getSequenceId(), embeddedRead.getName(), q, reverseE, s, e, alignedReads));
					} catch (InterruptedException e) {
						//TODO: Better handling
						e.printStackTrace();
					}
					//if (totalReads%1000==0) log.info("Path "+pathIdx+". Aligning. Processed reads: "+totalReads+" alignments: "+alignedReads.size()+" unaligned: "+unalignedReadIds.size());
				}
			} else {
				List<AssemblyEmbedded> embeddedList = graph.getAllEmbedded(readIndex);
				if(pathIdx == debugIdx) System.err.println("WARN: Unaligned consensus backbone read "+nextVertex.getRead().getName()+" to final consensus. Unaligned embedded reads: "+embeddedList.size());
				if (buildUnalignedReadRecords) {
					for(AssemblyEmbedded embedded:embeddedList) {
						QualifiedSequence seqEmb = graph.getSequence(embedded.getSequenceId());
						ReadAlignment unalignedReadRecord = buildUnalignedReadRecord(embedded.getSequenceId(),seqEmb.getName(),seqEmb.getCharacters());
						synchronized (alignedReads) {
							alignedReads.add(unalignedReadRecord);
						}
					}
				}	
			}
			lastVertex = nextVertex;
		}
		try {
			poolAlign.terminatePool();
		} catch (InterruptedException e) {
			// TODO Better handling
			e.printStackTrace();
		}
		alignedReads.addAll(alignInternalPaths(graph, path,numThreads));
		Collections.sort(alignedReads, ReadAlignmentPositionComparator.getInstance() );
		usedMemory = (runtime.totalMemory()-runtime.freeMemory())/1000000;
		log.info("Processed path "+pathIdx+". Length: "+path.getPathLength()+" Total reads: "+totalReads+" alignments: "+alignedReads.size()+" Memory (Mbp): "+usedMemory);
		return alignedReads;
	}
	private Map<Integer, Long> selectKmers(Map<Integer, Long> kmersSubject, int startConsensus, int endConsensus) {
		Map<Integer, Long> answer = new LinkedHashMap<Integer, Long>();
		for(Map.Entry<Integer, Long> entry:kmersSubject.entrySet()) {
			int pos = entry.getKey();
			if(pos>=startConsensus && pos<=endConsensus) answer.put(entry.getKey(), entry.getValue());
		}
		return answer;
	}
	private ReadAlignment alignReadProcess(int pathIdx, String consensus, Map<Integer, Long> kmersSubject,
			int readId, String readName, CharSequence sequence, boolean reverse, int startConsensus, int endConsensus, List<ReadAlignment> alignedReads) {
		int debugReadIdx = -1;
		//MinimizersTableReadAlignmentAlgorithm aligner = factory.requestLongReadsAligner(MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		LongReadsUngappedSearchHitsClusterAligner aligner = new LongReadsUngappedSearchHitsClusterAligner(LongReadsUngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		Map<Integer, Long> selKmersSubject = selectKmers(kmersSubject,startConsensus,endConsensus);
		//if(readId == 61) System.out.println("Query: "+sequence+"\nsubject: "+consensus.substring(startConsensus, endConsensus));
		ReadAlignment aln = alignRead(aligner,pathIdx, consensus, sequence, selKmersSubject);
		if(aln!=null) {
			aln.setReadName(readName);
			aln.setReadNumber(readId);
			aln.setNegativeStrand(reverse);
			aln.setReadCharacters(sequence);
			//aln.setReadCharacters(null);
			if(readId == debugReadIdx) System.out.println("Read name: "+readName+"First alignment attempt: "+aln+" mismatches: "+aln.getNumMismatches());
			synchronized (alignedReads) {
				alignedReads.add(aln);
			}
			return aln;
		}
		if(!haploid) {
			int n = sequence.length()/2;
			CharSequence subqueryLeft = sequence.subSequence(0, n);
			CharSequence subqueryRight = sequence.subSequence(n,sequence.length());
			ReadAlignment alnLeft = alignRead(aligner, pathIdx, subqueryLeft, readName, kmersSubject);
			if(readId==debugReadIdx) System.out.println("Read name: "+readName+" Subsequence length: "+subqueryLeft.length()+" Left aln: "+alnLeft);
			ReadAlignment alnRight = alignRead(aligner, pathIdx, subqueryRight, readName, kmersSubject);
			if(readId==debugReadIdx) System.out.println("Read name: "+readName+" Subsequence length: "+subqueryRight.length()+" Right aln: "+alnRight);
			ReadAlignment selected = alnLeft;
			if(selected == null || (alnRight!=null && alnRight.length()>alnLeft.length())) {
				selected = alnLeft;
			}
			if(selected!=null && selected == alnLeft) {
				selected.setReadName(readName);
				selected.setReadNumber(readId);
				selected.setNegativeStrand(reverse);
				selected.setReadCharacters(sequence);
				//TODO: Do this better
				selected.setCigarString(selected.getCigarString()+(sequence.length()-n)+"S");
				synchronized (alignedReads) {
					alignedReads.add(selected);
				}
				return selected;
			} else if (selected!=null) {
				selected.setReadName(readName);
				selected.setReadNumber(readId);
				selected.setNegativeStrand(reverse);
				selected.setReadCharacters(sequence);
				//TODO: Do this better
				selected.setCigarString(""+n+"S"+selected.getCigarString());
				synchronized (alignedReads) {
					alignedReads.add(selected);
				}
				return selected;
			}
		}
		ReadAlignment aln2 = alignRead(aligner,pathIdx, consensus, sequence, kmersSubject);
		if(aln2!=null) {
			System.err.println("WARN: Alignment found for previously unaligned read "+readId+" "+readName+" reverse: "+reverse+". Given limits: "+startConsensus+" "+endConsensus+" aln: "+aln2);
			aln2.setReadName(readName);
			aln2.setReadNumber(readId);
			aln2.setNegativeStrand(reverse);
			aln2.setReadCharacters(sequence);
			synchronized (alignedReads) {
				alignedReads.add(aln2);
			}
			return aln2;
		}
		if(buildUnalignedReadRecords && aln==null) {
			ReadAlignment unalignedReadRecord = buildUnalignedReadRecord(readId,readName,sequence);
			synchronized (alignedReads) {
				alignedReads.add(unalignedReadRecord);
			}
		}
		return aln;
	}
	private ReadAlignment buildUnalignedReadRecord(int readId, String readName, CharSequence sequence) {
		ReadAlignment aln = new ReadAlignment(0,0,0,sequence.length(),ReadAlignment.FLAG_READ_UNMAPPED);
		aln.setReadCharacters(sequence);
		return aln;
	}
	private ReadAlignment alignRead(LongReadsUngappedSearchHitsClusterAligner aligner, int subjectIdx, CharSequence subject, CharSequence read, Map<Integer, Long> codesSubject) {
		String readStr = read.toString();
		
		Map<Integer, Long> codesQuery = kmerCodesTable.computeSequenceCodesAsMap(readStr, 0, read.length());
		//if(read.length()==14871) System.out.println("Number of codes query: "+codesQuery.size());
		UngappedSearchHitsCluster bestCluster = PairwiseAlignerDynamicKmers.findBestKmersCluster(subject.length(), codesSubject, read.length(), codesQuery, KMER_LENGTH_LOCAL_ALN);
		if(bestCluster==null) return null;
		ReadAlignment aln;
		synchronized (aligner) {
			aln = aligner.buildAlignment(readStr, subject, bestCluster);
		}
		//if(read.length()==14871) System.out.println("Best cluster kmers: "+bestCluster.getNumDifferentKmers()+" alignment "+aln);
		if(!evaluateAlignment(aln)) aln = null;
		return aln;
	}
	private boolean evaluateAlignment(ReadAlignment aln) {
		if(aln==null) return false;
		Map<Integer,GenomicVariant> indelCalls = aln.getIndelCalls();
		if(!haploid) return true;
		if(indelCalls==null) return true;
		for(GenomicVariant indelCall:indelCalls.values()) {
			if(indelCall.length()>100) return false;
		}
		return true;
	}
	public List<CalledGenomicVariant> callIndels (String consensus, List<ReadAlignment> alignments, int normalPloidy) {
		String sequenceName = "";
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
			//TODO: check if it is worth to return homozygous reference calls
			if(!call.isUndecided() && !call.isHomozygousReference()) answer.add(call);
		}
		return answer;
	}
	private List<GenomicRegion> calculateActiveSegments(String sequenceName, List<ReadAlignment> alignments) {
		//Extract indel calls adding one bp on the sides for insertions
		List<GenomicRegion> rawRegions = new ArrayList<GenomicRegion>();
		for(ReadAlignment aln:alignments) {
			Map<Integer,GenomicVariant> indelCalls = aln.getIndelCalls();
			if(indelCalls==null) continue;
			for(GenomicVariant indelCall:indelCalls.values()) {
				if(indelCall.length()>10) {
					//if(indelCall.length()>100) System.out.println("WARN: Long indel from alignment: "+aln + "coordinates: "+indelCall.getFirst()+"-"+indelCall.getLast()+" Ignoring.");
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
				//if(nextRegion.length()>20) System.out.println("Adding long region "+nextRegion.getSequenceName()+":"+nextRegion.getFirst()+"-"+nextRegion.getLast()+" "+nextRegion.length()+" support: "+countSupport);
				if(countSupport>=5) mergedRegions.add(nextRegion);
				nextRegion = (GenomicRegionImpl)rawRegion;
				countSupport=1;
			}
		}
		//if(nextRegion.length()>20) System.out.println("Adding long region "+nextRegion.getSequenceName()+":"+nextRegion.getFirst()+"-"+nextRegion.getLast()+" "+nextRegion.length()+" support: "+countSupport);
		if(countSupport>=5) mergedRegions.add(nextRegion);
		return mergedRegions;
	}

	private String calculateLocalConsensus(int first, int last, List<ReadAlignment> alignments, int firstIdxAln, String consensusAllele) {
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

	private CalledGenomicVariant buildCall(String sequenceName, int first, String currentConsensus, String localConsensus, String altConsensus) {
		List<String> alleles = new ArrayList<String>(2);
		alleles.add(currentConsensus);
		boolean hetero = altConsensus!=null && !altConsensus.equals(localConsensus);
		if(localConsensus!=null && !localConsensus.equals(currentConsensus)) alleles.add(localConsensus);
		if(hetero && !altConsensus.equals(currentConsensus)) {
			alleles.add(altConsensus);
		}
		GenomicVariantImpl variant = new GenomicVariantImpl(sequenceName, first, alleles);
		CalledGenomicVariantImpl call;
		if (hetero) call = new CalledGenomicVariantImpl(variant, CalledGenomicVariant.GENOTYPE_HETERO);
		else if(alleles.size()==1) call = new CalledGenomicVariantImpl(variant, CalledGenomicVariant.GENOTYPE_HOMOREF);
		else call = new CalledGenomicVariantImpl(variant, CalledGenomicVariant.GENOTYPE_HOMOALT);
		return call;
	}
	private String calculateVariantSegment(ReadAlignment alignment, GenomicVariant indelReadCall, CalledGenomicVariant calledVariant, int normalPloidy) {
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
	private List<ReadAlignment> alignInternalPaths(AssemblyGraph graph, AssemblyPath path, int numThreads) {
		int debugIdx = -1;
		int pathIdx = path.getPathId();
		List<AssemblyPath> internalPaths = path.getAlternativeSmallPaths();
		if(internalPaths.size()==0 ) return new ArrayList<>();
		log.info("Aligning internal paths for path: "+path.getPathId()+" number of paths: "+internalPaths.size());
		ReadsAligner aligner = new ReadsAligner();
		ReferenceGenome genome = new ReferenceGenome(new QualifiedSequence("", path.getConsensus()));
		aligner.setGenome(genome);
		//log.info("Aligning internal paths for path: "+path.getPathId()+" loaded kmer codes table");
		List<ReadAlignment> alignedReads = new ArrayList<ReadAlignment>();
		int totalReads = 0;
		for(AssemblyPath internalPath: internalPaths) {
			List<AssemblyEdge> edges = internalPath.getEdges();
			int n = edges.size();
			for(int j = 0; j < n; j++) {
				AssemblyEdge edge = edges.get(j);
				if(!edge.isSameSequenceEdge()) continue;
				totalReads++;
				QualifiedSequence read = edge.getVertex1().getRead();
				List<ReadAlignment> alns = aligner.alignRead(read);
				if(pathIdx == debugIdx) System.err.println("Aligning internal paths for path: "+path.getPathId()+" Next read: "+read.getName()+" Alignments: "+alns);
				if(alns.size()==0) continue;
				ReadAlignment aln = alns.get(0);
				aln.setReadNumber(edge.getVertex1().getSequenceIndex());
				alignedReads.add(aln);
				List<AssemblyEmbedded> embeddedList = graph.getAllEmbedded(edge.getVertex1().getSequenceIndex());
				for(AssemblyEmbedded embedded:embeddedList) {
					QualifiedSequence embeddedRead = embedded.getRead();
					totalReads++;
					List<ReadAlignment> alnsE = aligner.alignRead(embeddedRead);
					if(alnsE.size()==0) continue;
					aln = alnsE.get(0);
					aln.setReadNumber(embedded.getSequenceId());
					alignedReads.add(aln);
				}
			}
		}
		log.info("Aligning internal paths for path: "+path.getPathId()+" Total reads internal paths: "+totalReads+" aligned: "+alignedReads.size());
		return alignedReads;
	}
}
