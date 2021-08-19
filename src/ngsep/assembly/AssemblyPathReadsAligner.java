package ngsep.assembly;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.PairwiseAlignerDynamicKmers;
import ngsep.alignments.UngappedSearchHitsCluster;
import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.GenomicVariant;

public class AssemblyPathReadsAligner {
	private Logger log = Logger.getLogger(AssemblyPathReadsAligner.class.getName());
	
	public static final int KMER_LENGTH_LOCAL_ALN = 31;
	
	//Parameters
	private boolean onlyGenerateConsensus = false;
	private boolean alignEmbedded = false;
	private int numThreads = 1;
	
	//Output
	private CharSequence consensus;
	private List<ReadAlignment> alignedReads;
	
	//private LongReadsAlignerFactory factory = new LongReadsAlignerFactory();
	
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
		//factory.setLog(log);
	}
	
	public boolean isOnlyGenerateConsensus() {
		return onlyGenerateConsensus;
	}
	public void setOnlyGenerateConsensus(boolean onlyGenerateConsensus) {
		this.onlyGenerateConsensus = onlyGenerateConsensus;
	}
	public boolean isAlignEmbedded() {
		return alignEmbedded;
	}
	public void setAlignEmbedded(boolean alignEmbedded) {
		this.alignEmbedded = alignEmbedded;
	}
	public CharSequence getConsensus() {
		return consensus;
	}
	public List<ReadAlignment> getAlignedReads() {
		return alignedReads;
	}
	
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
		//factory.setNumThreads(numThreads);
	}
	public void alignPathReads(AssemblyGraph graph, AssemblyPath path) {
		int debugIdx = -1;
		int n = path.getPathLength();
		int pathIdx = path.getPathId();
		consensus = null;
		alignedReads = new ArrayList<ReadAlignment>();
		Runtime runtime = Runtime.getRuntime();
		long usedMemory = (runtime.totalMemory()-runtime.freeMemory())/1000000;
		log.info("Building consensus for path "+pathIdx+" "+path.getSequenceName()+" with length: "+n+" Memory (Mbp): "+usedMemory);
		List<AssemblyEdge> edges = path.getEdges();
		Map<Integer,Integer> pathVerticesEnds = new HashMap<Integer, Integer>();
		StringBuilder rawConsensus = new StringBuilder();
		AssemblyVertex lastVertex = path.getVertexLeft();
		//Build consensus first
		//MinimizersTableReadAlignmentAlgorithm aligner = factory.requestLongReadsAligner(MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm(MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		int totalReads = 0;
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
				if(pathIdx == debugIdx) System.err.println("Added sequence: "+lastVertex.getRead().getName()+" reverse: "+reverse);
			} else if(!edge.isSameSequenceEdge()) {
				// Augment consensus with the next path read
				CharSequence nextPathSequence = nextVertex.getRead().getCharacters();
				boolean reverse = !nextVertex.isStart();
				if(reverse) nextPathSequence = DNAMaskedSequence.getReverseComplement(nextPathSequence);
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) printAllOverlappingSeqs(graph,path,j,vertexPreviousEdge);
				if(pathIdx == debugIdx && j<10) System.err.println("Aligning next path read "+nextVertex.getRead().getName()+". Reverse "+reverse+ " edge: "+edge);
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
					if(pathIdx == debugIdx && j<10) System.err.println("Sequence length: "+nextPathSequence.length()+" subject length: "+rawConsensus.length()+" Soft clip end: "+alnRead.getSoftClipEnd()+" pos aln: "+posAlnRead+" pos subject: "+lastPosSubject+" tail: "+tailSubject+" aln: "+alnRead);
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
					if(pathIdx == debugIdx && j<10) System.err.println("Enlarging consensus. Start new sequence: "+startSuffixQuery+" length of segment to append: "+remainingSegment.length());
					if(startRemove>0) rawConsensus.delete(startRemove, rawConsensus.length());
					rawConsensus.append(remainingSegment.toUpperCase());
				}
				pathVerticesEnds.put(nextVertex.getSequenceIndex(), rawConsensus.length());
			}
			lastVertex = nextVertex;
		}
		consensus = new DNASequence(rawConsensus);
		usedMemory = (runtime.totalMemory()-runtime.freeMemory())/1000000;
		log.info("Consensus built for path "+pathIdx+". Length: "+rawConsensus.length()+" Memory (Mbp): "+usedMemory+" Aligning reads.");
		if(onlyGenerateConsensus) return;
		
		ThreadPoolManager poolAlign = new ThreadPoolManager(numThreads, Math.max(numThreads, 1000));
		lastVertex = path.getVertexLeft();
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
			int endConsensusPathVertex = Math.min(rawConsensus.length(), pathVerticesEnds.get(readIndex));
			int startConsensusPathVertex = Math.max(0, endConsensusPathVertex-seq.length()-100);
			Map<Integer, Long> kmersSubject = KmersExtractor.extractDNAKmerCodes(rawConsensus, KMER_LENGTH_LOCAL_ALN, startConsensusPathVertex,endConsensusPathVertex);
			totalReads++;
			//Synchronic call to calculate actual backbone read ends
			ReadAlignment alnRead = alignReadProcess(pathIdx, rawConsensus, kmersSubject, readIndex, read.getName(), seq, reverse, startConsensusPathVertex,endConsensusPathVertex);
			
			if(pathIdx == debugIdx) System.out.println("Consensus length: "+rawConsensus.length()+" Limits consensus: "+startConsensusPathVertex+" "+endConsensusPathVertex+" Next path read: "+read.getName()+" sequence: "+seq.length()+" alignment: "+alnRead);
			
			if(!alignEmbedded) {
				lastVertex = nextVertex;
				continue;
			}
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
						poolAlign.queueTask(()->alignReadProcess(pathIdx, rawConsensus, kmersSubject, embedded.getSequenceId(), embeddedRead.getName(), q, reverseE, s, e));
					} catch (InterruptedException e) {
						//TODO: Better handling
						e.printStackTrace();
					}
					//if (totalReads%1000==0) log.info("Path "+pathIdx+". Aligning. Processed reads: "+totalReads+" alignments: "+alignedReads.size()+" unaligned: "+unalignedReadIds.size());
				}
			} else {
				List<AssemblyEmbedded> embeddedList = graph.getAllEmbedded(readIndex);
				if(pathIdx == debugIdx && j<10) System.err.println("WARN: Unaligned consensus backbone read "+nextVertex.getRead().getName()+" to final consensus. Unaligned embedded reads: "+embeddedList.size());
			}
			lastVertex = nextVertex;
		}
		try {
			poolAlign.terminatePool();
		} catch (InterruptedException e) {
			// TODO Better handling
			e.printStackTrace();
		}
		usedMemory = (runtime.totalMemory()-runtime.freeMemory())/1000000;
		log.info("Processed path "+pathIdx+". Length: "+path.getPathLength()+" Total reads: "+totalReads+" alignments: "+alignedReads.size()+" Memory (Mbp): "+usedMemory);
	}
	private Map<Integer, Long> selectKmers(Map<Integer, Long> kmersSubject, int startConsensus, int endConsensus) {
		Map<Integer, Long> answer = new LinkedHashMap<Integer, Long>();
		for(Map.Entry<Integer, Long> entry:kmersSubject.entrySet()) {
			int pos = entry.getKey();
			if(pos>=startConsensus && pos<=endConsensus) answer.put(entry.getKey(), entry.getValue());
		}
		return answer;
	}
	private ReadAlignment alignReadProcess(int pathIdx, StringBuilder subject, Map<Integer, Long> kmersSubject,
			int readId, String readName, CharSequence sequence, boolean reverse, int startConsensus, int endConsensus) {
		//MinimizersTableReadAlignmentAlgorithm aligner = factory.requestLongReadsAligner(MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm(MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		Map<Integer, Long> selKmersSubject = selectKmers(kmersSubject,startConsensus,endConsensus);
		ReadAlignment aln = alignRead(aligner,pathIdx, subject, sequence, selKmersSubject);
		if(aln!=null) {
			aln.setReadName(readName);
			aln.setReadNumber(readId);
			aln.setNegativeStrand(reverse);
			aln.setReadCharacters(sequence);
			//aln.setReadCharacters(null);
			synchronized (alignedReads) {
				alignedReads.add(aln);
			}
		} else {
			ReadAlignment aln2 = alignRead(aligner,pathIdx, subject, sequence, kmersSubject);
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
		}
		return aln;
	}
	private ReadAlignment alignRead(MinimizersTableReadAlignmentAlgorithm aligner, int subjectIdx, CharSequence subject, CharSequence read, Map<Integer, Long> codesSubject) {
		String readStr = read.toString();
		Map<Integer, Long> codesQuery = KmersExtractor.extractDNAKmerCodes(readStr, KMER_LENGTH_LOCAL_ALN, 0, read.length());
		//System.out.println("Number of unique k-mers read: "+uniqueKmersRead.size());
		UngappedSearchHitsCluster bestCluster = PairwiseAlignerDynamicKmers.findBestKmersCluster(subject.length(), codesSubject, read.length(), codesQuery, KMER_LENGTH_LOCAL_ALN);
		if(bestCluster==null) return null;
		ReadAlignment aln;
		synchronized (aligner) {
			aln = aligner.buildCompleteAlignment(subjectIdx, subject, readStr, bestCluster);
		}
		if(!evaluateAlignment(aln)) return null;
		//System.out.println("Number of clusters: "+clusters.size()+" best cluster kmers: "+bestCluster.getNumDifferentKmers()+" first "+bestCluster.getFirst()+" last "+bestCluster.getLast());
		return aln;
	}
	private boolean evaluateAlignment(ReadAlignment aln) {
		if(aln==null) return false;
		Map<Integer,GenomicVariant> indelCalls = aln.getIndelCalls();
		if(indelCalls==null) return true;
		for(GenomicVariant indelCall:indelCalls.values()) {
			if(indelCall.length()>100) return false;
		}
		return true;
	}
}
