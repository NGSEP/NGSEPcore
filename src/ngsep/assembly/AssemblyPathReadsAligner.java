package ngsep.assembly;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.alignments.PairwiseAlignerDynamicKmers;
import ngsep.alignments.UngappedSearchHitsCluster;
import ngsep.alignments.LongReadsAlignerFactory;
import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;

public class AssemblyPathReadsAligner {
	private Logger log = Logger.getLogger(AssemblyPathReadsAligner.class.getName());
	
	public static final int KMER_LENGTH_LOCAL_ALN = 31;
	
	//Parameters
	private boolean onlyGenerateConsensus = false;
	private boolean alignEmbedded = false;
	private int numThreads = 1;
	
	//Output
	private StringBuilder consensus;
	private List<ReadAlignment> alignedReads;
	private Set<Integer> unalignedReadIds;
	
	private LongReadsAlignerFactory factory = new LongReadsAlignerFactory();
	
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
		factory.setLog(log);
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
	public StringBuilder getConsensus() {
		return consensus;
	}
	public List<ReadAlignment> getAlignedReads() {
		return alignedReads;
	}
	
	public Set<Integer> getUnalignedReadIds() {
		return unalignedReadIds;
	}
	
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
		factory.setNumThreads(numThreads);
	}
	public void alignPathReads(AssemblyGraph graph, AssemblyPath path) {
		int debugIdx = -1;
		int n = path.getPathLength();
		int pathIdx = path.getPathId();
		log.info("Aligning reads for path "+pathIdx+" "+path.getSequenceName()+" with length: "+n);
		List<AssemblyEdge> edges = path.getEdges();
		StringBuilder rawConsensus = new StringBuilder();
		AssemblyVertex lastVertex = path.getVertexLeft();
		MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm(MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		ThreadPoolManager poolAlign = new ThreadPoolManager(numThreads, Math.max(numThreads, 100));
		alignedReads = new ArrayList<ReadAlignment>();
		unalignedReadIds = new HashSet<>();
		int totalReads = 0;
		for(int j = 0; j < n; j++) {
			AssemblyEdge edge = edges.get(j);
			AssemblyVertex nextVertex = edge.getConnectingVertex(lastVertex);
			if(nextVertex == null) throw new RuntimeException("Inconsistency found in path. Previouus vertex: "+lastVertex+" edge: "+edge);
			if(j == 0) {
				CharSequence seq = lastVertex.getRead().getCharacters();
				boolean reverse = !lastVertex.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
				rawConsensus.append(seq);
				if(pathIdx == debugIdx) System.err.println("Added sequence: "+lastVertex.getRead().getName()+" reverse: "+reverse);
			} else if(!edge.isSameSequenceEdge()) {
				// Augment consensus with the next path read
				CharSequence nextPathSequence = nextVertex.getRead().getCharacters();
				boolean reverse = !nextVertex.isStart();
				if(reverse) nextPathSequence = DNAMaskedSequence.getReverseComplement(nextPathSequence);
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) printAllOverlappingSeqs(graph,path,j,vertexPreviousEdge);
				if(pathIdx == debugIdx && j<10) System.err.println("Aligning next path read "+nextVertex.getRead().getName()+". Reverse "+reverse+ " edge: "+edge);
				int startSuffixConsensus = Math.max(0, rawConsensus.length()-edge.getOverlap()-30);
				Map<Integer, Long> kmersSubject = KmersExtractor.extractDNAKmerCodes(rawConsensus, KMER_LENGTH_LOCAL_ALN, startSuffixConsensus,rawConsensus.length());
				String prefixQuery = nextPathSequence.subSequence(0, Math.min(nextPathSequence.length(), edge.getOverlap()+30)).toString();
				ReadAlignment alnRead = alignRead(aligner, pathIdx, rawConsensus, prefixQuery, kmersSubject);
				int startRemove = -1;
				int startSuffixQuery;
				if(alnRead!=null) {
					alnRead.setReadName(nextVertex.getRead().getName());
					int posAlnRead = prefixQuery.length()-1-alnRead.getSoftClipEnd();
					int lastPosSubject = alnRead.getReferencePositionAlignedRead(posAlnRead);
					int tailSubject = rawConsensus.length()-lastPosSubject-1;
					if(pathIdx == debugIdx && j<10) System.err.println("Sequence length: "+nextPathSequence.length()+" subject length: "+rawConsensus.length()+" Soft clip end: "+alnRead.getSoftClipEnd()+" pos aln: "+posAlnRead+" pos subject: "+lastPosSubject+" aln: "+alnRead);
					//if(alnRead.getSoftClipEnd()>0 && lastPosSubject>=0 && tailSubject>50) System.err.println("Large subject tail not aligned. Tail length "+tailSubject+" new sequence suffix: "+(lastPosSubject+1)+" Next path alignment: "+alnRead+" Tail of subject: "+rawConsensus.substring(lastPosSubject+1)+" end read: "+nextPathSequence.subSequence(posAlnRead+1, nextPathSequence.length()));
					//Just in case cycle but if the read aligns this should not enter
					while(posAlnRead>0 && lastPosSubject<0) {
						if(pathIdx == debugIdx && j<10) System.err.println("Negative pos subject: "+lastPosSubject+" for read position: "+posAlnRead+" read length: "+nextPathSequence.length());
						posAlnRead--;
						lastPosSubject = alnRead.getReferencePositionAlignedRead(posAlnRead);
						tailSubject = rawConsensus.length()-lastPosSubject-1;
					}
					if(lastPosSubject>=0) {
						startSuffixQuery = posAlnRead + 1;
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
					if(pathIdx == debugIdx && j<10) System.err.println("Consensus backbone read "+nextVertex.getRead().getName()+" did not align to last consensus. Using overlap: "+edge.getOverlap());
					startSuffixQuery = edge.getOverlap();
				}
				if(pathIdx == debugIdx && j<10) System.err.println("Start suffix: "+startSuffixQuery+" next seq len: "+nextPathSequence.length()+" start remove previous: "+startRemove);
				if(startSuffixQuery<nextPathSequence.length()) {
					String remainingSegment = nextPathSequence.subSequence(startSuffixQuery, nextPathSequence.length()).toString();
					if(pathIdx == debugIdx && j<10) System.err.println("Enlarging consensus. Start new sequence: "+startSuffixQuery+" segment length: "+remainingSegment.length());
					if(startRemove>0) rawConsensus.delete(startRemove, rawConsensus.length());
					rawConsensus.append(remainingSegment.toUpperCase());
				}
			}
			if(!onlyGenerateConsensus && edge.isSameSequenceEdge()) {
				//Align to consensus next path read and its embedded sequences
				int readIndex = lastVertex.getSequenceIndex();
				QualifiedSequence read = lastVertex.getRead(); 
				CharSequence seq = read.getCharacters();
				boolean reverse = !lastVertex.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
				String seqStr = seq.toString();
				int startConsensus = Math.max(0, rawConsensus.length()-seq.length());
				Map<Integer, Long> kmersSubject = KmersExtractor.extractDNAKmerCodes(rawConsensus, KMER_LENGTH_LOCAL_ALN, startConsensus,rawConsensus.length());
				totalReads++;
				try {
					final int s = startConsensus;
					poolAlign.queueTask(()->alignReadProcess(pathIdx, rawConsensus, kmersSubject, readIndex, read.getName(), seqStr, reverse, s,rawConsensus.length()));
				} catch (InterruptedException e) {
					//TODO: Better handling
					e.printStackTrace();
				}
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) System.out.println("Consensus length: "+rawConsensus.length()+" Vertex: "+vertexNextEdge.getUniqueNumber()+" sequence: "+read.length()+" alignment: "+alnRead);
				if (totalReads%1000==0) log.info("Path "+pathIdx+". Aligning. Processed reads: "+totalReads+" alignments: "+alignedReads.size()+" unaligned: "+unalignedReadIds.size());
				if(alignEmbedded) {
					List<AssemblyEmbedded> embeddedList = graph.getAllEmbedded(readIndex);
					//List<AssemblyEmbedded> embeddedList = graph.getEmbeddedByHostId(vertexPreviousEdge.getSequenceIndex());
					for(AssemblyEmbedded embedded:embeddedList) {
						QualifiedSequence embeddedRead = embedded.getRead(); 
						CharSequence embeddedSeq = embeddedRead.getCharacters();
						boolean reverseE = (reverse!=embedded.isReverse());
						if(reverseE) embeddedSeq = DNAMaskedSequence.getReverseComplement(embeddedSeq);
						String embeddedString = embeddedSeq.toString();
						int endConsensus = rawConsensus.length();
						if(reverse) endConsensus -= embedded.getHostStart();
						else endConsensus -= (seq.length()-embedded.getHostEnd());
						startConsensus = endConsensus-embeddedRead.getLength();
						//if(embedded.getSequenceId()==1940) System.out.println("Consensus length: "+rawConsensus.length()+" limits: "+startConsensus+" "+endConsensus+" reverseEmb: "+reverseE+" host: "+readIndex+" "+read.getName()+" Reverse host: "+reverse+" rel: "+embedded);
						totalReads++;
						try {
							final int s = startConsensus;
							final int e = endConsensus;
							poolAlign.queueTask(()->alignReadProcess(pathIdx, rawConsensus, kmersSubject, embedded.getSequenceId(), embeddedRead.getName(), embeddedString, reverseE, s, e));
						} catch (InterruptedException e) {
							//TODO: Better handling
							e.printStackTrace();
						}
						if (totalReads%1000==0) log.info("Path "+pathIdx+". Aligning. Processed reads: "+totalReads+" alignments: "+alignedReads.size()+" unaligned: "+unalignedReadIds.size());
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
		consensus = rawConsensus;
		log.info("Processed path "+pathIdx+". Length: "+path.getPathLength()+" Total reads: "+totalReads+" alignments: "+alignedReads.size()+" unaligned: "+unalignedReadIds.size());
	}
	private Map<Integer, Long> selectKmers(Map<Integer, Long> kmersSubject, int startConsensus, int endConsensus) {
		Map<Integer, Long> answer = new LinkedHashMap<Integer, Long>();
		for(Map.Entry<Integer, Long> entry:kmersSubject.entrySet()) {
			int pos = entry.getKey();
			if(pos>=startConsensus && pos<=endConsensus) answer.put(entry.getKey(), entry.getValue());
		}
		return answer;
	}
	private void alignReadProcess(int pathIdx, StringBuilder rawConsensus, Map<Integer, Long> kmersSubject,
			int readId, String readName, CharSequence embeddedSeq, boolean reverse, int startConsensus, int endConsensus) {
		MinimizersTableReadAlignmentAlgorithm aligner = factory.requestLongReadsAligner(MinimizersTableReadAlignmentAlgorithm.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		Map<Integer, Long> selKmersSubject = selectKmers(kmersSubject,startConsensus,endConsensus);
		ReadAlignment aln = alignRead(aligner,pathIdx, rawConsensus, embeddedSeq, selKmersSubject);
		if(aln!=null) {
			aln.setReadName(readName);
			aln.setReadNumber(readId);
			aln.setNegativeStrand(reverse);
			synchronized (alignedReads) {
				alignedReads.add(aln);
			}
		} else {
			ReadAlignment aln2 = alignRead(aligner,pathIdx, rawConsensus, embeddedSeq, kmersSubject);
			if(aln2!=null) System.err.println("WARN: Alignment found for previously unaligned read "+readId+" "+readName+" reverse: "+reverse+". Given limits: "+startConsensus+" "+endConsensus+" aln: "+aln2);
			synchronized (unalignedReadIds) {
				unalignedReadIds.add(readId);
			}
			
		}
	}
	public ReadAlignment alignRead(MinimizersTableReadAlignmentAlgorithm aligner, int subjectIdx, CharSequence subject, CharSequence read, Map<Integer, Long> codesSubject) {
		Map<Integer, Long> codesQuery = KmersExtractor.extractDNAKmerCodes(read.toString(), KMER_LENGTH_LOCAL_ALN, 0, read.length());
		//System.out.println("Number of unique k-mers read: "+uniqueKmersRead.size());
		UngappedSearchHitsCluster bestCluster = PairwiseAlignerDynamicKmers.findBestKmersCluster(subject.length(), codesSubject, read.length(), codesQuery, KMER_LENGTH_LOCAL_ALN);
		if(bestCluster==null) return null;
		ReadAlignment aln;
		synchronized (aligner) {
			aln = aligner.buildCompleteAlignment(subjectIdx, subject, read, bestCluster);
		}
		//System.out.println("Number of clusters: "+clusters.size()+" best cluster kmers: "+bestCluster.getNumDifferentKmers()+" first "+bestCluster.getFirst()+" last "+bestCluster.getLast());
		return aln;
	}
}
