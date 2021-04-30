package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.sequences.AbstractLimitedSequence;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.UngappedSearchHit;

public class AssemblyPathReadsAligner {
	private Logger log = Logger.getLogger(AssemblyPathReadsAligner.class.getName());
	
	public static final int KMER_LENGTH_LOCAL_ALN = KmersExtractor.DEF_KMER_LENGTH;
	
	//Parameters
	private boolean onlyGenerateConsensus = false;
	private boolean alignEmbedded = false;
	
	//Output
	private StringBuilder consensus;
	private List<ReadAlignment> alignedReads;
	private Set<Integer> unalignedReadIds;
	
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
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
	public void alignPathReads(AssemblyGraph graph, AssemblyPath path) {
		int debugIdx = -1;
		int n = path.getPathLength();
		int pathIdx = path.getPathId();
		log.info("Aligning reads for path "+pathIdx+" "+path.getSequenceName()+" with length: "+n);
		List<AssemblyEdge> edges = path.getEdges();
		StringBuilder rawConsensus = new StringBuilder();
		AssemblyVertex lastVertex = path.getVertexLeft();
		MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm();
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
				Map<Long, Integer> uniqueKmersSubject = KmersExtractor.extractLocallyUniqueKmerCodes(rawConsensus, KMER_LENGTH_LOCAL_ALN, startSuffixConsensus,rawConsensus.length());
				String prefixQuery = nextPathSequence.subSequence(0, Math.min(nextPathSequence.length(), edge.getOverlap()+30)).toString();
				ReadAlignment alnRead = alignRead(aligner, pathIdx, rawConsensus, prefixQuery, uniqueKmersSubject);
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
				Map<Long, Integer> uniqueKmersSubject = KmersExtractor.extractLocallyUniqueKmerCodes(rawConsensus, KMER_LENGTH_LOCAL_ALN, Math.max(0, rawConsensus.length()-seq.length()),rawConsensus.length());
				totalReads++;
				ReadAlignment alnRead = alignRead(aligner, pathIdx, rawConsensus, seq, uniqueKmersSubject);
				if (alnRead!=null) {
					alnRead.setReadName(read.getName());
					alnRead.setReadNumber(readIndex);
					alnRead.setNegativeStrand(reverse);
					alignedReads.add(alnRead);
					if(alnRead.getSoftClipEnd()>10) {
						//log.warning("Weird alignment of consensus backbone read. Alignment: "+alnRead+" soft clip: "+alnRead.getSoftClipEnd()+" consensus end: "+rawConsensus.substring(rawConsensus.length()-alnRead.getSoftClipEnd()-10)+" soft clipped sequence: "+alnRead.getReadCharacters().subSequence(alnRead.getReadLength()-alnRead.getSoftClipEnd()-10, alnRead.getReadLength()) );
					}
				}
				else {
					if(pathIdx == debugIdx) System.err.println("Backbone read: "+read.getName()+" could not be aligned to extended consensus");
					unalignedReadIds.add(readIndex);
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
						totalReads++;
						ReadAlignment alnEmbedded = alignRead(aligner,pathIdx, rawConsensus, embeddedSeq, uniqueKmersSubject);
						if(alnEmbedded!=null) {
							alnEmbedded.setReadName(embeddedRead.getName());
							alnEmbedded.setReadNumber(embedded.getSequenceId());
							alnEmbedded.setNegativeStrand(reverseE);
							alignedReads.add(alnEmbedded);
						} else {
							if(pathIdx == debugIdx) System.err.println("Embedded read: "+read.getName()+" could not be aligned to extended consensus");
							unalignedReadIds.add(embedded.getSequenceId());
						}
						if (totalReads%1000==0) log.info("Path "+pathIdx+". Aligning. Processed reads: "+totalReads+" alignments: "+alignedReads.size()+" unaligned: "+unalignedReadIds.size());
					}
				}
			}
			lastVertex = nextVertex;
		}
		consensus = rawConsensus;
		log.info("Processed path "+pathIdx+". Length: "+path.getPathLength()+" Total reads: "+totalReads+" alignments: "+alignedReads.size()+" unaligned: "+unalignedReadIds.size());
	}
	public ReadAlignment alignRead(MinimizersTableReadAlignmentAlgorithm aligner, int subjectIdx, CharSequence subject, CharSequence read, int start, int end) {
		Map<Long, Integer> uniqueCodesSubject = KmersExtractor.extractLocallyUniqueKmerCodes(subject, KMER_LENGTH_LOCAL_ALN, start,end);
		//System.out.println("Number of unique k-mers subject: "+uniqueKmersSubject.size());
		return alignRead(aligner, subjectIdx, subject, read, uniqueCodesSubject);
	}
	public ReadAlignment alignRead(MinimizersTableReadAlignmentAlgorithm aligner, int subjectIdx, CharSequence subject, CharSequence read, Map<Long, Integer> uniqueCodesSubject) {
		Map<Integer, Long> codesQuery = KmersExtractor.extractDNAKmerCodes(read, KMER_LENGTH_LOCAL_ALN, 0, read.length());
		//System.out.println("Number of unique k-mers read: "+uniqueKmersRead.size());
		List<UngappedSearchHit> initialKmerHits = alignKmerCodes(-1,subject.length(), uniqueCodesSubject, codesQuery);
		if(initialKmerHits.size()==0) return null;
		List<KmerHitsCluster> clusters = KmerHitsCluster.clusterRegionKmerAlns(read.length(), subject.length(), initialKmerHits, 0);
		//printClusters(clusters);
		if(clusters.size()>1) {
			Collections.sort(clusters, (o1,o2)->o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
			KmerHitsCluster c1 = clusters.get(0);
			KmerHitsCluster c2 = clusters.get(1);
			int overlap = GenomicRegionSpanComparator.getInstance().getSpanLength(c1.getSubjectPredictedStart(), c1.getSubjectPredictedEnd(), c2.getSubjectPredictedStart(), c2.getSubjectPredictedEnd());
			int c1Length = c1.getSubjectPredictedEnd()-c1.getSubjectPredictedStart();
			int c2Length = c2.getSubjectPredictedEnd()-c2.getSubjectPredictedStart();
			if((overlap <0.9*c1Length || overlap < 0.9*c2Length) && c1.getNumDifferentKmers()<0.9*initialKmerHits.size()) {
				return null;
			}	
		} else if (clusters.size()==0) return null;
		KmerHitsCluster bestCluster = clusters.get(0);
		//System.out.println("Number of clusters: "+clusters.size()+" best cluster kmers: "+bestCluster.getNumDifferentKmers()+" first "+bestCluster.getFirst()+" last "+bestCluster.getLast());
		return aligner.buildCompleteAlignment(subjectIdx, subject, read, bestCluster);
	}
	private List<UngappedSearchHit> alignKmerCodes(int subjectIdx, int subjectLength, Map<Long, Integer> uniqueCodesSubject, Map<Integer, Long> codesQuery) {
		List<UngappedSearchHit> initialKmerHits = new ArrayList<UngappedSearchHit>();
		for(int i:codesQuery.keySet()) {
			Long codeRead = codesQuery.get(i);
			Integer subjectPos = uniqueCodesSubject.get(codeRead);
			if(subjectPos==null) continue;
			CharSequence kmerRead = new String(AbstractLimitedSequence.getSequence(codeRead, KMER_LENGTH_LOCAL_ALN, DNASequence.EMPTY_DNA_SEQUENCE));
			UngappedSearchHit hit = new UngappedSearchHit(kmerRead, subjectIdx , subjectPos);
			hit.setQueryIdx(i);
			initialKmerHits.add(hit);
		}
		return initialKmerHits;
	}
}
