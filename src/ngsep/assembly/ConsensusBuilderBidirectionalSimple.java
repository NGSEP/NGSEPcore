package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

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

public class ConsensusBuilderBidirectionalSimple implements ConsensusBuilder {
	
	private static int kmerLength = KmersExtractor.DEF_KMER_LENGTH;
	
	@Override
	public List<QualifiedSequence> makeConsensus(AssemblyGraph graph) 
	{
		//List of final contigs
		List<QualifiedSequence> consensusList = new ArrayList<QualifiedSequence>();
		List<List<AssemblyEdge>> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			List<AssemblyEdge> path = paths.get(i);
			String sequenceName = "Contig_"+(i+1);
			CharSequence consensusPath = makeConsensus (graph, path, i, sequenceName);
			consensusList.add(new QualifiedSequence(sequenceName,consensusPath));
		}
		
		return consensusList;
	}
	
	private CharSequence makeConsensus(AssemblyGraph graph, List<AssemblyEdge> path, int sequenceIdx, String sequenceName) 
	{
		StringBuilder consensus = new StringBuilder();
		AssemblyVertex lastVertex = null;
		MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm();
		String pathS = "";
		if(path.size()==1) {
			consensus.append(path.get(0).getVertex1().getRead());
			return consensus;
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
			} else if (lastVertex == edge.getVertex1()) {
				vertexPreviousEdge = edge.getVertex1();
				vertexNextEdge = edge.getVertex2();
			} else if (lastVertex == edge.getVertex2()) {
				vertexPreviousEdge = edge.getVertex2();
				vertexNextEdge = edge.getVertex1();
			} else {
				throw new RuntimeException("Inconsistency found in path");
			}
			if(j == 0) {
				pathS = pathS.concat(vertexPreviousEdge.getUniqueNumber() + ",");
				CharSequence seq = vertexPreviousEdge.getRead().getCharacters();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
				consensus.append(seq);
			} else if(vertexPreviousEdge.getRead()!=vertexNextEdge.getRead()) {
				// Augment consensus with the next path read
				CharSequence nextPathSequence = vertexNextEdge.getRead().getCharacters();
				boolean reverse = !vertexNextEdge.isStart();
				if(reverse) nextPathSequence = DNAMaskedSequence.getReverseComplement(nextPathSequence);
				//if (rawConsensus.length()>490000 && rawConsensus.length()<530000) printAllOverlappingSeqs(graph,path,j,vertexPreviousEdge);
				
				ReadAlignment alnRead = alignRead(aligner, sequenceIdx, consensus, nextPathSequence.toString(), Math.max(0, consensus.length()-nextPathSequence.length()),consensus.length(), 0.5);
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
						startSuffix = posAlnRead + (consensus.length()-lastPosSubject+1);
					} else {
						startSuffix = edge.getOverlap();
					}
					//System.out.println("Calculated overlap from alignment: "+startSuffix+" alignment: "+alnRead+" edge: "+edge );
				} else {
					startSuffix = edge.getOverlap();
				}
				if(startSuffix<nextPathSequence.length()) {
					pathS = pathS.concat(vertexNextEdge.getUniqueNumber() + ",");
					String remainingSegment = nextPathSequence.subSequence(startSuffix, nextPathSequence.length()).toString();
					//if (consensus.length()>490000 && consensus.length()<510000) System.out.println("Consensus length: "+consensus.length()+" Vertex: "+vertexNextEdge.getUniqueNumber()+" read length: "+seq.length()+" overlap: "+edge.getOverlap()+" remaining: "+remainingSegment.length());
					consensus.append(remainingSegment.toUpperCase());
				}
			}
			lastVertex = vertexNextEdge;
		}
		System.out.println(pathS);
		return consensus;
	}

	public static ReadAlignment alignRead(MinimizersTableReadAlignmentAlgorithm aligner, int subjectIdx, CharSequence subject, CharSequence read, int start, int end, double minQueryCoverage) {
		Map<Long, Integer> uniqueCodesSubject = KmersExtractor.extractLocallyUniqueKmerCodes(subject,kmerLength, start,end);
		//System.out.println("Number of unique k-mers subject: "+uniqueKmersSubject.size());
		return alignRead(aligner, subjectIdx, subject, read, uniqueCodesSubject, minQueryCoverage);
	}
	public static ReadAlignment alignRead(MinimizersTableReadAlignmentAlgorithm aligner, int subjectIdx, CharSequence subject, CharSequence read, Map<Long, Integer> uniqueCodesSubject, double minQueryCoverage) {
		Map<Long, Integer> uniqueCodesRead = KmersExtractor.extractLocallyUniqueKmerCodes(read,kmerLength,0,read.length());
		//System.out.println("Number of unique k-mers read: "+uniqueKmersRead.size());
		List<UngappedSearchHit> initialKmerHits = alignUniqueKmerCodes(-1,subject.length(),uniqueCodesSubject, uniqueCodesRead);
		if(initialKmerHits.size()==0) return null;
		List<KmerHitsCluster> clusters = KmerHitsCluster.clusterRegionKmerAlns(read.length(), subject.length(), initialKmerHits, minQueryCoverage);
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
	private static List<UngappedSearchHit> alignUniqueKmerCodes(int subjectIdx, int subjectLength, Map<Long, Integer> uniqueCodesSubject, Map<Long, Integer> uniqueCodesQuery) {
		List<UngappedSearchHit> initialKmerHits = new ArrayList<UngappedSearchHit>();
		for(Long codeRead:uniqueCodesQuery.keySet()) {
			Integer subjectPos = uniqueCodesSubject.get(codeRead);
			if(subjectPos==null) continue;
			CharSequence kmerRead = new String(AbstractLimitedSequence.getSequence(codeRead, 15, DNASequence.EMPTY_DNA_SEQUENCE));
			UngappedSearchHit hit = new UngappedSearchHit(kmerRead, subjectIdx , subjectPos);
			hit.setQueryIdx(uniqueCodesQuery.get(codeRead));
			initialKmerHits.add(hit);
		}
		return initialKmerHits;
	}
}
