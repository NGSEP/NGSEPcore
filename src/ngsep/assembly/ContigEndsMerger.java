package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.alignments.UngappedSearchHitsCluster;
import ngsep.alignments.UngappedSearchHitsClusterBuilder;
import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.ShortKmerCodesTable;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.io.FastaSequencesHandler;

public class ContigEndsMerger {

	private Logger log = Logger.getLogger(ContigEndsMerger.class.getName());
	private static final int END_LENGTH = 50000;
	private MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm();
	public static void main(String[] args) throws Exception {
		ContigEndsMerger instance = new ContigEndsMerger();
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> contigs = handler.loadSequences(args[0]);
		List<QualifiedSequence> answer = instance.mergeContigs(contigs);
		handler.saveSequences(answer, System.out, 100);
	}
	
	public List<QualifiedSequence> mergeContigs(List<QualifiedSequence> contigs) {
		
		ShortKmerCodesTable table = new ShortKmerCodesTable(25, 40);
		Map<Integer,String> contigEndsMap = new HashMap<Integer, String>();
		List<QualifiedSequence> contigsForGraph = new ArrayList<QualifiedSequence>(contigs.size());
		List<QualifiedSequence> smallContigs = new ArrayList<QualifiedSequence>(contigs.size());
		int i=0;
		for (QualifiedSequence contig:contigs) {
			CharSequence seq = contig.getCharacters();
			if(seq.length()<3*END_LENGTH) {
				smallContigs.add(contig);
				continue;
			}
			contigsForGraph.add(contig);
			String startSeq = seq.subSequence(0, END_LENGTH).toString();
			contigEndsMap.put(2*i, startSeq);
			table.addSequence(2*i, startSeq);
			String endSeq = seq.subSequence(seq.length()-END_LENGTH,seq.length()).toString();
			contigEndsMap.put(2*i+1, endSeq);
			table.addSequence(2*i+1, endSeq);
			i++;
		}
		if(contigsForGraph.size()<2) {
			log.info("Not enough large contigs to join");
			return contigs;
		}
		log.info("Built kmer codes table");
		AssemblyGraph graph = new AssemblyGraph(contigsForGraph);
		for(Map.Entry<Integer, String> entry:contigEndsMap.entrySet()) {
			int j = entry.getKey();
			String seq = entry.getValue();
			QualifiedSequence contig = contigsForGraph.get(j/2);
			log.info("Processing end: "+j+" of contig "+contig.getName()+" contig length: "+contig.getLength()+" query length: "+seq.length());
			Map<Integer,List<UngappedSearchHit>> hitsForward = table.match(j, seq);
			int maxCount = filterHits(j, hitsForward);
			for(int subjectId:hitsForward.keySet()) System.err.println("Next subject: "+subjectId+" contig name: "+contigsForGraph.get(subjectId/2).getName()+" count forward: "+getQueryIdxsCount(hitsForward.get(subjectId))+" maxcount: "+maxCount);
			String seqRev = DNAMaskedSequence.getReverseComplement(seq).toString();
			Map<Integer,List<UngappedSearchHit>> hitsReverse = table.match(j, seqRev);
			maxCount = Math.max(maxCount, filterHits(j, hitsReverse));
			for(int subjectId:hitsReverse.keySet()) System.err.println("Next subject: "+subjectId+" contig name: "+contigsForGraph.get(subjectId/2).getName()+" count reverse: "+getQueryIdxsCount(hitsReverse.get(subjectId))+" maxcount: "+maxCount);
			if(hitsForward.size()<20) buildEdges(graph,j,seq,contigEndsMap,hitsForward,false,maxCount/4);
			if(hitsReverse.size()<20) buildEdges(graph,j,seqRev,contigEndsMap,hitsReverse,true,maxCount/4);
			
			log.info("Processed end: "+j+" of contig "+contig.getName()+" contig length: "+contig.getLength()+" query length: "+seq.length()+" matches: "+hitsForward.size()+" "+hitsReverse.size());
		}
		LayoutBuilderGreedyMaxOverlap builder = new LayoutBuilderGreedyMaxOverlap();
		builder.findPaths(graph);
		List<AssemblyPath> paths = graph.getPaths();
		Set<Integer> sequenceIdsInPaths = new HashSet<Integer>();
		System.err.println("Built paths. Number: "+paths.size());
		for(AssemblyPath path:paths) {
			System.err.println("Path: "+path.getPathLength());
			List<AssemblyEdge> edges = path.getEdges();
			for(AssemblyEdge edge:edges) {
				if(!edge.isSameSequenceEdge()) System.err.println("Next edge: "+edge);
				else sequenceIdsInPaths.add(edge.getVertex1().getSequenceIndex());
			}
		}
		ConsensusBuilderBidirectionalSimple consensus = new ConsensusBuilderBidirectionalSimple();
		consensus.setSequenceNamePrefix("SuperContig");
		List<QualifiedSequence> answer = new ArrayList<QualifiedSequence>();
		answer.addAll(consensus.makeConsensus(graph));
		for(int k=0;k<contigsForGraph.size();k++) {
			if(!sequenceIdsInPaths.contains(k)) answer.add(contigsForGraph.get(k));
		}
		answer.addAll(smallContigs);
		return answer;
	}

	private int filterHits(int queryIdx, Map<Integer, List<UngappedSearchHit>> hits) {
		int maxCount = 0;
		Map<Integer,Integer> subjectIdxsCounts = new HashMap<Integer, Integer>();
		for(Map.Entry<Integer, List<UngappedSearchHit>> entry:hits.entrySet()) {
			if (entry.getKey()==queryIdx) continue;
			List<UngappedSearchHit> hitsSubject = entry.getValue();
			int count = getQueryIdxsCount(hitsSubject);
			subjectIdxsCounts.put(entry.getKey(),count);
			maxCount = Math.max(maxCount, count);
		}
		double limit = 0.5*maxCount;
		for(Map.Entry<Integer, Integer> entry:subjectIdxsCounts.entrySet()) {
			int count = entry.getValue();
			if(count<limit) hits.remove(entry.getKey());
		}
		return maxCount;
	}

	private int getQueryIdxsCount (List<UngappedSearchHit> hits) {
		Set<Integer> queryStarts = new HashSet<Integer>();
		for(UngappedSearchHit hit:hits) queryStarts.add(hit.getQueryIdx());
		return queryStarts.size();
	}
	private void buildEdges(AssemblyGraph graph, int queryEndIdx,String queryEnd, Map<Integer,String> contigEndsMap, Map<Integer, List<UngappedSearchHit>> hits, boolean revQuery, int countLimit) {
		int debugIdx = -1;
		//int debugIdx = 48;
		int querySeqId = queryEndIdx/2;
		int queryLength = graph.getSequenceLength(querySeqId);
		boolean queryStartSeq = (queryEndIdx%2==0);
		for(Map.Entry<Integer, List<UngappedSearchHit>> entry:hits.entrySet()) {
			int subjectEndIdx = entry.getKey();
			int subjectSeqId = subjectEndIdx/2;
			
			if(subjectSeqId>=querySeqId) continue;
			String subjectEnd = contigEndsMap.get(subjectEndIdx);
			List<UngappedSearchHit> hitsSubject = entry.getValue();
			if(queryEndIdx == debugIdx) System.err.println("Query: "+querySeqId+" subject: "+subjectSeqId+" end: "+subjectEndIdx+" qstart: "+queryStartSeq+" hits: "+hitsSubject.size());
			if (hitsSubject.size() < 20) continue;
			List<UngappedSearchHitsCluster> clusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(END_LENGTH, END_LENGTH, hitsSubject, 0);
			Collections.sort(clusters, (c1,c2)->c2.getNumDifferentKmers()-c1.getNumDifferentKmers());
			if(queryEndIdx == debugIdx) {
				System.err.println("Query: "+querySeqId+" subject end: "+subjectEndIdx+" clusters: "+clusters.size());
				for(UngappedSearchHitsCluster cluster:clusters) System.err.println("Query: "+querySeqId+" subject end: "+subjectEndIdx+" next cluster "+cluster.getPredictedOverlap()+" estimated: "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" kmers: "+cluster.getNumDifferentKmers());
			}
			for(int i=0;i<clusters.size() && i<2;i++) {
				UngappedSearchHitsCluster cluster = clusters.get(i);
				if(queryEndIdx == debugIdx) System.err.println("Query: "+querySeqId+" subject: "+subjectSeqId+" end: "+subjectEndIdx+" cluster: "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evSub "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" kmers "+cluster.getNumDifferentKmers());
				int numKmers = cluster.getNumDifferentKmers(); 
				if(numKmers<countLimit) break;
				ReadAlignment aln = aligner.buildCompleteAlignment(subjectEndIdx, subjectEnd, queryEnd, cluster);
				if(aln==null) continue;
				if(queryEndIdx == debugIdx) System.err.println("Query: "+querySeqId+" subject: "+subjectSeqId+" end: "+subjectEndIdx+" mismatches: "+aln.getNumMismatches()+" indel length: "+aln.getTotalLengthIndelCalls()+" CSK "+aln.getCoverageSharedKmers()+" overlap: "+cluster.getPredictedOverlap()+" alignment: "+aln);
				int overlap1 = cluster.getPredictedOverlap();
				if(aln.getTotalLengthIndelCalls()>0.01*overlap1) continue;
				if(aln.getCoverageSharedKmers()<0.3*overlap1) continue;	
				int softClipStart = aln.getSoftClipStart();
				int softClipEnd = aln.getSoftClipEnd();
				
				if(subjectEndIdx%2==0 && queryStartSeq==revQuery && softClipStart>1000 && softClipEnd<100) {
					//Query before subject
					AssemblyVertex subjectVertex = graph.getVertex(subjectSeqId, true);
					AssemblyVertex queryVertex = graph.getVertex(querySeqId, queryStartSeq);
					int overlap = aln.getLast();
					AssemblyEdge edge = new AssemblyEdge(queryVertex, subjectVertex, overlap);
					edge.setVertex1EvidenceStart(aln.getAlignedReadPosition(aln.getFirst()));
					edge.setVertex1EvidenceEnd(aln.getAlignedReadPosition(aln.getLast()));
					edge.setVertex2EvidenceStart(aln.getFirst()-1);
					edge.setVertex2EvidenceEnd(aln.getLast());
					edge.setCoverageSharedKmers(aln.getCoverageSharedKmers());
					edge.setWeightedCoverageSharedKmers(aln.getWeightedCoverageSharedKmers());
					edge.setNumMismatches(aln.getNumMismatches());
					edge.setNumIndels(aln.getTotalLengthIndelCalls());
					System.err.println("Adding edge. Num mismatches: "+aln.getNumMismatches()+" num indels: "+aln.getTotalLengthIndelCalls()+" edge: "+edge);
					graph.addEdge(edge);
				}
				if(subjectEndIdx%2==1 && queryStartSeq!=revQuery && softClipStart<100 && softClipEnd>1000) {
					//Query after subject
					AssemblyVertex subjectVertex = graph.getVertex(subjectSeqId, false);
					AssemblyVertex queryVertex = graph.getVertex(querySeqId, queryStartSeq);
					int overlap = 100000-aln.getFirst();
					AssemblyEdge edge = new AssemblyEdge(subjectVertex, queryVertex, overlap);
					edge.setVertex1EvidenceStart(aln.getFirst()-1);
					edge.setVertex1EvidenceEnd(aln.getLast());
					edge.setVertex2EvidenceStart(aln.getAlignedReadPosition(aln.getFirst()));
					edge.setVertex2EvidenceEnd(aln.getAlignedReadPosition(aln.getLast())+1);
					edge.setCoverageSharedKmers(aln.getCoverageSharedKmers());
					edge.setWeightedCoverageSharedKmers(aln.getWeightedCoverageSharedKmers());
					edge.setNumMismatches(aln.getNumMismatches());
					edge.setNumIndels(aln.getTotalLengthIndelCalls());
					System.err.println("Adding edge. Num mismatches: "+aln.getNumMismatches()+" num indels: "+aln.getTotalLengthIndelCalls()+" edge: "+edge);
					graph.addEdge(edge);
				}
			}
			
		}
		
	}

}
