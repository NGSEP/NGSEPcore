package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.io.FastaSequencesHandler;

public class ContigEndsMerger {

	private MinimizersTableReadAlignmentAlgorithm aligner = new MinimizersTableReadAlignmentAlgorithm();
	public static void main(String[] args) throws Exception {
		ContigEndsMerger instance = new ContigEndsMerger();
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList contigs = handler.loadSequences(args[0]);
		List<QualifiedSequence> answer = instance.mergeContigs(contigs);
		handler.saveSequences(answer, System.out, 100);
	}
	
	public List<QualifiedSequence> mergeContigs(List<QualifiedSequence> contigs) {
		
		MinimizersTable table = new MinimizersTable(15, 20);
		Map<Integer,String> contigEndsMap = new HashMap<Integer, String>();
		List<QualifiedSequence> contigsForGraph = new ArrayList<QualifiedSequence>(contigs.size());
		List<QualifiedSequence> smallContigs = new ArrayList<QualifiedSequence>(contigs.size());
		int i=0;
		for (QualifiedSequence contig:contigs) {
			CharSequence seq = contig.getCharacters();
			if(seq.length()<200000) {
				smallContigs.add(contig);
				continue;
			}
			contigsForGraph.add(contig);
			String startSeq = seq.subSequence(0, 100000).toString();
			contigEndsMap.put(2*i, startSeq);
			table.addSequence(2*i, startSeq);
			String endSeq = seq.subSequence(seq.length()-100000,seq.length()).toString();
			contigEndsMap.put(2*i+1, endSeq);
			table.addSequence(2*i+1, endSeq);
			i++;
		}
		if(contigsForGraph.size()<2) {
			System.err.println("Not enough large contigs to join");
			return contigs;
		}
		System.err.println("Built minimizers table");
		AssemblyGraph graph = new AssemblyGraph(contigsForGraph);
		for(Map.Entry<Integer, String> entry:contigEndsMap.entrySet()) {
			int j = entry.getKey();
			String seq = entry.getValue();
			Map<Integer,List<UngappedSearchHit>> hitsForward = table.match(j, seq);
			buildEdges(graph,j,seq,contigEndsMap, hitsForward,false);
			seq = DNAMaskedSequence.getReverseComplement(seq).toString();
			Map<Integer,List<UngappedSearchHit>> hitsReverse = table.match(j, seq);
			buildEdges(graph,j,seq,contigEndsMap, hitsReverse,true);
			System.err.println("Processed contig end: "+j);
		}
		LayoutBuilderGreedyMaxOverlap builder = new LayoutBuilderGreedyMaxOverlap();
		builder.findPaths(graph);
		List<List<AssemblyEdge>> paths = graph.getPaths();
		Set<Integer> sequenceIdsInPaths = new HashSet<Integer>();
		System.err.println("Built paths. Number: "+paths.size());
		for(List<AssemblyEdge> path:paths) {
			System.err.println("Path: "+path.size());
			for(AssemblyEdge edge:path) {
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

	private void buildEdges(AssemblyGraph graph, int queryEndIdx,String queryEnd, Map<Integer,String> contigEndsMap, Map<Integer, List<UngappedSearchHit>> hits, boolean revQuery) {
		int debugLength = -1;
		int querySeqId = queryEndIdx/2;
		int queryLength = graph.getSequenceLength(querySeqId);
		boolean queryStartSeq = (queryEndIdx%2==0);
		for(Map.Entry<Integer, List<UngappedSearchHit>> entry:hits.entrySet()) {
			int subjectEndIdx = entry.getKey();
			int subjectSeqId = subjectEndIdx/2;
			
			if(subjectSeqId>=querySeqId) continue;
			String subjectEnd = contigEndsMap.get(subjectEndIdx);
			List<UngappedSearchHit> hitsSubject = entry.getValue();
			if(queryLength == debugLength) System.err.println("Query: "+querySeqId+" subject: "+subjectSeqId+" end: "+subjectEndIdx+" qstart: "+queryStartSeq+" hits: "+hitsSubject.size());
			if (hitsSubject.size() < 20) continue;
			List<KmerHitsCluster> clusters = KmerHitsCluster.clusterRegionKmerAlns(100000, 100000, hitsSubject, 0);
			Collections.sort(clusters, (c1,c2)->c2.getNumDifferentKmers()-c1.getNumDifferentKmers());
			int maxNumDifKmer = -1;
			if(queryLength == debugLength) System.err.println("Query: "+querySeqId+" subject: "+subjectSeqId+" end: "+subjectEndIdx+" clusters: "+clusters.size());
			for(KmerHitsCluster cluster:clusters) {
				if(queryLength == debugLength) System.err.println("Query: "+querySeqId+" subject: "+subjectSeqId+" end: "+subjectEndIdx+" cluster: "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evSub "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" kmers "+cluster.getNumDifferentKmers());
				int numKmers = cluster.getNumDifferentKmers(); 
				if(numKmers<20) break;
				if(maxNumDifKmer==-1) maxNumDifKmer = numKmers;
				else if(numKmers<0.5*maxNumDifKmer) break;
				ReadAlignment aln = aligner.buildCompleteAlignment(subjectEndIdx, subjectEnd, queryEnd, cluster);
				if(aln==null) continue;
				if(queryLength == debugLength) System.err.println("Query: "+querySeqId+" subject: "+subjectSeqId+" end: "+subjectEndIdx+" mismatches: "+aln.getNumMismatches()+" CSK "+aln.getCoverageSharedKmers()+" overlap: "+cluster.getPredictedOverlap()+" alignment: "+aln);
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