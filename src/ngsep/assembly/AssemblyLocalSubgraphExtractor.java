package ngsep.assembly;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.assembly.io.AssemblyGraphFileHandler;
import ngsep.sequences.QualifiedSequence;

public class AssemblyLocalSubgraphExtractor {

	private List<Integer> selectedSequences;
	private int maxLevel = 2;
	public static void main(String[] args) throws Exception {
		AssemblyLocalSubgraphExtractor instance = new AssemblyLocalSubgraphExtractor();
		instance.loadSelectedSequences(args[1]);
		List<QualifiedSequence> sequences = AssemblyGraphFileHandler.loadSequenceNamesFromGraphFile(args[0]);
		AssemblyGraph graph = AssemblyGraphFileHandler.load(sequences, args[0]);
		updateSequenceNames (graph, args[2]);
		instance.buildSubgraphs(graph, args[3]);	
	}
	
	private void loadSelectedSequences(String filename) throws IOException {
		selectedSequences = new ArrayList<Integer>();
		try(FileReader reader = new FileReader(filename);
			BufferedReader in = new BufferedReader(reader)) {
			String line = in.readLine();
			while(line!=null) {
				selectedSequences.add(Integer.parseInt(line));
				line = in.readLine();
			}
		}
	}
	
	private static void updateSequenceNames(AssemblyGraph graph, String alignmentsFile) throws IOException {
		List<QualifiedSequence> sequences = graph.getSequences();
		Map<String,Integer> seqIds = new HashMap<String, Integer>();
		for(int i=0;i<graph.getNumSequences();i++) {
			QualifiedSequence seq = graph.getSequence(i);
			if(seqIds.containsKey(seq.getName())) {
				System.err.println("Duplicated read name: "+seq.getName());
			} else {
				seqIds.put(seq.getName(), i);
			}
		}
		try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(alignmentsFile)) {
			reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_ALIGNMENT_NAME);
			reader.setFilterFlags(ReadAlignment.FLAG_SECONDARY);
			Iterator<ReadAlignment> it = reader.iterator();
			while (it.hasNext()) {
				ReadAlignment aln = it.next();
				if(aln.getReadLength()<Assembler.DEF_MIN_READ_LENGTH) continue;
				if(aln.isReadUnmapped()) continue;
				Integer idx = seqIds.get(aln.getReadName());
				if(idx==null) {
					System.err.println("Aligned read: "+aln.getReadName()+" not found in graph");
					continue;
				}
				else aln.setReadNumber(idx);
				
				
				//sequences.add(new QualifiedSequence(aln.getReadName(), characters));
				
				String newReadName = aln.getSequenceName()+"_"+aln.getFirst()+"_"+(aln.isNegativeStrand()?1:0);
				//aln.setReadName(newReadName);
				sequences.get(idx).setName(newReadName);
			}
		}
	}
	private void buildSubgraphs(AssemblyGraph graph, String outPrefix) throws IOException {
		AssemblyGraph unfiltered = buildSubgraph(graph);
		AssemblyGraphFileHandler.save(unfiltered, outPrefix+"_unfiltered.graph.gz");
		System.out.println("Unfiltered subgraph");
		//logSubgraph(unfiltered);
		graph.removeVerticesChimericReads();
		graph.updateScores();
		(new AssemblySequencesRelationshipFilter()).filterEdgesAndEmbedded(graph, 0.3);
		AssemblyGraph filtered = buildSubgraph(graph);
		AssemblyGraphFileHandler.save(filtered, outPrefix+"_filtered.graph.gz");
		System.out.println("Filtered subgraph");
		logSubgraph(filtered);
	}
	private void logSubgraph(AssemblyGraph graph) {
		System.out.println("EDGES");
		for(AssemblyEdge edge:graph.getEdges()) System.out.println(edge);
		System.out.println();
		System.out.println("EMBEDDED");
		for(AssemblyEmbedded embedded:graph.getAllEmbedded()) {
			System.out.println(embedded);
		}
		
	}
	private AssemblyGraph buildSubgraph(AssemblyGraph graph) {
		Set<Integer> allSeqIds = new HashSet<Integer>();
		for(int sourceSeqId:selectedSequences) {
			allSeqIds.add(sourceSeqId);
			Map<Integer,Integer> levels = new HashMap<Integer, Integer>();
			Queue<Integer> agenda = new LinkedList<Integer>();
			agenda.add(sourceSeqId);
			levels.put(sourceSeqId,0);
			while(agenda.size()>0) {
				int next = agenda.remove();
				int level = levels.get(next);
				if(level == maxLevel) continue;
				AssemblyVertex v1 = graph.getVertex(next, true);
				if(v1==null) continue;
				List<AssemblyEdge> edges = graph.getEdges(v1);
				for(AssemblyEdge edge:edges) {
					int seqId1 = edge.getVertex1().getSequenceIndex();
					if(!levels.containsKey(seqId1)) {
						allSeqIds.add(seqId1);
						agenda.add(seqId1);
						levels.put(seqId1,level+1);
					}
					int seqId2 = edge.getVertex2().getSequenceIndex();
					if(!levels.containsKey(seqId2)) {
						allSeqIds.add(seqId2);
						agenda.add(seqId2);
						levels.put(seqId2,level+1);
					}
				}
				AssemblyVertex v2 = graph.getVertex(next, false);
				List<AssemblyEdge> edges2 = graph.getEdges(v2);
				for(AssemblyEdge edge:edges2) {
					int seqId1 = edge.getVertex1().getSequenceIndex();
					if(!levels.containsKey(seqId1)) {
						allSeqIds.add(seqId1);
						agenda.add(seqId1);
						levels.put(seqId1,level+1);
					}
					int seqId2 = edge.getVertex2().getSequenceIndex();
					if(!levels.containsKey(seqId2)) {
						allSeqIds.add(seqId2);
						agenda.add(seqId2);
						levels.put(seqId2,level+1);
					}
				}
			}
		}
		List<Integer> allSeqIdsList = new ArrayList<Integer>(allSeqIds);
		for(int seqId:allSeqIdsList) {
			List<AssemblyEmbedded> embRels = graph.getEmbeddedByHostId(seqId);
			for(AssemblyEmbedded embedded: embRels) allSeqIds.add(embedded.getSequenceId());
			List<AssemblyEmbedded> embRels2 = graph.getEmbeddedBySequenceId(seqId);
			for(AssemblyEmbedded embedded: embRels2) allSeqIds.add(embedded.getHostId());
		}
		AssemblyGraph subgraph = graph.buildSubgraph(allSeqIds);
		return subgraph;
	}

}
