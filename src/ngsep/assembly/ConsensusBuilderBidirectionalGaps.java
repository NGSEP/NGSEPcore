package ngsep.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.alignments.ReadAlignment;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;


public class ConsensusBuilderBidirectionalGaps implements ConsensusBuilder {
	int match;
	int gap;
	int mismatch;
	boolean startConsensus = true;
	
	public ConsensusBuilderBidirectionalGaps() 
	{
		match = 1;
		mismatch = -1;
		gap = -2;		
	}
	
	@Override
	public List<CharSequence> makeConsensus(AssemblyGraph graph) 
	{
		//List of final contigs
		List<CharSequence> consensusList = new ArrayList<CharSequence>();
		List<List<AssemblyEdge>> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			List<AssemblyEdge> path = paths.get(i);
			CharSequence consensusPath = makeConsensus (graph, path);
			consensusList.add(consensusPath);
		}
		
		return consensusList;
	}
	
	private CharSequence makeConsensus(AssemblyGraph graph, List<AssemblyEdge> path) 
	{
		StringBuilder consensus = new StringBuilder();
		List<byte[]> consensusCounts = new ArrayList<byte[]>();
		AssemblyVertex lastVertex = null;
		String pathS = "";
		if(path.size()==1) {
			consensus.append(path.get(0).getVertex1().getRead());
			return consensus;
		}
		int currentPos = 0;
		for(int j = 0; j < path.size(); j++)
		{
			//Needed to find which is the origin vertex
			AssemblyEdge edge = path.get(j);
			AssemblyVertex vertexPreviousEdge;
			AssemblyVertex vertexNextEdge;
			//If the first edge is being checked, compare to the second edge to find the origin vertex
			if(j == 0)
			{
				AssemblyEdge nextEdge = path.get(j + 1);
				vertexNextEdge = edge.getSharedVertex(nextEdge);
				if(vertexNextEdge== null) throw new RuntimeException("Inconsistency found in first edge of path");
				vertexPreviousEdge = edge.getVertex1();
				if(vertexPreviousEdge == vertexNextEdge) vertexPreviousEdge = edge.getVertex2();
			}
			else if (lastVertex == edge.getVertex1())
			{
				vertexPreviousEdge = edge.getVertex1();
				vertexNextEdge = edge.getVertex2();
			}
			else if (lastVertex == edge.getVertex2())
			{
				vertexPreviousEdge = edge.getVertex2();
				vertexNextEdge = edge.getVertex1();
			}
			else 
			{
				throw new RuntimeException("Inconsistency found in path");
			}
			
			
			//If it's the first edge in the path, add the first read base count to the consensus
			if(j == 0) 
			{
				pathS = pathS.concat(vertexPreviousEdge.getIndex() + ",");
				int seqId = vertexPreviousEdge.getIndex();
				String seq = vertexPreviousEdge.getRead().toString();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
				
				List<AssemblyEmbedded> embeddedReads = graph.getEmbedded(seqId);
				Map<Integer,byte[]> consensusEmbedded = calculateEmbeddedConsensusCounts(seq, embeddedReads, reverse);
				for(int i = 0; i < seq.length(); i++) {
					byte [] countsPos;
					byte [] countsE = consensusEmbedded.get(i);
					if(countsE==null) countsPos = new byte[DNASequence.BASES_STRING.length()];
					else countsPos = Arrays.copyOf(countsE, countsE.length);
					consensusCounts.add(countsPos);
					int idx = DNASequence.BASES_STRING.indexOf(seq.charAt(i));
					if(idx>=0) countsPos[idx]++;
					
				}
				currentPos = seq.length();
			} else if(vertexPreviousEdge.getRead()!=vertexNextEdge.getRead()) {
				int seqId = vertexNextEdge.getIndex();
				List<AssemblyEmbedded> embeddedReads = graph.getEmbedded(seqId);
				
				String seq = vertexNextEdge.getRead().toString();
				boolean reverse = !vertexNextEdge.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
				
				Map<Integer,byte[]> consensusEmbedded = calculateEmbeddedConsensusCounts(seq, embeddedReads, reverse);
				//TODO: Actually build alignment
				currentPos = currentPos - edge.getOverlap();
				for(int i = 0; i < seq.length(); i++) {
					byte [] countsPos;
					byte [] countsE = consensusEmbedded.get(i);
					if(countsE==null) countsPos = new byte[DNASequence.BASES_STRING.length()];
					else countsPos = countsE;
					byte [] countsTotal;
					if(consensusCounts.size() == currentPos) {
						countsTotal = Arrays.copyOf(countsPos, countsPos.length);
						consensusCounts.add(countsTotal);
					}
					else {
						countsTotal = consensusCounts.get(currentPos);
						for(int k = 0; k < countsTotal.length; k++) {
							//Add the base count of the embedded count to the global consensus for each base position
							countsTotal[k] += countsPos[k];
						}
					}
					int idx = DNASequence.BASES_STRING.indexOf(seq.charAt(i));
					if(idx>=0) countsTotal[idx]++;
					currentPos++;
				}
				
				if(seq.length() <= edge.getOverlap()) 
				{
					System.err.println("Non-embedded edge has overlap: " + edge.getOverlap() + " and length: " + seq.length());
				} 
			}
			lastVertex = vertexNextEdge;
		}
		
		//Find the most common occurrence in each position and add the most common character to the consensus
		for(byte[] counts : consensusCounts)
		{
			int indexMax = NumberArrays.getIndexMaximum(counts);
			if(indexMax>=0 && counts[indexMax]>0) {
				char base = DNASequence.BASES_STRING.charAt(indexMax);
				consensus.append(base);
			}
		}
		return consensus.toString();
	}
	
	private Map<Integer,byte[]> calculateEmbeddedConsensusCounts(String read, List<AssemblyEmbedded> embeddedReads, boolean reverse) {
		Map<Integer,byte[]> consensusCounts = new TreeMap<Integer,byte[]>();
		//Align the read a with its embedded reads
		EmbeddedReadsAlignerConstantGap embeddedAligner = new EmbeddedReadsAlignerConstantGap(match, gap, mismatch);
		List<ReadAlignment> embeddedAlns = embeddedAligner.alignEmbeddedReads(read, embeddedReads, reverse);
		for(ReadAlignment aln:embeddedAlns) {
			String embeddedRead = aln.getReadCharacters().toString();
			//Get each character of the aligned read
			for(int l = 0; l < embeddedRead.length(); l++)
			{
				char embeddedChar = embeddedRead.charAt(l);
				int subjectPos = aln.getReferencePosition(l);
				if(subjectPos<0) continue;
				byte [] countsPos = consensusCounts.get(subjectPos);
				if(countsPos==null) {
					countsPos = new byte[DNASequence.BASES_STRING.length()];
					consensusCounts.put(subjectPos,countsPos);
				}
				int idx = DNASequence.BASES_STRING.indexOf(embeddedChar);
				if(idx>=0) countsPos[idx]++;
			}
		}
		return consensusCounts;
	}
}
