package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;

public class ConsensusBuilderBidirectionalSimple implements ConsensusBuilder {
	int match;
	int openGap;
	int extGap;
	int mismatch;
	int windowSize;
	int tolerance;
	AlignmentAffineGap aligner;
	boolean startConsensus = true;
	
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
		AssemblyVertex lastVertex = null;
		String pathS = "";
		for(int j = 0; j < path.size(); j++)
		{
			//Needed to find which is the origin vertex
			AssemblyEdge edge = path.get(j);
			AssemblyVertex a = edge.getVertex1();
			AssemblyVertex b = edge.getVertex2();
			//If the first edge is being checked, compare to the second edge to find the origin vertex
			if(j == 0)
			{
				AssemblyEdge nextEdge = path.get(j + 1);
				//The common vertex is the second vertex of the path
				if(nextEdge.getVertex1().getIndex() == edge.getVertex1().getIndex() || nextEdge.getVertex2().getIndex() == edge.getVertex1().getIndex())
				{
					a = edge.getVertex2();
					b = edge.getVertex1();
				}
			}
			else if(lastVertex == b)
			{
				//The common vertex is the first vertex to be compared
				a = edge.getVertex2();
				b = edge.getVertex1();
			}
			if(j > 0 && lastVertex !=a) 
			{
				throw new RuntimeException("Inconsistency found in path");
			}
			if(j == 0) 
			{
				pathS = pathS.concat(a.getIndex() + ",");
				consensus.append(a.isStart() ? a.getRead().toString(): reverseComplement(a.getRead().toString()));
			} 
			else if(a.getRead()!=b.getRead())
			{
				//If the second string isn't start, then the reverse complement is added to the consensus
				String nextSequence = b.isStart() ? b.getRead().toString(): reverseComplement(b.getRead().toString());
					
				if(nextSequence.length() - edge.getOverlap() > tolerance) 
				{
					pathS = pathS.concat(b.getIndex() + ",");
					String overlapSegment = nextSequence.substring(0, edge.getOverlap());
					String remainingSegment = nextSequence.substring(edge.getOverlap());
					consensus.append(remainingSegment);
				} 
				else 
				{
					System.err.println("Non embedded edge has overlap: "+edge.getOverlap()+ " and length: "+nextSequence.length());
				}
			}
			lastVertex = b;
		}
		System.out.println(pathS);
		return consensus;
	}
	
	private String reverseComplement(String s)
	{
		StringBuilder complementaryStrand = new StringBuilder();
		for(int i = 0; i < s.length(); i++)
		{
			complementaryStrand.append(complementaryBase(s.charAt(i)));
		}
		return complementaryStrand.reverse().toString();
	}
	
	private char complementaryBase(char b)
	{
		char complementaryBase;
		if(b == 'A')
			complementaryBase = 'T';
		else if(b == 'T')
			complementaryBase = 'A';
		else if(b == 'C')
			complementaryBase = 'G';
		else if(b == 'G')
			complementaryBase = 'C';
		else if (b == '-')
			complementaryBase = '-';
		else
			complementaryBase = 'N';
		return complementaryBase;
	}
}
