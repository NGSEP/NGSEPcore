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
	
	public ConsensusBuilderBidirectionalSimple(int match, int openGap, int extGap, int mismatch, int windowSize, double Rate_of_changes, double Rate_of_cuts, double Rate_of_cover) 
	{
		double rate_of_error = Rate_of_changes + Rate_of_cuts - Rate_of_changes * Rate_of_cuts;
		this.match = match;
		this.openGap = openGap;
		this.extGap = extGap;
		this.mismatch = mismatch;
		this.windowSize = windowSize;
		this.tolerance = (int) (Rate_of_cuts * (11.51292546 / rate_of_error));;
		aligner = new AlignmentAffineGap(match, openGap, extGap, mismatch);
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
		AssemblyVertex lastVertex = null;
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
				consensus.append(a.isStart() ? a.getRead().toString(): reverseComplement(a.getRead().toString()));
			} 
			else if(a.getRead()!=b.getRead())
			{
				//If the second string isn't start, then the reverse complement is added to the consensus
				String nextSequence = b.isStart() ? b.getRead().toString(): reverseComplement(b.getRead().toString());
					
				if(nextSequence.length() - edge.getOverlap() > tolerance) 
				{
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