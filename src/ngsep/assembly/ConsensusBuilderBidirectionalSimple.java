package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;

public class ConsensusBuilderBidirectionalSimple implements ConsensusBuilder {
	int tolerance;
	
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
		if(path.size()==1) {
			consensus.append(path.get(0).getVertex1().getRead());
			return consensus;
		}
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
			if(j == 0) 
			{
				pathS = pathS.concat(vertexPreviousEdge.getUniqueNumber() + ",");
				CharSequence seq = vertexPreviousEdge.getRead();
				boolean reverse = !vertexPreviousEdge.isStart();
				if(reverse) seq = DNAMaskedSequence.getReverseComplement(seq);
			} 
			else if(vertexPreviousEdge.getRead()!=vertexNextEdge.getRead())
			{
				//If the second string isn't start, then the reverse complement is added to the consensus
				CharSequence seq = vertexNextEdge.getRead();
				boolean reverse = !vertexNextEdge.isStart();
				if(reverse) seq = DNASequence.getReverseComplement(seq);
				if(seq.length() - edge.getOverlap() > tolerance) 
				{
					pathS = pathS.concat(vertexNextEdge.getUniqueNumber() + ",");
					//String overlapSegment = nextSequence.substring(0, edge.getOverlap());
					String remainingSegment = seq.subSequence(edge.getOverlap(),seq.length()).toString();
					//if (consensus.length()>490000 && consensus.length()<510000) System.out.println("Consensus length: "+consensus.length()+" Vertex: "+vertexNextEdge.getUniqueNumber()+" read length: "+seq.length()+" overlap: "+edge.getOverlap()+" remaining: "+remainingSegment.length());
					consensus.append(remainingSegment);
				} 
				else 
				{
					System.err.println("Non embedded edge has overlap: "+edge.getOverlap()+ " and length: "+seq.length());
				}
			}
			lastVertex = vertexNextEdge;
		}
		System.out.println(pathS);
		return consensus;
	}
}
