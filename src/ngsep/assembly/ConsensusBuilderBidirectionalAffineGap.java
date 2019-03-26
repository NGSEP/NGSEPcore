package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;

public class ConsensusBuilderBidirectionalAffineGap implements ConsensusBuilder {
	int match = 2;
	int openGap = 20;
	int extGap = 1;
	int mismatch = 8;
	boolean startConsensus = true;
	
	@Override
	public List<CharSequence> makeConsensus(AssemblyGraph graph) 
	{
		List<CharSequence> consensusList = new ArrayList<CharSequence>();
		for(int i = 0; i < graph.getPaths().size(); i++)
		{
			List<AssemblyEdge> path = graph.getPaths().get(i);
			String consensus = "";
			startConsensus = true;
			for(int j = 0; j < path.size(); j++)
			{
				AssemblyEdge previousEdge = null;
				if(j > 0)
					previousEdge = path.get(j - 1);
				AssemblyEdge edge = path.get(j);
				AssemblyVertex a = edge.getVertex1();
				AssemblyVertex b = edge.getVertex2();
				if(previousEdge == null && path.size() > j + 1)
				{
					previousEdge = path.get(j + 1);
					if(previousEdge.getVertex1().getIndex() == edge.getVertex1().getIndex() || previousEdge.getVertex2().getIndex() == edge.getVertex1().getIndex())
					{
						a = edge.getVertex2();
						b = edge.getVertex1();
					}
				}
				else if (previousEdge != null)
				{
					if(previousEdge.getVertex1().getIndex() == edge.getVertex2().getIndex() || previousEdge.getVertex2().getIndex() == edge.getVertex2().getIndex())
					{
						a = edge.getVertex2();
						b = edge.getVertex1();
					}
				}
				String s1 = a.isStart() ? a.getRead().toString() : complementaryStrand(a.getRead().toString());
				String s2 = b.isStart() ? b.getRead().toString() : complementaryStrand(b.getRead().toString());
				AlignmentAffineGap alignment = new AlignmentAffineGap(match, openGap, extGap, mismatch);
				String[] alignments = alignment.getAlignment(s1, s2);	
				System.out.println(alignments[0]);
				System.out.println(alignments[1]);
				//consensus = consensus.concat(joinedString(graph.getEmbedded(a.getIndex()), graph.getEmbedded(j+1), alignmentOrig));
			}
			consensusList.add(consensus);
		}	
		return consensusList;
	}
	
	private String joinedString(List<AssemblyEmbedded> embedded1, List<AssemblyEmbedded> embedded2, String[] alignment)
	{
		StringBuilder finalString = new StringBuilder();
		int i = 0;
		if(!startConsensus)
		{
			for(i = alignment[0].length() - 1; i >= 0; i--)
			{
				char b = alignment[0].charAt(i);
				if(b != '-')
					break;
			}
			i++;
		}
		else
			startConsensus = false;
		
		for (; i < alignment[0].length(); i++) 
		{
			char a = alignment[0].charAt(i);
			char b = alignment[1].charAt(i);
			if(a == '-' && b != '-')
			{
				finalString.append(b);
			}
			else if(a != '-' && b == '-')
			{
				finalString.append(a);
			}
			else if(a == b)
			{
				finalString.append(a);
			}
			else if(a != b)
			{
				if(embedded1 != null && embedded2 != null)
					finalString.append(embedded1.size() > embedded2.size() ? a : b);
				else if (embedded1 != null)
					finalString.append(a);
				else 
					finalString.append(b);
			}
		}
		boolean stop = false;
		return finalString.toString();
	}
	
	private void printAlignmentMatrix(int[][] matrix, String s1, String s2)
	{
		System.out.print("\t-\t");
		for (int i = 0; i < s2.length(); i++) {
			System.out.print(s2.charAt(i) + "\t");
		}
		System.out.println();
		for (int i = 0; i < matrix.length; i++) {
			if(i == 0)
				System.out.print("-\t");
			else 
				System.out.print(s1.charAt(i - 1) + "\t");
		    for (int j = 0; j < matrix[i].length; j++) {
		        System.out.print(matrix[i][j] + "\t");
		    }
		    System.out.println();
		}
	}
	
	private String complementaryStrand(String s)
	{
		StringBuilder complementaryStrand = new StringBuilder();
		for(int i = 0; i < s.length(); i++)
		{
			complementaryStrand.append(complementaryBase(s.charAt(i)));
		}
		return complementaryStrand.toString();
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
