package ngsep.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.sun.xml.internal.ws.util.StringUtils;

public class ConsensusBuilderBidirectionalAffineGap implements ConsensusBuilder {
	int match = 2;
	int openGap = -20;
	int extGap = -1;
	int mismatch = -8;
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
			for(AssemblyEdge edge:path)
			{
				// FIXME: Look for path direction to check from which to which vertex
				AssemblyVertex a = edge.getVertex1();
				AssemblyVertex b = edge.getVertex2();
				String s1 = a.getRead().toString();
				String s2 = b.getRead().toString();
				//Align original
				List<int[][]> matrixOrig = alignmentMatrixAffineGap(s1, s2);
				String[] alignmentOrig = sequencesAlignment(matrixOrig, s1, s2);				
				//consensus = consensus.concat(joinedString(graph.getEmbedded(a.getIndex()), graph.getEmbedded(j+1), alignmentOrig));
			}
			consensusList.add(consensus);
		}	
		return consensusList;
	}
	
	private String joinedString(List<AssembyEmbedded> embedded1, List<AssembyEmbedded> embedded2, String[] alignment)
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
	
	private int match(char a, char b)
	{
		if (a == b)
			return match;
		else
			return mismatch;
	}
	
	private int initX(int i, int j)
	{
		if (i > 0 && j == 0)
			return Integer.MIN_VALUE;
		else
		{
			if (j > 0)
				return openGap + (extGap * j);
			else
				return 0;
		}
	}
	
	private int initY(int i, int j)
	{
		if (j > 0 && i == 0)
			return Integer.MIN_VALUE;
		else
		{
			if (i > 0)
				return openGap + (extGap * i);
			else
				return 0;
		}
	}
	
	private int initM(int i, int j)
	{
		if (i == 0 && j == 0)
			return 0;
		else
		{
			if (j == 0 || i == 0)
				return Integer.MIN_VALUE;
			else
				return 0;
		}
	}
	
	private List<int[][]> alignmentMatrixAffineGap(String s1, String s2)
	{
		List<int[][]> matrices = new ArrayList<int[][]>();
		int[][] x = new int[s1.length() + 1][s2.length() + 1];
		int[][] y = new int[s1.length() + 1][s2.length() + 1];
		int[][] m = new int[s1.length() + 1][s2.length() + 1];
		
		for(int i = 0; i < x.length; i++)
		{
			for (int j = 0; j < x[0].length; j++)
			{
				x[i][j] = initX(i, j);
				y[i][j] = initY(i, j);
				m[i][j] = initM(i, j);
			}
		}
		
		for(int i = 1; i < s1.length(); i++)
		{
			for(int j = 1; j < s2.length(); j++)
			{
				x[i][j] = Math.max(openGap + extGap + m[i][j-1], Math.max(extGap + x[i][j-1], openGap + extGap + y[i][j-1]));
				y[i][j] = Math.max(openGap + extGap + m[i-1][j], Math.max(openGap + extGap + x[i-1][j], extGap + y[i-1][j]));
				m[i][j] = Math.max(match(s1.charAt(i - 1), s2.charAt(j - 1)) + m[i-1][j-1], Math.max(x[i][j], y[i][j]));
			}
		}
		matrices.add(x);
		matrices.add(y);
		matrices.add(m);
		return matrices;
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
	
	private String[] sequencesAlignment(List<int[][]> matrices, String s1, String s2)
	{	
		String[] sequences = {"", ""};
		int i = s1.length();
		int j = s2.length();
		int[][] x = matrices.get(0);
		int[][] y = matrices.get(1);
		int[][] m = matrices.get(2);
		while(i > 0 || j > 0)
		{
			System.out.println(j);
			if(i > 0 && j > 0 && m[i][j] == (m[i-1][j-1] + match(s1.charAt(i - 1), s2.charAt(j - 1))))
			{
				System.out.println("A");
				sequences[0] += s1.charAt(i - 1);
				sequences[1] += s2.charAt(j - 1);
				i--;
				j--;
			}
			else if (j > 0 && m[i][j] == y[i][j])
			{
				System.out.println("B");
				sequences[0] += '-';
				sequences[1] += s2.charAt(j-1);
				i--;
			}
			else if (j > 0 && m[i][j] == y[i][j])
			{
				System.out.println("C");
				sequences[0] += s1.charAt(i-1);
				sequences[1] += '-';
				j--;
			}
		}
		System.out.println(sequences[0]);
		System.out.println(sequences[1]);
		System.out.println();
		return sequences;
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
