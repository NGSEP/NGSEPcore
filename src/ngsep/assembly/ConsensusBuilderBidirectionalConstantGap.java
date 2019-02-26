package ngsep.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.sun.xml.internal.ws.util.StringUtils;

public class ConsensusBuilderBidirectionalConstantGap implements ConsensusBuilder {
	int match = 5;
	int gap = -2;
	int mismatch = -1;
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
				int[][] matrixOrig = alignmentMatrixConstantGap(s1, s2);
				String[] alignmentOrig = sequencesAlignment(matrixOrig, s1, s2);
				int scoreOrig = matrixOrig[s1.length()][s2.length()];
				//Align complementary
				String s3 = complementaryStrand(s2);
				int[][] matrixComp = alignmentMatrixConstantGap(s1, s3);
				String[] alignmentComp = sequencesAlignment(matrixComp, s1, s3);
				int scoreComp = matrixComp[s1.length()][s3.length()];
				scoreComp = Integer.MIN_VALUE;
				/*if(scoreOrig > scoreComp)
				{
					consensus = consensus.concat(joinedString(graph.getEmbedded(j), graph.getEmbedded(j+1), alignmentOrig));
				}
				else
				{
					consensus = consensus.concat(joinedString(graph.getEmbedded(j), graph.getEmbedded(j+1), alignmentComp));
				}*/
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
	
	private int[][] alignmentMatrixConstantGap(String s1, String s2)
	{
		int[][] matrix = new int[s1.length() + 1][s2.length() + 1];
		for(int i = 0; i < s1.length() + 1; i++)
		{
			matrix[i][0] = i * gap;
		}
		for(int i = 0; i < s2.length() + 1; i++)
		{
			matrix[0][i] = i * gap;
		}
		for(int i = 1; i < s1.length() + 1; i++)
		{
			for(int j = 1; j < s2.length() + 1; j++)
			{
				int diag = matrix[i-1][j-1] + (s1.charAt(i - 1) == s2.charAt(j - 1) ? match : gap);
				int left = matrix[i-1][j] + gap;
				int up = matrix[i][j-1] + gap;
				matrix[i][j] = Math.max(Math.max(diag, left), up);
			}
		}
		
		//printAlignmentMatrix(matrix, s1, s2);
		
		return matrix;
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
	
	private String[] sequencesAlignment(int[][] matrix, String s1, String s2)
	{
		String[] sequences = {"", ""};
		int i = s1.length();
		int j = s2.length();
		while(i > 0 && j > 0)
		{
			int diag = matrix[i - 1][j - 1];
			int left = matrix[i - 1][j];
			int up = matrix[i][j - 1];
			

			if (j > 0 && up >= diag && up >= left)
			{
				sequences[0] = "-" + sequences[0];
				sequences[1] = s2.charAt(j - 1) + sequences[1];
				j--;
			}
			else if(i > 0 && j > 0 && diag >= up && diag >= left)
			{
				sequences[0] = s1.charAt(i - 1) + sequences[0];
				sequences[1] = s2.charAt(j - 1) + sequences[1];
				i--;
				j--;
			}
			else if(i > 0 && left >= diag && left >= up)
			{
				sequences[0] = s1.charAt(i - 1) + sequences[0];
				sequences[1] = "-" + sequences[1];
				i--;
			}
			
		}
		if(i == 0 && j > 0)
		{
			char[] gaps = new char[j];
			Arrays.fill(gaps, '-');
			sequences[0] = new String(gaps) + sequences[0];
			sequences[1] = s1.substring(0, j) + sequences[1];
		}
		else if (i > 0 && j == 0)
		{
			char[] gaps = new char[i];
			Arrays.fill(gaps, '-');
			sequences[1] = new String(gaps) + sequences[1];
			sequences[0] = s1.substring(0, i) + sequences[0];
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
