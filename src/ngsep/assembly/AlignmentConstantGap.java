package ngsep.assembly;

import java.util.Arrays;

public class AlignmentConstantGap {
	public static final char GAP_CHARACTER = '-';
	private int match;
	private int gap;
	private int mismatch;
	
	public AlignmentConstantGap(int match, int gap, int mismatch) 
	{
		this.match = match;           
		this.mismatch = mismatch;   
		this.gap = gap;             
	}
	
	
	public String[] getAlignment(String s1, String s2)
	{
		int[][] matrixOrig = alignmentMatrixConstantGap(s1, s2);
		String[] alignmentOrig = sequencesAlignment(matrixOrig, s1, s2);
		return alignmentOrig;
	}
	
	private String[] sequencesAlignment(int[][] matrix, String s1, String s2)
	{
		String[] sequences = {"", ""};
		int i = s1.length();
		int j = s2.length();
		while(i > 0 && j > 0)
		{
			int current = matrix[i][j];
			int diag = matrix[i - 1][j - 1];
			int left = matrix[i - 1][j];
			int up = matrix[i][j - 1];
			
			if( (diag == current - mismatch && s1.charAt(i - 1) != s2.charAt(j - 1)) 
					|| (diag == current - match && s1.charAt(i - 1) == s2.charAt(j - 1)))
			{
				sequences[0] = s1.charAt(i - 1) + sequences[0];
				sequences[1] = s2.charAt(j - 1) + sequences[1];
				i--;
				j--;
			}
			else if (up == current - gap)
			{
				sequences[0] = GAP_CHARACTER + sequences[0];
				sequences[1] = s2.charAt(j - 1) + sequences[1];
				j--;
			}
			else if( left == current - gap)
			{
				sequences[0] = s1.charAt(i - 1) + sequences[0];
				sequences[1] = GAP_CHARACTER + sequences[1];
				i--;
			}
			else
			{
				throw new RuntimeException("Inconsistency in alignment matrix");
			}
		}
		if(j > 0)
		{
			char[] gaps = new char[j];
			Arrays.fill(gaps, GAP_CHARACTER);
			sequences[0] = new String(gaps) + sequences[0];
			sequences[1] = s2.substring(0, j) + sequences[1];
		}
		else if (i > 0)
		{
			char[] gaps = new char[i];
			Arrays.fill(gaps, GAP_CHARACTER);
			sequences[1] = new String(gaps) + sequences[1];
			sequences[0] = s1.substring(0, i) + sequences[0];
		}
		//printAlignmentMatrix(matrix, s1, s2);
		
		return sequences;
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
				int diag = matrix[i-1][j-1] + (s1.charAt(i - 1) == s2.charAt(j - 1) ? match : mismatch);
				int left = matrix[i-1][j] + gap;
				int up = matrix[i][j-1] + gap;
				matrix[i][j] = Math.max(Math.max(diag, left), up);
			}
		}		
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
}
