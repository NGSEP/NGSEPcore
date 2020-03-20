/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.sequences;

/**
 * Performs pairwise alignment using the affine gap method.
 * Adapted from https://www.itu.dk/~sestoft/bsa/Match2.java
 * @author David Guevara
 */
public class PairwiseAlignmentAffineGap {
	
	int match;
	int openGap;
	int extGap;
	int mismatch;
	int[][] x;
	int[][] y;
	int[][] m;
	Traceback[][][] b;
	
	public PairwiseAlignmentAffineGap(int match, int openGap, int extGap, int mismatch) 
	{
		System.out.println("Creating pairwise aligner");
		this.match = match;
		this.openGap = openGap;
		this.extGap = extGap;
		this.mismatch = mismatch;
		x = new int [2100][2100];
		y = new int [2100][2100];
		m = new int [2100][2100];
		b = new Traceback [3][2100][2100];
	}
	
	public String[] getAlignment(String s1, String s2) 
	{		
		initMatrices(s1, s2);
	    calculateMatrices(s1, s2);
	    int k = getMaxK(s1.length(), s2.length());
    	Traceback root = new Traceback(k, s1.length(), s2.length());	    
        return getAlignedStrings(root, s1, s2);
	}
	
	private void initMatrices(String s1, String s2)
	{
		if(x==null || x.length<s1.length()+1 || x[0].length < s2.length() +1 ) {
			System.out.println("Resizing matrices to "+(s1.length() + 1)+" - "+(s2.length() +1));
			x = new int[s1.length() + 1][s2.length() + 1];
			y = new int[s1.length() + 1][s2.length() + 1];
			m = new int[s1.length() + 1][s2.length() + 1];
			b = new Traceback[3][s1.length() + 1][s2.length() + 1];
		}
		
		
		for (int i = 1; i < x.length; i++) 
		{
			x[i][0] = - openGap - extGap * (i - 1);
	    	y[i][0] = s1.length() * -openGap * 1000;
	    	m[i][0] = y[i][0];
	    	b[1][i][0] = new Traceback(1, i-1, 0);
	    }
	    for (int i = 1; i < x[0].length; i++) 
	    {
	    	y[0][i] = - openGap - extGap * (i - 1);
	        x[0][i] = s2.length() * -openGap * 1000;
	        m[0][i] = x[0][i];
	    	b[2][0][i] = new Traceback(2, 0, i-1);
	    }
	}
	
	private void calculateMatrices(String s1, String s2)
	{
		for (int i = 1; i <= s1.length(); i++)
	    {
	    	for (int j = 1; j <= s2.length(); j++)
	    	{
	    		int matchScore = getMatchScore(s1.charAt(i - 1), s2.charAt(j - 1));
	    		m[i][j] = Math.max(m[i-1][j-1] + matchScore, Math.max(x[i-1][j-1] + matchScore, y[i-1][j-1] + matchScore));
	    		if(m[i][j] == m[i-1][j-1] + matchScore)
	    		{
	    			b[0][i][j] = new Traceback(0, i-1, j-1);
	    		}
	    		else if(m[i][j] == x[i-1][j-1] + matchScore)
	    		{
	    			b[0][i][j] = new Traceback(1, i-1, j-1);
	    		}
	    		else if(m[i][j] == y[i-1][j-1] + matchScore)
	    		{
	    			b[0][i][j] = new Traceback(2, i-1, j-1);
	    		}
	    		x[i][j] = Math.max(m[i-1][j] - openGap, Math.max(x[i-1][j] - extGap, y[i-1][j] - openGap));
	    		if(x[i][j] == m[i-1][j] - openGap)
	    		{
	    			b[1][i][j] = new Traceback(0, i-1, j);
	    		}
	    		else if(x[i][j] == x[i-1][j] - extGap)
	    		{
	    			b[1][i][j] = new Traceback(1, i-1, j);
	    		}
	    		else if(x[i][j] == y[i-1][j] - openGap)
	    		{
	    			b[1][i][j] = new Traceback(2, i-1, j);
	    		}
	    		y[i][j] = Math.max(m[i][j-1] - openGap, Math.max(x[i][j-1] - openGap, y[i][j-1] - extGap));
	    		if(y[i][j] == m[i][j-1] - openGap)
	    		{
	    			b[2][i][j] = new Traceback(0, i, j-1);
	    		}
	    		else if(y[i][j] == x[i][j-1] - openGap)
	    		{
	    			b[2][i][j] = new Traceback(1, i, j-1);
	    		}
	    		else if(y[i][j] == y[i][j-1] - extGap)
	    		{
	    			b[2][i][j] = new Traceback(2, i, j-1);
	    		}
	    	}
	    }
//		System.out.println("X");
//		printAlignmentMatrix(x, s1, s2);
//		System.out.println("Y");
//		printAlignmentMatrix(y, s1, s2);
//		System.out.println("M");
//		printAlignmentMatrix(m, s1, s2);
	}

	private int getMatchScore(char a, char b)
	{
		if (a == b)
			return match;
		else
			return -mismatch;
	}
	
	private Traceback getNextTraceback(Traceback t)
	{
		return b[t.k][t.i][t.j];
	}
	
	private int getMaxK(int i, int j)
	{
		int k = 0;
	    int val = m[i][j];
    	if (val < x[i][j]) {
    		k = 1;
    		val = x[i][j];
    	}
    	if (val < y[i][j]) {
    		k = 2;
    	}
    	return k;
	}
	
	private String[] getAlignedStrings(Traceback root, String s1, String s2)
	{
		StringBuffer sb1 = new StringBuffer();
    	StringBuffer sb2 = new StringBuffer();
	    String[] seqs = new String[2];
        seqs[0] = "";
        seqs[1] = "";
        Traceback tb = root;
        int i = tb.i;
        int j = tb.j; 
        while((tb = getNextTraceback(tb)) != null) {
	        if (i == tb.i) {
	        	sb1.append(LimitedSequence.GAP_CHARACTER);
	        } else {
	        	sb1.append(s1.charAt(i - 1));
	        }	
	        if (j == tb.j) {
	        	sb2.append(LimitedSequence.GAP_CHARACTER);
	        }
	        else {
	        	sb2.append(s2.charAt(j - 1));
	        }
	        i = tb.i;
	        j = tb.j;
        }
         
        seqs[0] = sb1.reverse().toString();
        seqs[1] = sb2.reverse().toString();
        
        return seqs;
	}
	
	public void printAlignmentMatrix(int[][] matrix, String s1, String s2)
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
class Traceback 
{
	int k;
	int i;
	int j;
	public Traceback(int k, int i, int j) 
	{
		this.k = k;
		this.i = i;
		this.j = j;
	}

}

