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
package ngsep.alignments;

import ngsep.sequences.LimitedSequence;

/**
 * Performs pairwise alignment using the simple gap scoring method.
 * Adapted from https://www.itu.dk/~sestoft/bsa/Match2.java
 * @author David Guevara
 */
public class PairwiseAlignerSimpleGap implements PairwiseAligner {
	
	private int match=1;
	private int openGap=2;
	private int mismatch=1;
	
	private boolean forceStart1 = true;
	private boolean forceStart2 = true;
	private boolean forceEnd1 = true;
	private boolean forceEnd2 = true;
	private boolean local = false;
	
	private int[][] matchScores;
	private int maxScore = 0;
	
	public PairwiseAlignerSimpleGap() {
		
	}
	public PairwiseAlignerSimpleGap(int capacity) 
	{
		matchScores = new int [capacity][capacity];
	}
	
	public int getMatch() {
		return match;
	}

	public void setMatch(int match) {
		this.match = match;
	}

	public int getOpenGap() {
		return openGap;
	}

	public void setOpenGap(int openGap) {
		this.openGap = openGap;
	}

	public int getMismatch() {
		return mismatch;
	}

	public void setMismatch(int mismatch) {
		this.mismatch = mismatch;
	}

	public boolean isForceStart1() {
		return forceStart1;
	}

	public void setForceStart1(boolean forceStart1) {
		this.forceStart1 = forceStart1;
	}

	public boolean isForceStart2() {
		return forceStart2;
	}

	public void setForceStart2(boolean forceStart2) {
		this.forceStart2 = forceStart2;
	}

	public boolean isForceEnd1() {
		return forceEnd1;
	}

	public void setForceEnd1(boolean forceEnd1) {
		this.forceEnd1 = forceEnd1;
	}

	public boolean isForceEnd2() {
		return forceEnd2;
	}

	public void setForceEnd2(boolean forceEnd2) {
		this.forceEnd2 = forceEnd2;
	}
	
	

	public boolean isLocal() {
		return local;
	}
	public void setLocal(boolean local) {
		this.local = local;
		forceStart1 = forceEnd1 = forceStart2 = forceEnd2 = false;
	}
	
	
	
	public int[][] getMatchScores() {
		return matchScores;
	}

	public int getMaxScore() {
		return maxScore;
	}
	
	public String[] calculateAlignment(CharSequence s1, CharSequence s2) 
	{		
		initMatrices(s1, s2);
	    calculateMatrices(s1, s2);	    
        return getAlignedStrings(s1, s2);
	}
	
	private void initMatrices(CharSequence s1, CharSequence s2)
	{
		if(matchScores==null || matchScores.length<s1.length()+1 || matchScores[0].length < s2.length() +1 ) {
			matchScores = new int[s1.length() + 1][s2.length() + 1];
		}
		
		matchScores[0][0] = 0;
		for (int i = 1; i < matchScores.length; i++) 
		{
			if (forceStart1) matchScores[i][0] = - openGap * i;
			else matchScores[i][0] = 0;
	    }
	    for (int i = 1; i < matchScores[0].length; i++) 
	    {
	    	if (forceStart2) matchScores[0][i] = - openGap * i;
	    	else matchScores[0][i] = 0;
	    }
	}
	
	private void calculateMatrices(CharSequence s1, CharSequence s2)
	{
		maxScore = 0;
		for (int i = 1; i <= s1.length(); i++)
	    {
	    	for (int j = 1; j <= s2.length(); j++)
	    	{
	    		int matchScore = getMatchScore(s1.charAt(i - 1), s2.charAt(j - 1));
	    		matchScores[i][j] = Math.max(matchScores[i-1][j-1] + matchScore, Math.max(matchScores[i-1][j] - openGap, matchScores[i][j-1] - openGap));
	    		if(local) matchScores[i][j] = Math.max(matchScores[i][j], 0);
	    		maxScore = Math.max(maxScore, matchScores[i][j]);
	    	}
	    }
	}

	private int getMatchScore(char a, char b)
	{
		if (a == b)
			return match;
		else
			return -mismatch;
	}
	
	
	
	private String[] getAlignedStrings(CharSequence s1, CharSequence s2)
	{
		StringBuffer sb1 = new StringBuffer();
		StringBuffer sb2 = new StringBuffer();
		int i = s1.length();
		int j = s2.length();
	    int val = matchScores[i][j];
	    if (local) {
	    	int maxI=0;
	    	int maxJ=0;
	    	for (i = 1; i <= s1.length(); i++)
		    {
		    	for (j = 1; j <= s2.length(); j++)
		    	{
		    		if(maxScore==matchScores[i][j]) {
		    			maxI = i;
		    			maxJ = j;
		    		}
		    	}
		    }
	    	i=maxI;
	    	j=maxJ;
	    }
	    else if (!forceEnd1) {
    		// Find better score over the last column
    		for (int h=i;h>=0;h--) {
    			int score = matchScores[h][s2.length()];
    			if (score>val) {
    				i=h;
    				val = score; 
    			}
    		}
    	}
	    else if (!forceEnd2) {
    		// Find better score over the last row
    		for (int h=j;h>=0;h--) {
    			int score = matchScores[s1.length()][h];
    			if (score>val) {
    				i=s1.length();
    				j=h;
    				val = score; 
    			}
    		}
    	}
	    if(!local) {
	    	for (int h = s1.length();h>i;h--) {
	    		sb1.append(s1.charAt(h - 1));
				sb2.append(LimitedSequence.GAP_CHARACTER);
	    	}
	    	for (int h = s2.length();h>j;h--) {
	    		sb1.append(LimitedSequence.GAP_CHARACTER);
				sb2.append(s2.charAt(j - 1));
	    	}
	    }
    	
    	
    	// Traceback cycle
		while(i>0 && j>0) {
			int matchScore = getMatchScore(s1.charAt(i - 1), s2.charAt(j - 1));
			int score = matchScores[i][j];
			if(local && score==0) break;
			if(score == matchScores[i-1][j-1] + matchScore) {
				sb1.append(s1.charAt(i - 1));
				sb2.append(s2.charAt(j - 1));
				i--;
				j--;
			} else if(score == matchScores[i-1][j] - openGap) {
				sb1.append(s1.charAt(i - 1));
				sb2.append(LimitedSequence.GAP_CHARACTER);
				i--;
			}
    		else if(score == matchScores[i][j-1] - openGap) {
    			sb1.append(LimitedSequence.GAP_CHARACTER);
    			sb2.append(s2.charAt(j - 1));
    			j--;
    		}
    		else throw new RuntimeException("Unexpected score error at "+i+" "+j);
        }
		if(!local) {
			while (i>0) {
				sb1.append(s1.charAt(i - 1));
				sb2.append(LimitedSequence.GAP_CHARACTER);
				i--;
			}
			while (j>0) {
				sb1.append(LimitedSequence.GAP_CHARACTER);
				sb2.append(s2.charAt(j - 1));
				j--;
			}
		}
        String[] seqs = new String[2]; 
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