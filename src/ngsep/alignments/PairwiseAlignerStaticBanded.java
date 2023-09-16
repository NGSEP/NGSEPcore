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

import java.util.HashMap;
import java.util.Map;

import ngsep.sequences.LimitedSequence; 

/**
 * 
 * @author Ana Sofia Castellanos
 *
 */
public class PairwiseAlignerStaticBanded implements PairwiseAligner {

	private int match =1; 

	private int mismatch = 1; 

	private int  indel = 2; 

	private int band = 20;
	
	public PairwiseAlignerStaticBanded () {}
	public PairwiseAlignerStaticBanded (int band) {
		this.band = band;
	}
   

	public int getBand() {
		return band;
	}
	public void setBand(int band) {
		this.band = band;
	}
	@Override
	public String[] calculateAlignment(CharSequence sequence1, CharSequence sequence2 ) {
		Map<Integer,Integer> dp = calculateHashMap(sequence1, sequence2);
		//return null;
		return alignSequences(dp, sequence1, sequence2);
    }
	
	private Map<Integer,Integer> calculateHashMap(CharSequence sequence1, CharSequence sequence2) {
		boolean debug = false;
		int n1 = sequence1.length();
		int n2 = sequence2.length();
		//int alignmentBand = band;
		int alignmentBand = Math.max(band, 5*Math.abs(n1-n2));
		if(debug) System.out.println("N1: "+n1+" N2: "+n2+" band: "+alignmentBand);
		Map<Integer,Integer> dp = new HashMap<>(n1*band);
		for (int row = 0; row <=n1; row++ ) {
			int firstCol = Math.max(0,row-alignmentBand);
			int lastCol = Math.min(row+alignmentBand, sequence2.length());
			for (int col = firstCol ; col<=lastCol; col++) {
				if (row == 0 && col == 0 ){
					dp.put(getHash(n1, 0, 0), 0);
				} else if (row == 0){
					dp.put(getHash(n1, row,col), dp.get(getHash(n1, row, col-1))+ insertionCost(sequence2.charAt(col-1)));
					if(debug && col==2) System.out.println("hash first row: "+getHash(n1, row, col)+" value: "+dp.get(getHash(n1, row, col)));
				} else if (col == 0){
					dp.put(getHash(n1, row,col), dp.get(getHash(n1, row-1, col))+deletionCost(sequence1.charAt(row-1)));
				} else {
					Integer diagonal = dp.get(getHash(n1, row-1, col-1));
					Integer up = dp.get(getHash(n1, row-1, col));
					Integer left = dp.get(getHash(n1, row, col-1));
					char c1 = sequence1.charAt(row-1);
					char c2 = sequence2.charAt(col-1);
					if(debug && diagonal == null && up==null && left == null) System.out.println("No previous data calculated");
					//if(debug && diagonal == null) System.out.println("Diagonal not calculated");
					//if(debug && up == null) System.out.println("Up not calculated");
					//if(debug && left == null) System.out.println("Left not calculated");
					int maxNum = -n1*n2;
					if(diagonal!=null) maxNum = diagonal+matchMismatchCost(c1, c2); 
					if(up!=null) maxNum = Math.max(maxNum, up+deletionCost(c1)); 
					if(left!=null) maxNum = Math.max(maxNum, left+insertionCost(c2));
					if(debug && getHash(n1,row,col) == 762) System.out.println("row: "+row+" col: "+col);
					dp.put(getHash(n1,row,col),maxNum);
					if(debug && row==1 && col==3) System.out.println("hash: "+getHash(n1, row, col)+" score: "+maxNum+" diag: "+diagonal+" left: "+left+" up: "+up+" c1: "+c1+" c2: "+c2);
                }
            }
        }
		if(debug) System.out.println("Final score: "+dp.get(getHash(n1, n1, n2)));
		return dp;
    }
	private String[] alignSequences(Map<Integer,Integer> dp, CharSequence sequence1, CharSequence sequence2) {
		int n1 = sequence1.length();
		int n2 = sequence2.length();
		int i = n1;
		int j = n2;
		
		StringBuffer ns1 = new StringBuffer(n1);
		StringBuffer ns2 = new StringBuffer(n2);
		while (j>0 || i>0){
			Integer actual = dp.get(getHash(n1, i,j));
			if(actual == null) {
				if(n1<n2) {
					ns1.append(LimitedSequence.GAP_CHARACTER); 
					ns2.append(sequence2.charAt(j-1)); 
					j--;
				} else {
					ns1.append(sequence1.charAt(i-1));
					ns2.append(LimitedSequence.GAP_CHARACTER);
					i--;
				}
			} else if (i==0 && j>0){
				ns1.append(LimitedSequence.GAP_CHARACTER); 
				ns2.append(sequence2.charAt(j-1)); 
				j--;
			} else if(j==0 && i>0 ){
				ns1.append(sequence1.charAt(i-1));
				ns2.append(LimitedSequence.GAP_CHARACTER);
				i--;
			} else {
				Integer diagonal = dp.get(getHash(n1, i-1, j-1));
				Integer up = dp.get(getHash(n1, i-1, j));
				Integer left = dp.get(getHash(n1, i, j-1));
				char c1 = sequence1.charAt(i-1);
				char c2 = sequence2.charAt(j-1);
				if(diagonal!=null && actual == diagonal+matchMismatchCost(c1,c2)) {
					ns1.append(c1);
					ns2.append(c2);
					j--; 
					i--;
				} else if (up!=null && actual == up+deletionCost(c1)) {
					ns1.append(sequence1.charAt(i-1)); 
					ns2.append(LimitedSequence.GAP_CHARACTER);
					i--;
				} else if (left!=null && actual == left+insertionCost(c1)) {
					ns1.append(LimitedSequence.GAP_CHARACTER);
					ns2.append(sequence2.charAt(j-1));
					j--;
				} else {
					throw new RuntimeException("Reached position without calculation. "+i+" "+j+" score: "+actual+" diag: "+diagonal+" left: "+left+" up: "+up+" c1: "+c1+" c2: "+c2);
				}
			}
		}
		String[] seqs = new String[2]; 
		seqs[0] = ns1.reverse().toString();
		seqs[1] = ns2.reverse().toString();
		return seqs;
	}

	private int insertionCost(char l1) {
		return -indel; 
	}
	private int deletionCost(char l1) {
		return -indel;
	}
	private int matchMismatchCost(char l1, char l2){
		if (l1 == l2){
			return match;
		} else {
			return -mismatch;
		}
	}
	public int getHash (int l1, int i, int j) {
		//return (""+i+"-"+j).hashCode();
		return j*(l1+1)+i;
	}    
}