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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Class that implements an FM-index to perform quick queries over large sequence databases
 * @author German Andrade
 * @author Jorge Duitama
 */
public class FMIndexSingleSequence implements Serializable 
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 5981359942407474671L;

	private static final char SPECIAL_CHARACTER = '$';
	private static final int DEFAULT_TALLY_DISTANCE = 100;
	private static final int DEFAULT_SUFFIX_FRACTION = 50;

	//Start position in the original sequence of some rows of the BW matrix representing a partial suffix array
	private Map<Integer,Integer> partialSuffixArray = new HashMap<>();

	//Ranks in the bwt for each character in the alphabet for some of the rows in the BW matrix
	private int [][] tallyIndexes;

	//1 of each tallyDistance is saved
	private int tallyDistance;

	// 1/suffixFraction indexes are saved
	private int suffixFraction;

	//Burrows Wheeler transform
	private char [] bwt;

	//For each character tells the first time it appears in the left column of the BW matrix
	private Map<Character,Integer> firstRowsInMatrix;

	//For each character tells the last time it appears in the left column of the BW matrix
	private Map<Character,Integer> lastRowsInMatrix;

	//Inferred alphabet of the sequence ordered lexicographical 
	private String alphabet;

	public FMIndexSingleSequence(CharSequence sequence) 
	{
		this(sequence,DEFAULT_TALLY_DISTANCE,DEFAULT_SUFFIX_FRACTION);
	}

	public FMIndexSingleSequence(CharSequence sequence, int tallyDistance, int suffixFraction) {
		this.tallyDistance = tallyDistance;
		this.suffixFraction = suffixFraction;
		calculate(sequence);
	}

	public int getTallyDistance() {
		return tallyDistance;
	}
	public void setTallyDistance(int tallyDistance) {
		this.tallyDistance = tallyDistance;
	}

	private void calculate(CharSequence sequence) 
	{
		List<Integer> suffixes = buildSuffixArray(sequence);
		buildBWT(sequence, suffixes);
		buildAlphabetAndCounts(sequence,suffixes);
		buildTally();
		createPartialSuffixArray(suffixes);

	}
	private List<Integer> buildSuffixArray(CharSequence sequence) {
		ArrayList<Integer> sufixes = new ArrayList<Integer>();
		for (int i = 0; i < sequence.length(); i++) {
			sufixes.add(i);
		}
		Collections.sort(sufixes, new SuffixCharSequencePositionComparator(sequence));
		return sufixes;
	}

	private void buildBWT(CharSequence sequence, List<Integer> suffixes) 
	{
		bwt = new char [sequence.length()+1];
		bwt[0] = sequence.charAt(sequence.length()-1);
		int j=1;
		for (int i:suffixes) {
			if(i>0) {
				bwt[j] = sequence.charAt(i-1);
			}
			else {
				bwt[j]= SPECIAL_CHARACTER;	
			}
			j++;
		}
	}

	private void buildAlphabetAndCounts(CharSequence seq, List<Integer> suffixArray) {
		Map<Character,Integer> counts= new TreeMap<>();
		firstRowsInMatrix = new TreeMap<>();
		lastRowsInMatrix = new TreeMap<>();
		char lastC = SPECIAL_CHARACTER;
		StringBuilder alpB = new StringBuilder();
		firstRowsInMatrix.put(lastC, 0);
		lastRowsInMatrix.put(lastC, 0);
		//iterate last column to know alphabet and counts...
		for (int i = 0; i < suffixArray.size(); i++) 
		{
			int j = suffixArray.get(i);
			char c = seq.charAt(j);
			Integer countC = counts.get(c); 
			if(countC == null) {
				counts.put(c,1);	
			} else {
				counts.put(c, countC+1);
			}
			if(lastC != c) {
				alpB.append(c);
				firstRowsInMatrix.put(c, i+1);
				lastRowsInMatrix.put(lastC, i+1);
			}
			lastC = c;
		}
		lastRowsInMatrix.put(lastC, suffixArray.size());
		alphabet = alpB.toString();

	}
	private void buildTally() 
	{
		int [] arr= new int[alphabet.length()];
		Arrays.fill(arr, 0);
		int tallyRows = bwt.length/tallyDistance;
		if(bwt.length%tallyDistance>0)tallyRows++;
		tallyIndexes = new int[tallyRows][arr.length];
		int j=0;
		for (int i=0;i<bwt.length;i++) {
			char c = bwt[i];
			if (c != SPECIAL_CHARACTER) {
				int indexC = alphabet.indexOf(c);
				if(indexC<0) throw new RuntimeException("Character "+c+" not found in the alphabet "+alphabet);
				arr[indexC]++;
			}
			if(i%tallyDistance==0) {
				int [] copy= Arrays.copyOf(arr, arr.length);
				tallyIndexes[j] = copy;
				j++;
			}
		}
	}


	private void createPartialSuffixArray(List<Integer> suffixes) 
	{
		partialSuffixArray = new HashMap<Integer,Integer>();
		int n = suffixes.size();
		for(int i=0;i<n;i++) 
		{
			int startSeq = suffixes.get(i);
			if(startSeq%suffixFraction==0) 
			{
				partialSuffixArray.put(i+1, startSeq);
			}
		}
	}

	public List<Integer> search (String searchSequence) 
	{
		int[] range = getRange(searchSequence);

		
		return getRealIndexes(range);
	}


	private List<Integer> getRealIndexes(int[] range) 
	{
		List<Integer> startIndexes=new ArrayList<>();
		if (range ==null) return startIndexes;
		//From this point is just transform the range into the real indexes in the sequence
		for (int i = range[0]; i <= range[1]; i++) 
		{
			int row = i;
			Integer begin=partialSuffixArray.get(row);
			int steps;
			for(steps =0; begin == null;steps++) 
			{
				row = lfMapping(row);
				begin=partialSuffixArray.get(row);
			}
			begin += steps;
			startIndexes.add(begin);
		}
		return startIndexes;
	}

	public int[] getRange(String query)
	{
		char actualChar = query.charAt(query.length()-1);

		Integer rowS=firstRowsInMatrix.get(actualChar);
		Integer rowF=lastRowsInMatrix.get(actualChar);
		if(rowS == null || rowF==null || rowS == -1 || rowF==-1 ) 
		{
			return null;
		}
		for(int j=query.length()-2;j>=0;j--) 
		{
			actualChar = query.charAt(j);
			if(alphabet.indexOf(actualChar)<0) return null;
			boolean add1 = (bwt[rowS]!=actualChar); 
			rowS = lfMapping(actualChar, rowS);
			if(add1) rowS++; 
			rowF = lfMapping(actualChar, rowF);
			if(rowS>rowF) 
			{
				return null;
			}
		}

		return new int[] {rowS, rowF };
	}


	public int getTallyOf(char c, int row)
	{
		int r = 0;

		int a = row/tallyDistance;
		int b = a+1;

		if( row-a*tallyDistance < b*tallyDistance-row || tallyIndexes.length<=b) {
			//Recalculate from top record
			r = tallyIndexes[a][alphabet.indexOf(c)];

			for (int j = a*tallyDistance+1; j <= row; j++) {
				char cA = bwt[j];
				if(cA==c) r++;
			}
		} else {
			//Recalculate from bottom record
			r = tallyIndexes[b][alphabet.indexOf(c)];
			for (int j = b*tallyDistance; j > row; j--) 
			{
				char cA = bwt[j];
				if(cA==c) r--;
			}
		}
		return r;
	}

	private int lfMapping(char c, int row) 
	{
		int rank = getTallyOf(c, row);
		return firstRowsInMatrix.get(c) + rank - 1;
	}

	private int lfMapping (int row)
	{
		char c = bwt[row];
		//		System.out.println(""+c);
		return lfMapping(c,row);
	}

	/*
	 * Methods for inexact matching
	 */
	/**
	 * Calculates the number of differences between w and x 
	 * @param w the string we are going to search
	 * @return
	 */
	int [] calculateD(String w)
	{
		int [] d= new int[w.length()];
		int z=0;
		int j=0;
		for (int i = 1; i <= d.length; i++) 
		{
			String sub = w.substring(j, i);
			if(search(sub).size()==0)
			{
				z++;
				j++;
			}
			d[i-1]=z;
		}
		return d;
	}



	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args)
	{
		CharSequence c = "abaaba";
		FMIndexSingleSequence f = new FMIndexSingleSequence(c);
		String query ="abc";
		//		System.out.println(query.substring(0, query.length()));

		System.out.println(f.search(query));
//		int [] d= f.calculateD(query);
//		for (int i = 0; i < d.length; i++) 
//		{
//			System.out.print(d[(i)]+"\t");
//		}
		
		/*
		int[]a= {1,2};
		@SuppressWarnings("unchecked")
		Set<int[]> set = new TreeSet(new Comparator<int[]>() 
		{
	        @Override
	        public int compare(int[] o1, int[] o2) {
	            String s1=o1[0]+";"+o2[1];
	            String s2=o2[0]+";"+o2[1];
	            return s1.compareTo(s2);
	        }
	    });
		set.add(a);
		
		int[]a1= {1,2};
		Set<int[]> set2 = new HashSet<int[]>();
		set2.add(a1);
		
		Set<int[]> set3 = new HashSet<>();
		set3.add(a);
		
		set.addAll(set2);
		set.addAll(set3);
		Iterator<int[]> it = set.iterator();
		while (it.hasNext()) {
			int[] object = (int[]) it.next();
			System.out.print("["+object[0]+" "+object[1]+"]\t");
		}
		 */

		f.searchInexact(query);
	}

	@SuppressWarnings("unchecked")
	private void searchInexact(String query) 
	{
		int[] d = calculateD(query);
		int z=1;
		Set<int[]> a= inexRecur(query,query.length()-1,z,1,bwt.length-1,d);
		System.out.println("found");
		Iterator i =a.iterator();
		while (i.hasNext()) {
			int[] object = (int[]) i.next();
			System.out.print("["+object[0]+" "+object[1]+"]\t");
			System.out.println(getRealIndexes(object));
		}
		
		
	}

	private Set<int[]> inexRecur(String w, int i, int z, int k, int l,int[]d) 
	{
		Set<int[]> arr = new TreeSet(new Comparator<int[]>() 
		{
	        @Override
	        public int compare(int[] o1, int[] o2) {
	            String s1=o1[0]+";"+o2[1];
	            String s2=o2[0]+";"+o2[1];
	            return s1.compareTo(s2);
	        }
	    });
		if(i<0)
		{
			int[] range = {k,l};
			arr.add(range);
			return arr;
		}
		if(z<d[i])
		{
			return arr;
		}
//		ArrayList
		arr.addAll(inexRecur(w, i-1, z-1, k, l, d));
		for (int j =0 ; j < alphabet.length(); j++) 
		{
			char b = alphabet.charAt(j);
			int[] rangeActual=getRange(b+"");
			k=rangeActual[0];
			l=rangeActual[1];
			if(k<=l)
			{
				arr.addAll(inexRecur(w, i, z-1, k, l, d));
				if(b==w.charAt(i))
				{
					arr.addAll(inexRecur(w, i-1, z, k, l, d));
				}
				else
				{
					arr.addAll(inexRecur(w, i-1, z-1, k, l, d));
				}
			}
		}

		return arr;
	}
}