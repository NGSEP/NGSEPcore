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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Class that implements an FM-index to perform quick queries over large
 * sequence databases
 * 
 * @author German Andrade
 * @author Jorge Duitama
 */
public class FMIndexSingleSequence implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 5981359942407474671L;

	private static final char SPECIAL_CHARACTER = '$';
	private static final int DEFAULT_TALLY_DISTANCE = 100;
	private static final int DEFAULT_SUFFIX_FRACTION = 50;

	// Start position in the original sequence of some rows of the BW matrix
	// representing a partial suffix array
	private Map<Integer, Integer> partialSuffixArray = new HashMap<>();

	// Ranks in the bwt for each character in the alphabet for some of the rows in
	// the BW matrix
	private int[][] tallyIndexes;

	// 1 of each tallyDistance is saved
	private int tallyDistance;

	// 1/suffixFraction indexes are saved
	private int suffixFraction;

	// Burrows Wheeler transform
	private char[] bwt;

	// For each character tells the first time it appears in the left column of the
	// BW matrix
	private Map<Character, Integer> firstRowsInMatrix;

	// For each character tells the last time it appears in the left column of the
	// BW matrix
	private Map<Character, Integer> lastRowsInMatrix;

	// Inferred alphabet of the sequence ordered lexicographical
	private String alphabet;

	private int maxDifferencesInexactSearch = 1;

	public FMIndexSingleSequence(CharSequence sequence) {
		this(sequence, DEFAULT_TALLY_DISTANCE, DEFAULT_SUFFIX_FRACTION);
	}

	public FMIndexSingleSequence(CharSequence sequence, int tallyDistance, int suffixFraction) {
		this.tallyDistance = tallyDistance;
		this.suffixFraction = suffixFraction;
		calculate2(sequence);
	}

	public int getTallyDistance() {
		return tallyDistance;
	}

	public void setTallyDistance(int tallyDistance) {
		this.tallyDistance = tallyDistance;
	}

	private void calculate(CharSequence sequence) {
		//List<Integer> suffixes = buildSuffixArray(sequence);
		int [] suffixes = buildSuffixArrayDC3(sequence);
		buildBWT(sequence, suffixes);
		buildAlphabetAndCounts(sequence, suffixes);
		buildTally();
		createPartialSuffixArray(suffixes);

	}
	
	private void calculate2(CharSequence sequence) {
		//List<Integer> suffixes = buildSuffixArray(sequence);
		SuffixArrayGenerator suffixArrayGenerator=new SuffixArrayGenerator(sequence);
		alphabet = suffixArrayGenerator.getAlphabet();
		firstRowsInMatrix = suffixArrayGenerator.getFirstRowsInMatrix();
		lastRowsInMatrix = suffixArrayGenerator.getLastRowsInMatrix();
		partialSuffixArray = suffixArrayGenerator.getPartialSuffixArray(suffixFraction);
		bwt = suffixArrayGenerator.getBWT();
		tallyIndexes = suffixArrayGenerator.getTallyIndexes(tallyDistance, bwt);
	}

	/**
	 * Builds a suffix array using the default merge sort algorithm
	 * @param sequence to build the array
	 * @return int [] indexes of the sorted suffixes
	 */
	public static int [] buildSuffixArrayMergeSort(CharSequence sequence) {
		ArrayList<Integer> sufixes = new ArrayList<Integer>();
		for (int i = 0; i < sequence.length(); i++) {
			sufixes.add(i);
		}
		Collections.sort(sufixes, new SuffixCharSequencePositionComparator(sequence));
		int [] answer = new int [sufixes.size()];
		for(int i=0;i<answer.length;i++) {
			answer[i] = sufixes.get(i);
		}
		return answer;
	}
	/**
	 * Builds a suffix array using the DC3 linear algorithm for sorting
	 * @param sequence to build the array
	 * @return int [] indexes of the sorted suffixes
	 */
	public static int [] buildSuffixArrayDC3(CharSequence sequence) {
		SuffixArrayGenerator suffixArrayGenerator = new SuffixArrayGenerator(sequence);
		return suffixArrayGenerator.getSA();
	}

	private void buildBWT(CharSequence sequence, int [] suffixes) {
		bwt = new char[sequence.length() + 1];
		bwt[0] = sequence.charAt(sequence.length() - 1);
		int j = 1;
		for (int i : suffixes) {
			if (i > 0) {
				bwt[j] = sequence.charAt(i - 1);
			} else {
				bwt[j] = SPECIAL_CHARACTER;
			}
			j++;
		}
	}

	private void buildAlphabetAndCounts(CharSequence seq, int [] suffixArray) {
		Map<Character, Integer> counts = new TreeMap<>();
		firstRowsInMatrix = new TreeMap<>();
		lastRowsInMatrix = new TreeMap<>();
		char lastC = SPECIAL_CHARACTER;
		StringBuilder alpB = new StringBuilder();
		firstRowsInMatrix.put(lastC, 0);
		lastRowsInMatrix.put(lastC, 0);
		// iterate last column to know alphabet and counts...
		for (int i = 0; i < suffixArray.length; i++) {
			int j = suffixArray[i];
			char c = seq.charAt(j);
			Integer countC = counts.get(c);
			if (countC == null) {
				counts.put(c, 1);
			} else {
				counts.put(c, countC + 1);
			}
			if (lastC != c) {
				alpB.append(c);
				firstRowsInMatrix.put(c, i + 1);
				lastRowsInMatrix.put(lastC, i);
				// System.out.println("Last row "+lastC+": "+i);
				// System.out.println("First row "+c+": "+(i+1));
			}
			lastC = c;
		}
		lastRowsInMatrix.put(lastC, suffixArray.length);
		// System.out.println("Last row "+lastC+": "+suffixArray.size());
		alphabet = alpB.toString();

	}

	private void buildTally() {
		int[] arr = new int[alphabet.length()];
		Arrays.fill(arr, 0);
		int tallyRows = bwt.length / tallyDistance;
		if (bwt.length % tallyDistance > 0)
			tallyRows++;
		tallyIndexes = new int[tallyRows][arr.length];
		int j = 0;
		for (int i = 0; i < bwt.length; i++) {
			char c = bwt[i];
			if (c != SPECIAL_CHARACTER) {
				int indexC = alphabet.indexOf(c);
				if (indexC < 0)
					throw new RuntimeException("Character " + c + " not found in the alphabet " + alphabet);
				arr[indexC]++;
			}
			if (i % tallyDistance == 0) {
				int[] copy = Arrays.copyOf(arr, arr.length);
				tallyIndexes[j] = copy;
				j++;
			}
		}
	}

	private void createPartialSuffixArray(int [] suffixes) {
		partialSuffixArray = new HashMap<Integer, Integer>();
		int n = suffixes.length;
		for (int i = 0; i < n; i++) {
			int startSeq = suffixes[i];
			if (startSeq % suffixFraction == 0) {
				partialSuffixArray.put(i + 1, startSeq);
			}
		}
	}

	public Set<Integer> search(String searchSequence) {
		return exactSearch(searchSequence);
		// return inexactSearchBWAAlgorithm(searchSequence);
	}

	public Set<Integer> exactSearch(String searchSequence) {
		int[] range = getRange(searchSequence);
		// if(range!=null) System.out.println("Search sequence: "+searchSequence+"
		// range: "+range[0]+"-"+range[1]);
		// else System.out.println("No hits for search sequence: "+searchSequence);

		return getRealIndexes(range);
	}

	public Set<Integer> inexactSearchBWAAlgorithm(String searchSequence) {
		int[] d = calculateD(searchSequence);

		List<int[]> ranges = inexactRecurrentSearch(searchSequence, searchSequence.length() - 1,
				maxDifferencesInexactSearch, 1, bwt.length - 1, d);
		Set<Integer> indexes = new TreeSet<>();
		for (int[] range : ranges) {
			indexes.addAll(getRealIndexes(range));
		}
		return indexes;
	}

	public static void main(String[] args) {
		System.out.println("Testing inexactSearch BWA");
		FMIndexSingleSequence f = new FMIndexSingleSequence(args[0]);

		String query = args[1];
		Set<Integer> set = f.inexactSearchBWAAlgorithm(query);
		Iterator<Integer> i = set.iterator();
		while (i.hasNext()) {
			System.out.println(i.next());
		}
	}

	/**
	 * 
	 * @param range
	 *            of rows in the FM index
	 * @return Set<Integer> Start positions in the subject sequence (values of the
	 *         suffix array)
	 */
	private Set<Integer> getRealIndexes(int[] range) {
		Set<Integer> startIndexes = new TreeSet<>();
		if (range == null)
			return startIndexes;
		// From this point is just transform the range into the real indexes in the
		// sequence
		for (int i = range[0]; i <= range[1]; i++) {
			int row = i;
			Integer begin = partialSuffixArray.get(row);
			int steps;
			for (steps = 0; begin == null; steps++) {
				row = lfMapping(row);
				begin = partialSuffixArray.get(row);
			}
			begin += steps;
			startIndexes.add(begin);
		}
		return startIndexes;
	}

	public int[] getRange(String query) {
		char actualChar = query.charAt(query.length() - 1);

		Integer rowS = firstRowsInMatrix.get(actualChar);
		Integer rowF = lastRowsInMatrix.get(actualChar);
		// System.out.println("Char: "+actualChar+" Range: "+rowS+"-"+rowF);
		if (rowS == null || rowF == null || rowS == -1 || rowF == -1) {
			return null;
		}
		for (int j = query.length() - 2; j >= 0; j--) {
			actualChar = query.charAt(j);
			if (alphabet.indexOf(actualChar) < 0)
				return null;
			rowS = lfMapping(actualChar, rowS, true);
			rowF = lfMapping(actualChar, rowF, false);
			if (rowS > rowF) {
				return null;
			}
			// System.out.println("Char: "+actualChar+" Range: "+rowS+"-"+rowF);
		}

		return new int[] { rowS, rowF };
	}

	public int getTallyOf(char c, int row) {
		int r = 0;

		int a = row / tallyDistance;
		int b = a + 1;

		if (row - a * tallyDistance < b * tallyDistance - row || tallyIndexes.length <= b) {
			// Recalculate from top record
			r = tallyIndexes[a][alphabet.indexOf(c)];

			for (int j = a * tallyDistance + 1; j <= row; j++) {
				char cA = bwt[j];
				if (cA == c)
					r++;
			}
		} else {
			// Recalculate from bottom record
			r = tallyIndexes[b][alphabet.indexOf(c)];
			for (int j = b * tallyDistance; j > row; j--) {
				char cA = bwt[j];
				if (cA == c)
					r--;
			}
		}
		return r;
	}

	/**
	 * Finds the row corresponding to the given character in the given row of the
	 * index, according to the tally indexes in that row
	 * 
	 * @param c
	 *            Character to query
	 * @param row
	 *            of the index to query
	 * @param firstIndexAfter
	 *            If true, calculates the rank of the character at or after the row
	 * @return int Row of the FM-index of the rank of the given character according
	 *         to the tally indexes at the given row
	 */
	private int lfMapping(char c, int row, boolean firstIndexAfter) {

		int rank = getTallyOf(c, row);
		// add1 is true when actualChar is different of bwt[rowS] because in this case,
		// the last appearance of actualChar before rowS is outside the range defined by
		// rowS, rowF
		boolean add1 = firstIndexAfter && (bwt[row] != c);
		// System.out.println("char: "+c+" row: "+row+" rank: "+rank+" first c:
		// "+firstRowsInMatrix.get(c));
		int newRank = firstRowsInMatrix.get(c) + rank - 1;
		if (add1)
			newRank++;
		return newRank;
	}

	private int lfMapping(int row) {
		char c = bwt[row];
		// System.out.println(""+c);
		return lfMapping(c, row, false);
	}

	/*
	 * Methods for inexact matching
	 */
	/**
	 * Calculates the number of differences between w and x
	 * 
	 * @param query
	 *            the string we are going to search
	 * @return
	 */
	int[] calculateD(String query) {
		int[] d = new int[query.length()];
		int z = 0;
		int j = 0;
		for (int i = 1; i <= d.length; i++) {
			if (j <= i) {
				String sub = query.substring(j, i);
				if (!sub.equals("")) {

					if (exactSearch(sub).size() == 0) {
						System.out.println("sub " + sub);
						z++;
						j = i + 1;
					}
				}
			}

			d[i - 1] = z;
		}
		return d;
	}

	/**
	 * Makes an inexact search over the index
	 * 
	 * @param query
	 * @param lastIdx
	 *            Last index of the query to search
	 * @param maxDiff
	 *            Maximum number of differences between the query and
	 * @param firstRow
	 *            First row of the FM index to look for
	 * @param lastRow
	 *            Last row of the FM index to look for
	 * @param d
	 * @return
	 */
	private List<int[]> inexactRecurrentSearch(String query, int lastIdxQuery, int maxDiff, int firstRow, int lastRow,
			int[] d) {
		// System.out.println("Recursi�n lastIdxQuery:"+lastIdxQuery+"
		// maxDiff:"+maxDiff+" firstRow:"+firstRow+" lastRow:"+lastRow);
		List<int[]> arr = new ArrayList<>();
		if (maxDiff < 0) {
			// System.out.println("sale");
			// Base case when the number of differences is larger than the maximum allowed
			return arr;
		}
		if (lastIdxQuery < 0) {
			// Base case an empty query
			// System.out.println("\tagrega {"+firstRow+","+lastRow+"} en
			// maxDiff:"+maxDiff);
			int[] range = { firstRow, lastRow };
			arr.add(range);
			return arr;
		}
		if (maxDiff < d[lastIdxQuery]) {
			// Base case when the number of differences is larger than the maximum allowed
			return arr;
		}
		// Insertion
		arr.addAll(inexactRecurrentSearch(query, lastIdxQuery - 1, maxDiff - 1, firstRow, lastRow, d));
		for (int j = 0; j < alphabet.length(); j++) {
			char b = alphabet.charAt(j);
			int temp = 0;
			if (firstRow != 1) {
				temp = firstRow;
			}
			int newFirst = lfMapping(b, temp, true);
			int newLast = lfMapping(b, lastRow, false);
			if (newFirst <= newLast) {
				// Deletion
				arr.addAll(inexactRecurrentSearch(query, lastIdxQuery, maxDiff - 1, newFirst, newLast, d));
				if (b == query.charAt(lastIdxQuery)) {
					// System.out.println("entra "+b);
					// Follow indexes for the matching character
					arr.addAll(inexactRecurrentSearch(query, lastIdxQuery - 1, maxDiff, newFirst, newLast, d));
				} else {
					// Follow indexes for the possible mismatch characters reducing in 1 the number
					// of differences
					arr.addAll(inexactRecurrentSearch(query, lastIdxQuery - 1, maxDiff - 1, newFirst, newLast, d));
				}
			}
		}

		return arr;
	}
	
	public String reverseBWT()
	{
		String r = "$";
		int rowi=0;
		while(bwt[rowi]!='$')
		{
			char c = bwt[rowi];
			r=c+r;
			rowi=lfMapping(c,  rowi, false);
		}
		return r;
	}

	public String getSequenceSubString(int first, int last) 
	{
		// TODO Auto-generated method stub 
		
		return reverseBWT().substring(first, last);
	}
}