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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Class that implements an FM-index to perform quick queries over large
 * sequence databases
 * 
 * @author German Andrade
 * @author Jorge Duitama
 * @author Juan Camilo Bojaca
 */
public class FMIndexSingleSequence implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 5981359942407474671L;

	/** Character to BWT */
	public static final char SPECIAL_CHARACTER = 0;
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
	private byte [] bwt;
	private int rowBWTSpecialCharacter;

	//For each character tells the number of times it appears
	private Map<Character, Integer> characterCounts;
	// For each character tells the first time it appears in the left column of the
	// BW matrix
	private Map<Character, Integer> firstRowsInMatrix;

	// For each character tells the last time it appears in the left column of the
	// BW matrix
	private Map<Character, Integer> lastRowsInMatrix;
	
	

	// Inferred alphabet of the sequence ordered lexicographical
	private String alphabet;
	
	private Map<Character, Integer> alphabetIndexes;

	private int maxDifferencesInexactSearch = 1;

	public FMIndexSingleSequence(CharSequence sequence) {
		this(sequence, DEFAULT_TALLY_DISTANCE, DEFAULT_SUFFIX_FRACTION);
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
	/**
	 * @return Length of the sequence represented by this FMIndex
	 */
	public int getSequenceLength() {
		return bwt.length-1;
	}

	private void calculate(CharSequence sequence) {
		countCharacters (sequence);
		buildCharacterFirstAndLastRows();
		SuffixArrayGenerator suffixArrayGenerator = new DC3SuffixArrayGenerator(sequence);
		//SuffixArrayGenerator suffixArrayGenerator = new CollectionsSortSuffixArrayGenerator(sequence);
		int [] sa = suffixArrayGenerator.getSuffixArray();
		//System.out.println("First pos SA: "+sa[0]+" "+sa[1]+" "+sa[2] );
		
		
		alphabetIndexes = new HashMap<>();
		for(int i=0;i<alphabet.length();i++) alphabetIndexes.put(alphabet.charAt(i), i);
		
		
		buildBWT(sequence, sa);
		
		createPartialSuffixArray(sa);
		buildTally();
		//printIndexInfo();
	}

	private void countCharacters(CharSequence sequence) {
		characterCounts = new HashMap<Character, Integer>();
		for(int i=0;i<sequence.length();i++) {
			char c = sequence.charAt(i);
			characterCounts.compute(c, (k, v) -> (v == null) ? 1 : v+1);
		}
		Set<Character> sortedAlphabet = new TreeSet<Character>();
		sortedAlphabet.addAll(characterCounts.keySet());
		StringBuilder alphB = new StringBuilder();
		for(char c:sortedAlphabet) alphB.append(c);
		alphabet = alphB.toString();
	}

	private void buildCharacterFirstAndLastRows() {
		firstRowsInMatrix = new HashMap<>();
		lastRowsInMatrix = new HashMap<>();
		
		firstRowsInMatrix.put(SPECIAL_CHARACTER, 0);
		lastRowsInMatrix.put(SPECIAL_CHARACTER, 0);
		int totalChars = 1;
		for(int i=0;i<alphabet.length();i++) {
			char c = alphabet.charAt(i);
			firstRowsInMatrix.put(c, totalChars);
	    	totalChars += characterCounts.get(c);
	    	lastRowsInMatrix.put(c, totalChars - 1);
		}
	}
	
	private void printIndexInfo() {
		System.out.println("Alphabet: "+alphabet);
		System.out.println("BWT: "+new String(bwt));
		System.out.println("Partial array: "+partialSuffixArray);
		System.out.println("First rows: "+firstRowsInMatrix);
		System.out.println("Last rows: "+lastRowsInMatrix);
	}
	
	private void buildBWT(CharSequence sequence, int [] sa) {
		bwt = new byte[sequence.length() + 1];
	
		if(sa[0]!=sequence.length()) throw new RuntimeException("Suffix array should have "+sequence.length()+" as first entry");
		//assert sa[0]==sequence.length();
		int j = 0;
		for (int i : sa) {
			if (i > 0) {
				bwt[j] = (byte)sequence.charAt(i - 1);
			} else {
				bwt[j] = SPECIAL_CHARACTER;
				rowBWTSpecialCharacter = j;
			}
			j++;
		}
	}

	private void buildTally() {
		int tallyRows = bwt.length / tallyDistance;
		if (bwt.length % tallyDistance > 0) tallyRows++;
		
		final int[] arr = new int[alphabet.length()];
		tallyIndexes = new int[tallyRows][alphabet.length()];
		

		int j = 0;
		for (int i = 0; i < bwt.length; i++) {
			char c = (char)bwt[i];
			if (c != SPECIAL_CHARACTER) {
				int indexC = alphabetIndexes.get(c);
				arr[indexC]++;
			}
			if (i % tallyDistance == 0) {
				System.arraycopy(arr, 0, tallyIndexes[j], 0, arr.length);
				j++;
			}
		}
	}

	private void createPartialSuffixArray(int [] sa) {
		partialSuffixArray = new HashMap<>();
		partialSuffixArray.put(0, sa[0]);
		//Partial suffix array module should be calculated on the suffix values (real sequence positions)
		for (int i = 1; i < sa.length-1; i ++) {
			int value = sa[i];
			if(value%suffixFraction==0) partialSuffixArray.put(i, sa[i]);
		}
		partialSuffixArray.put(sa.length-1, sa[sa.length-1]);
	}

	/**
	 * Searches the given sequence in this FMIndex
	 * @param searchSequence Sequence to search
	 * @return Set<Integer> Set of start positions for the given sequence
	 */
	public Set<Integer> ungappedSearch(String searchSequence) {
		//printIndexInfo();
		return exactSearch(searchSequence);
		// return inexactSearchBWAAlgorithm(searchSequence);
	}

	public Set<Integer> exactSearch(String searchSequence) {
		int[] range = getRange(searchSequence);
		if(range == null) {
			//System.out.println("No hits for search sequence: "+searchSequence);
			return new TreeSet<>();
		}
		//System.out.println("Search sequence: "+searchSequence+"range: "+range[0]+"-"+range[1]);
		return getSequenceIndexes(range[0],range[1]);
	}
	
	/**
	 * Looks for the range of row indexes in this index having matches to the given query
	 * @param query sequence
	 * @return int [] array with two indexes, the first and last row of this index having exact matches to the given query
	 * null if the sequence can not be found
	 */
	public int[] getRange(String query) {
		char actualChar = query.charAt(query.length() - 1);

		Integer rowS = firstRowsInMatrix.get(actualChar);
		Integer rowF = lastRowsInMatrix.get(actualChar);
		//System.out.println("Char: "+actualChar+" Range: "+rowS+"-"+rowF);
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
			 //System.out.println("Char: "+actualChar+" Range: "+rowS+"-"+rowF);
		}

		return new int[] { rowS, rowF };
	}

	/**
	 * Provides the start indexes in the original sequence corresponding to the given start 
	 * @param firstRow of this index
	 * @param lastRow of this index
	 * @return Set<Integer> Start positions in the subject sequence (values of the suffix array)
	 */
	public Set<Integer> getSequenceIndexes(int firstRow, int lastRow) {
		Set<Integer> startIndexes = new TreeSet<>();
		// From this point is just transform the range into the real indexes in the
		// sequence
		for (int i = firstRow; i <= lastRow; i++) {
			int row = i;
			Integer begin = partialSuffixArray.get(row);
			int steps;
			for (steps = 0; begin == null; steps++) {
				//System.out.println("Next row: "+row+" bwt: "+((char)bwt[row])+" steps: "+steps);
				row = lfMapping(row);
				begin = partialSuffixArray.get(row);
			}
			begin += steps;
			startIndexes.add(begin);
		}
		return startIndexes;
	}

	/**
	 * Returns the tally count for the given character in the given row of this index 
	 * @param c character to count. c must belong to the alphabet of this index
	 * @param row to query
	 * @return int count of appearances of the character c in the bwt up to the given row
	 */
	public int getTallyCount(char c, int row) {
		if(c==SPECIAL_CHARACTER) {
			return (row>=rowBWTSpecialCharacter)?1:0;
		}
		int r = 0;

		int a = row / tallyDistance;
		int b = a + 1;

		if (row - a * tallyDistance < b * tallyDistance - row || tallyIndexes.length <= b) {
			// Recalculate from top record
			r = tallyIndexes[a][alphabet.indexOf(c)];

			for (int j = a * tallyDistance + 1; j <= row; j++) {
				char cA = (char)bwt[j];
				if (cA == c)
					r++;
			}
		} else {
			// Recalculate from bottom record
			r = tallyIndexes[b][alphabet.indexOf(c)];
			for (int j = b * tallyDistance; j > row; j--) {
				char cA = (char)bwt[j];
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
	 * @param c Character to query
	 * @param row of the index to query
	 * @param firstIndexAfter If true, calculates the rank of the character at or after the row
	 * @return int Row of the FM-index of the rank of the given character according
	 *         to the tally indexes at the given row
	 */
	private int lfMapping(char c, int row, boolean firstIndexAfter) {

		int rank = getTallyCount(c, row);
		// add1 is true when actualChar is different of bwt[rowS] because in this case,
		// the last appearance of actualChar before rowS is outside the range defined by
		// rowS, rowF
		boolean add1 = firstIndexAfter && (bwt[row] != c);
		// System.out.println("char: "+c+" row: "+row+" rank: "+rank+" first c: "+firstRowsInMatrix.get(c));
		int newRank = firstRowsInMatrix.get(c) + rank - 1;
		if (add1) newRank++;
		return newRank;
	}

	private int lfMapping(int row) {
		char c = (char)bwt[row];
		// System.out.println(""+c);
		return lfMapping(c, row, false);
	}
	
	public static void main(String[] args) {
		System.out.println("Testing inexactSearch BWA");
		FMIndexSingleSequence f = new FMIndexSingleSequence(args[0]);

		String query = args[1];
		f.printIndexInfo();
		f.printSA();
		
		Set<Integer> set = f.ungappedSearch(query);
		//Set<Integer> set = f.inexactSearchBWAAlgorithm(query);
		
		System.out.println("Result indexes: "+set);
	}
	
	private void printSA() {
		Set<Integer> allIndexes = getSequenceIndexes(1, 1);
		for(int i:allIndexes) {
			System.out.println(""+i);
		}
		
	}
	
	
	/*
	 * Methods for inexact matching
	 */
	


	
	
	

	public Set<Integer> inexactSearchBWAAlgorithm(String searchSequence) {
		int[] d = calculateD(searchSequence);

		List<int[]> ranges = inexactRecurrentSearch(searchSequence, searchSequence.length() - 1,
				maxDifferencesInexactSearch, 1, bwt.length - 1, d);
		Set<Integer> indexes = new TreeSet<>();
		for (int[] range : ranges) {
			indexes.addAll(getSequenceIndexes(range[0], range[1]));
		}
		return indexes;
	}
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
		// System.out.println("Recursion lastIdxQuery:"+lastIdxQuery+"
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
}