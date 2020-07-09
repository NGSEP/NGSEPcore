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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.HashMap;
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
	
	//Maximum hits to return per query
	private int maxHitsQuery = 100000;
	
	

	// Inferred alphabet of the sequence ordered lexicographical
	private String alphabet;
	
	private Map<Character, Integer> alphabetIndexes;

	//Used for loading
	private FMIndexSingleSequence () {
		characterCounts = new HashMap<Character, Integer>();
		firstRowsInMatrix = new HashMap<Character, Integer>();
		lastRowsInMatrix = new HashMap<Character, Integer>();
		alphabetIndexes = new HashMap<Character, Integer>();
		
	}
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
	
	public int getMaxHitsQuery() {
		return maxHitsQuery;
	}

	public void setMaxHitsQuery(int maxHitsQuery) {
		this.maxHitsQuery = maxHitsQuery;
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
		alphabetIndexes = new HashMap<>();
		for(int i=0;i<alphabet.length();i++) alphabetIndexes.put(alphabet.charAt(i), i);
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
		if(query.length()==0) return null;
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
		for (int i = firstRow; i <= lastRow && startIndexes.size()<maxHitsQuery; i++) {
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
	
	public void save (PrintStream out) {
		out.println("#INDEX\t"+alphabet+"\t"+suffixFraction+"\t"+tallyDistance+"\t"+rowBWTSpecialCharacter+"\t"+maxHitsQuery+"\t"+bwt.length);
		for (int i=0;i<alphabet.length();i++) {
			char c = alphabet.charAt(i);
			out.println(""+c+"\t"+characterCounts.get(c)+"\t"+firstRowsInMatrix.get(c)+"\t"+lastRowsInMatrix.get(c)+"\t"+alphabetIndexes.get(c));
		}
		out.println("#PartialSuffixArray");
		for(int key:partialSuffixArray.keySet()) {
			int value = partialSuffixArray.get(key);
			out.println(""+key+"\t"+value);
		}
		out.println("#BWT");
		StringBuffer buffer = new StringBuffer(10000);
		int i=0;
		while(i<bwt.length) {
			buffer.append((char)bwt[i]);
			i++;
			if(i%10000==0 || i==bwt.length) {
				out.println(buffer.toString());
				if(i<bwt.length) buffer = new StringBuffer(10000);
			}
		}
		out.println("#END");
	}
	public static FMIndexSingleSequence load (BufferedReader reader) throws IOException {
		String line = reader.readLine();
		if(line == null) return null;
		if(!line.startsWith("#INDEX")) throw new IOException("#INDEX header not found. Line: "+line);
		String [] items = line.split("\t");
		FMIndexSingleSequence index = new FMIndexSingleSequence();
		index.alphabet = items[1];
		index.suffixFraction = Integer.parseInt(items[2]);
		index.tallyDistance = Integer.parseInt(items[3]);
		index.rowBWTSpecialCharacter = Integer.parseInt(items[4]);
		index.maxHitsQuery = Integer.parseInt(items[5]);
		int bwtLength = Integer.parseInt(items[6]);
		
		for (int i=0;i<index.alphabet.length();i++) {
			char c = index.alphabet.charAt(i);
			line = reader.readLine();
			if(line==null) throw new IOException("Unexpected end of file reading character counts.");
			items = line.split("\t");
			if(items[0].length()!=1 || c!=items[0].charAt(0)) throw new IOException("Inconsistency found reading line for character "+c+". Line: "+line);
			index.characterCounts.put(c, Integer.parseInt(items[1]));
			index.firstRowsInMatrix.put(c, Integer.parseInt(items[2]));
			index.lastRowsInMatrix.put(c, Integer.parseInt(items[3]));
			index.alphabetIndexes.put(c, Integer.parseInt(items[4]));
		}
		line = reader.readLine();
		if(line==null) throw new IOException("Unexpected end of file reading suffix array.");
		if(!line.startsWith("#PartialSuffixArray")) throw new IOException("#PartialSuffixArray section not found. Line: "+line);
		line = reader.readLine();
		while (line!=null && !line.equals("#BWT")) {
			items = line.split("\t");
			index.partialSuffixArray.put(Integer.parseInt(items[0]), Integer.parseInt(items[1]));
			line = reader.readLine();
		}
		if(line == null) throw new IOException("Unexpected end of file reading suffix array.");
		index.bwt = new byte[bwtLength];
		line = reader.readLine();
		int i=0;
		while (line!=null && !line.equals("#END")) {
			for(int j=0;j<line.length();j++) {
				if(i>=bwtLength)  throw new IOException("Inconsistent bwt length: "+bwtLength);
				index.bwt[i] = (byte) line.charAt(j);
				i++;
			}
			line = reader.readLine();
		}
		if(line == null) throw new IOException("Unexpected end of file reading bwt.");
		index.buildTally();
		return index;
		
	}
}