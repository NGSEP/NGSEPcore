package ngsep.sequences;

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

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * @author Juan Camilo Bojaca
 *
 */
public class SuffixArrayGenerator {
	// -------------------------------------------------------------------------------
	// CONSTANT VALUES
	// -------------------------------------------------------------------------------
	private static final char SPECIAL_CHARACTER = '$';
	/**
	 * models the radix factor (bits used to sort)
	 */
	private final static int Cte_Radix = 8;
	/**
	 * represents the maximum achievable value with CteRadix bites
	 */
	private final static int Cte_Radix_Bit = (int) Math.pow(2, Cte_Radix) - 1;
	/**
	 * contains the number of occurrences of each hexadecimal digit
	 */
	private final static int[] contSort = new int[Cte_Radix_Bit + 2];
	/**
	 * the ASCII interval that it handles
	 */
	private final static int AlfabetEnd = 127;

	private final static byte[] map = new byte[AlfabetEnd];
	// -------------------------------------------------------------------------------
	// VARAIBLES
	// -------------------------------------------------------------------------------
	/**
	 * the suffix Array
	 */
	private CharSequence sequence;
	private int[] suffixArray;
	private String alphabet;
	private Map<Character, Integer> firstRowsInMatrix;
	private Map<Character, Integer> lastRowsInMatrix;

	/**
	 * is used to sort
	 */
	private final boolean[] repeated;
	private final byte[] auxSort2;
	private final int[] auxSort, sort1, sort2, partialSA;
	private final Stack loStack = new Stack(), hiStack = new Stack(), loStackAux = new Stack(),
			hiStackAux = new Stack();

	// -------------------------------------------------------------------------------
	// CONSTRUCTOR && PUBLIC METHODS
	// -------------------------------------------------------------------------------
	/**
	 * 
	 * @param charSequence
	 *            the sequence to which the suffix array is calculated
	 */
	public SuffixArrayGenerator(CharSequence charSequence) {
		long ini = System.currentTimeMillis();

		sequence = charSequence;
		getMap(charSequence);
		int[] data = transform(charSequence);

		partialSA = new int[data.length];
		suffixArray = new int[data.length];
		repeated = new boolean[(int) ((2 * data.length) / 3)];
		auxSort = new int[(int) ((2 * data.length) / 3)];
		auxSort2 = new byte[(int) ((2 * data.length) / 3)];
		sort1 = new int[(int) ((2 * data.length) / 3)];
		sort2 = new int[((data.length - 1) / 3) + 1];

		// map[map.length - 1] is the maximum possible value in data
		getSuffix(data, round(highestOneBitPos(map[map.length - 1])));
		System.out.println("time sort: " + ((System.currentTimeMillis() - ini) / (double) 1000) + "("
				+ charSequence.length() + ")");
		for (int i = 0; i < suffixArray.length; i++)
			partialSA[suffixArray[i]] = i;
	}

	public int[] getSA() {
		int[] answer = new int[suffixArray.length - 1];
		System.arraycopy(suffixArray, 1, answer, 0, answer.length);
		return answer;
	}

	public char[] getBWT() {
		char[] bwt = new char[suffixArray.length];
		bwt[partialSA[0]] = SPECIAL_CHARACTER;
		for (int i = 1; i < suffixArray.length; i++) {
			bwt[partialSA[i]] = sequence.charAt(i - 1);
		}
		return bwt;
	}

	public Map<Integer, Integer> getPartialSuffixArray(final int suffixFraction) {
		Map<Integer, Integer> partialSuffixArray = new HashMap<Integer, Integer>();
		for (int i = 0; i < suffixArray.length; i += suffixFraction)
			partialSuffixArray.put(partialSA[i], i);
		return partialSuffixArray;
	}

	public int[][] getTallyIndexes(int tallyDistance, char[] bwt) {
		int[] arr = new int[alphabet.length() + 1];
		int tallyRows = (bwt.length + (tallyDistance - 1)) / tallyDistance;
		int[][] tallyIndexes = new int[tallyRows][alphabet.length()];

		arr[map[bwt[0]]]++;
		System.arraycopy(arr, 1, tallyIndexes[0], 0, arr.length - 1);

		for (int j = 1; j < tallyRows; j++) {
			for (int i = (j - 1) * tallyDistance + 1; i <= j * tallyDistance; i++)
				arr[map[bwt[i]]]++;
			System.arraycopy(arr, 1, tallyIndexes[j], 0, arr.length - 1);
		}
		return tallyIndexes;
	}

	public String getAlphabet() {
		return alphabet;
	}

	public Map<Character, Integer> getFirstRowsInMatrix() {
		return firstRowsInMatrix;
	}

	public Map<Character, Integer> getLastRowsInMatrix() {
		return lastRowsInMatrix;
	}

	// -------------------------------------------------------------------------------
	// PRIVATE METHODS
	// -------------------------------------------------------------------------------

	/**
	 * 
	 * @param data
	 *            the chain to get the suffix array
	 * @param MaxBit
	 *            the maximum possible multiple Cte_Radix shift for at least one
	 *            element of data to be different from 0
	 * @return the suffix Array of data
	 */
	private void getSuffix(final int[] data, final int MaxBit) {
		// -----------------------------------------------------------------------------
		// STEP #1: create C (sort1)
		// -----------------------------------------------------------------------------
		int m = moduleThree(sort1, 1, 0, data.length);
		int sort1Size = moduleThree(sort1, 2, m, data.length);
		// -----------------------------------------------------------------------------
		// STEP #2: radix sort C = b1b2
		// -----------------------------------------------------------------------------
		changeSortInterval(0, sort1Size - 1);

		for (int d = 0; !loStack.isEmpty() && d < 3; d++)
			for (int bit = MaxBit; !loStack.isEmpty() && bit >= 0; bit -= Cte_Radix)
				sort(sort1, data, d, bit);
		// -----------------------------------------------------------------------------
		// STEP #2.1: recursion
		// -----------------------------------------------------------------------------
		if (!loStack.isEmpty()) {
			calculateRepeatedIndexes(loStack, hiStack, sort1Size);
			int[] r = new int[sort1Size + 1];
			int maxValueR = calculateR(r, m, sort1Size);

			getSuffix(r, round(highestOneBitPos(maxValueR)));

			// using the Suffix array of r, finishes order sort1
			for (int j = 1; j <= sort1Size; j++) {
				int idx = suffixArray[j];
				sort1[j - 1] = (idx < m) ? idx * 3 + 1 : (idx - m) * 3 + 2;
			}
		}
		calculatePartialSA(sort1Size, data.length);
		// -----------------------------------------------------------------------------
		// STEP #3: Sorting the non sample Suffices b0
		// -----------------------------------------------------------------------------
		int sort2Size = moduleThree(sort2, 0, 0, data.length);
		changeSortInterval(0, sort2Size - 1);

		for (int bit = MaxBit; !loStack.isEmpty() && bit >= 0; bit -= Cte_Radix)
			sort(sort2, data, 0, bit);
		for (int bit = round(highestOneBitPos(sort1Size)); !loStack.isEmpty() && bit >= 0; bit -= Cte_Radix)
			sort(sort2, partialSA, 0, bit);
		// -----------------------------------------------------------------------------
		// STEP #4: Merge
		// -----------------------------------------------------------------------------
		merge(sort1Size, sort2Size, data);

	}

	/**
	 * define the interval to order
	 * 
	 * @param lo
	 *            start
	 * @param hi
	 *            end
	 */
	private void changeSortInterval(final int lo, final int hi) {
		loStack.add(lo);
		hiStack.add(hi);
	}

	/**
	 * 
	 * @param array
	 *            the array of indexes to sort
	 * @param loStack
	 *            , hiStack contains the intervals to be sorted, ends with the
	 *            intervals where there are repeated elements
	 * @param data
	 *            , d , bit every element i of array is ordered according to the
	 *            function: (data[i + d] >>> bit) & Cte_Radix_Bit
	 * 
	 *            this extracts a hexadecimal number according to the bit shift of
	 *            the d-th data element after i.
	 */
	private void sort(int[] array, final int[] data, final int d, final int bit) {
		radixCount(array, data, d, bit);
	}

	/**
	 * 
	 * @param array
	 *            the array of indexes to sort
	 * @param loStack
	 *            , hiStack contains the intervals to be sorted, ends with the
	 *            intervals where there are repeated elements
	 * @param data
	 *            , d , bit every element i of array is ordered according to the
	 *            function: (data[i + d] >>> bit) & Cte_Radix_Bit
	 * 
	 *            this extracts a hexadecimal number according to the bit shift of
	 *            the d-th data element after i.
	 */
	private void radixCount(int[] array, final int[] data, final int d, final int bit) {
		while (!loStack.isEmpty()) {
			int lo = loStack.pop();
			int hi = hiStack.pop();

			Arrays.fill(contSort, 0);

			for (int i = lo; i <= hi; i++) {
				int ind = Cte_Radix_Bit & (data[array[i] + d] >>> bit);
				auxSort2[i] = (byte) ind;
				contSort[ind]++;
			}

			contSort[0] += lo - 1;
			for (int i = 1, prev = 0, cont = contSort[0]; prev != hi; i++) {
				prev = cont;
				contSort[i] = cont += contSort[i];
			}

			for (int i = lo; i <= hi; i++) {
				auxSort[contSort[Byte.toUnsignedInt(auxSort2[i])]--] = array[i];
			}

			System.arraycopy(auxSort, lo, array, lo, (hi - lo + 1));

			for (int act, prev = contSort[0], i = 1; prev != hi; prev = act, i++) {
				act = contSort[i];
				if (prev + 1 < act) {
					loStackAux.add(prev + 1);
					hiStackAux.add(act);
				}
			}
		}
		Stack.exchange(loStack, loStackAux);
		Stack.exchange(hiStack, hiStackAux);
	}

	/**
	 * 
	 * @param sort1Size
	 *            suffix array of c
	 * @param sort2Size
	 *            suffix array of b0
	 * @param partialSA
	 *            partial suffix matrix with the c indexes
	 * @param data
	 *            the original data
	 * @return the SA
	 */
	private int[] merge(final int sort1Size, final int sort2Size, final int[] data) {
		int index1 = 0, index2 = 0, indexAns = 0;
		while (index1 != sort1Size && index2 != sort2Size) {
			if (compare(sort2[index2], sort1[index1], data) < 0)
				suffixArray[indexAns++] = sort2[index2++];
			else
				suffixArray[indexAns++] = sort1[index1++];
		}
		while (index1 != sort1Size)
			suffixArray[indexAns++] = sort1[index1++];
		while (index2 != sort2Size)
			suffixArray[indexAns++] = sort2[index2++];
		return suffixArray;
	}

	/**
	 * 
	 * @param r
	 *            for every element it contains the ranking obtained when ordering
	 *            the corresponding c(original values ​​of sort1) element
	 * @param m
	 *            in sort[m] starts the module 2 values
	 * @param sort1Size
	 *            the sorted indexes
	 * @param repeated
	 *            array that represents the repeated elements
	 * @return the maximum value in r
	 */
	private int calculateR(int[] r, final int m, final int sort1Size) {
		int d = 0;
		for (int i = 0; i < sort1Size; ++i) {
			int idx = sort1[i];
			if (!repeated[i])
				++d;
			r[idx % 3 == 1 ? idx / 3 : idx / 3 + m] = d;
		}
		return d;
	}

	/**
	 * 
	 * change the format of repeated intervals to a Boolean array
	 * 
	 * @param loStack,
	 *            hiStack repeated intervals
	 * @param size
	 *            of the array
	 * @return repeated intervals in Boolean array
	 */
	private void calculateRepeatedIndexes(Stack loStack, Stack hiStack, final int size) {
		Arrays.fill(repeated, 0, size, false);

		while (!loStack.isEmpty()) {
			int lo = loStack.pop();
			int hi = hiStack.pop();
			for (int i = lo + 1; i <= hi; i++)
				repeated[i] = true;
		}
	}

	/**
	 * @param sort1Size
	 *            c sorted
	 * @param size
	 *            size of the data
	 * @return a partial suffix matrix with the c indexes
	 */
	private void calculatePartialSA(final int sort1Size, final int dataSize) {
		int j = 1;
		for (int i = 0; i < sort1Size; i++)
			partialSA[sort1[i] - 1] = j++;
	}

	/**
	 * specific comparison for the merge
	 * 
	 * @param valueC
	 *            value in array c
	 * @param valueB
	 *            value in array B
	 * @param data
	 *            the original data
	 * @param partialSA
	 *            partial suffix matrix with the c indexes
	 * @return valueB - valueC
	 */
	private int compare(int valueB, int valueC, final int[] data) {
		int ans;
		if (valueC % 3 == 1) {
			ans = data[valueB] - data[valueC];
			if (ans == 0)
				ans = partialSA[valueB] - partialSA[valueC];
		} else {
			ans = data[valueB] - data[valueC];
			if (ans == 0)
				ans = data[valueB + 1] - data[valueC + 1];
			if (ans == 0)
				ans = partialSA[valueB + 1] - partialSA[valueC + 1];
		}
		return ans;
	}

	/**
	 * 
	 * @param charSequence
	 *            the CharSequence with the letters
	 * @return an array where each letter in the sequence is assigned a value
	 *         according to its lexicographical order, starts from 1.
	 * 
	 *         the letter is indexed by the value ASCII - alphabet Start
	 */
	private void getMap(final CharSequence charSequence) {
		int[] map2 = new int[map.length];
		// mark the appearances whit 1
		for (int i = 0; i < charSequence.length(); i++)
			map2[charSequence.charAt(i)]++;

		StringBuilder str = new StringBuilder();
		firstRowsInMatrix = new TreeMap<>();
		lastRowsInMatrix = new TreeMap<>();
		firstRowsInMatrix.put(SPECIAL_CHARACTER, 0);
		lastRowsInMatrix.put(SPECIAL_CHARACTER, 0);

		byte cont = 1;
		int acum = 1;
		for (int i = 0; i < map.length; i++)
			if (map2[i] != 0) {
				char c = (char) i;
				firstRowsInMatrix.put(c, acum);
				acum += map2[i];
				lastRowsInMatrix.put(c, acum - 1);
				map[i] = cont++;
				str.append(c);
			}
		alphabet = str.toString();
	}

	/**
	 * 
	 * @param charSequence
	 *            the original char sequence
	 * @param map
	 *            the map with the value for each element
	 * @return use the map to transform the char sequence in a numbers array, the
	 *         last position represent the cut character $==0
	 */
	private static int[] transform(CharSequence charSequence) {
		int[] data = new int[charSequence.length() + 1];

		for (int i = 0; i < charSequence.length(); i++)
			data[i] = map[charSequence.charAt(i)];
		return data;
	}

	/**
	 * fill the array from startPoint with the values ​​whose module equals
	 * startValue that are less than maxValue
	 * 
	 * @param array
	 *            the array where the values ​​go
	 * @param startPoint
	 *            the point from where the array begins to fill
	 * @param startValue
	 *            the value to which the module has to be equal
	 * @param maxValue
	 *            the maximum possible value
	 * @return the maximum position that was filled +1
	 */
	private static int moduleThree(int[] array, final int startValue, final int startPoint, final int maxValue) {
		int C = startPoint;
		for (int i = startValue; i < maxValue; i += 3)
			array[C++] = i;
		return C;
	}

	/**
	 * this method uses the principle of the binary search
	 * 
	 * @param i
	 *            the number
	 * @return the position of the most significant bit
	 */
	private static int highestOneBitPos(final int i) {
		int lo = 0, hi = 31;
		int m = 0, ans;
		while (lo + 1 != hi) {
			m = (hi + lo) / 2;
			ans = i >>> m;
			if (ans == 0)
				hi = m;
			else
				lo = m;
		}
		return hi;
	}

	/**
	 * 
	 * @param i
	 *            the number to round
	 * @return the smallest multiple of nearest Cte_Radix Cte_Radix
	 */
	private static int round(final int i) {
		return (i / Cte_Radix) * Cte_Radix;
	}

	// -------------------------------------------------------------------------------
	// PRIVATE CLASS
	// -------------------------------------------------------------------------------
	/**
	 * This class represents an int stack
	 * 
	 * @author jc.bojaca
	 */
	private static class Stack {
		private Element apunt;

		public Stack() {
			apunt = null;
		}

		public boolean isEmpty() {
			return apunt == null;
		}

		public int pop() {
			int lo = apunt.val;
			apunt = apunt.sig;
			return lo;
		}

		public void add(int i) {
			Element temp = new Element(i);
			temp.sig = apunt;
			apunt = temp;
		}

		public static void exchange(Stack a, Stack b) {
			Element temp = a.apunt;
			a.apunt = b.apunt;
			b.apunt = temp;
		}

		private static class Element {
			Element sig;
			int val;

			public Element(int val) {
				this.val = val;
			}
		}
	}
}
