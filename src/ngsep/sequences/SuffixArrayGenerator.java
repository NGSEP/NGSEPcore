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

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;

/**
 * 
 * this is a test of generating the SA with the DC3 algorithm
 * 
 * 
 * @author Juan Camilo Bojaca
 *
 */
public class SuffixArrayGenerator {
	// -------------------------------------------------------------------------------
	// CONSTANT VALUES
	// -------------------------------------------------------------------------------
	/**
	 * models the radix factor (hex: 4 == log(16))
	 */
	private final static int Cte_Radix = 4;
	/**
	 * represents the maximum achievable value with CteRadix bites
	 */
	private final static int Cte_Radix_Bit = (int) Math.pow(2, Cte_Radix) - 1;
	/**
	 * the ASCII interval that it handles
	 */
	private final static int AlfabetEnd = 126;
	// -------------------------------------------------------------------------------
	// VARAIBLES
	// -------------------------------------------------------------------------------
	/**
	 * the suffix Array
	 */
	private int[] ans;

	// -------------------------------------------------------------------------------
	// CONSTRUCTOR && PUBLIC METHODS
	// -------------------------------------------------------------------------------
	/**
	 * 
	 * @param charSequence
	 *            the sequence to which the suffix array is calculated
	 */
	public SuffixArrayGenerator(CharSequence charSequence) {
		NumberFormat n = new DecimalFormat("0.000");
		long ini = System.currentTimeMillis();
		byte[] map = getMap(charSequence);
		int[] data = transform(charSequence, map);

		int[] sort1 = new int[(int) ((2 * data.length) / 3)];
		int[] sort2 = new int[((data.length - 1) / 3) + 1];
		int[] auxSort = new int[data.length];
		final int[] contSort = new int[Cte_Radix_Bit + 3];

		// map[map.length - 1] is the maximum possible value in data
		ans = getSuffix(data, round(highestOneBitPos(map[map.length - 1])), sort1, sort2, auxSort, contSort);
		System.out.println(n.format((System.currentTimeMillis() - ini) / (double) 1000));
	}

	/**
	 * 
	 * @return the suffix Array
	 */
	public int[] getSA() {
		int[] answer = new int[ans.length - 1];
		for (int i = 1; i < ans.length; i++) {
			answer[i - 1] = ans[i];
		}
		return answer;
	}

	/**
	 * 
	 * @param data
	 *            the chain to get the suffix array
	 * @param MaxBit
	 *            the maximum possible multiple Cte_Radix shift for at least one
	 *            element of data to be different from 0
	 * @return the suffix Array of data
	 */
	private static int[] getSuffix(final int[] data, final int MaxBit, final int[] sort1, final int[] sort2,
			final int[] auxSort, final int[] contSort) {
		// -----------------------------------------------------------------------------
		// STEP #1: create C (sort1)
		// -----------------------------------------------------------------------------
		// in sort[m] starts the module 2 values
		int m = moduleThree(sort1, 1, 0, data.length);
		int c = moduleThree(sort1, 2, m, data.length);
		// -----------------------------------------------------------------------------
		// STEP #2: radix sort C = b1b2
		// -----------------------------------------------------------------------------
		Stack loStack = new Stack(0), hiStack = new Stack(c - 1);

		int d = 0, bit = MaxBit;
		while (!loStack.isEmpty() && d < 3) {
			radixSort(sort1, loStack, hiStack, data, d, bit, auxSort, contSort);

			if (bit != 0)
				bit -= Cte_Radix;
			else {
				d++;
				bit = MaxBit;
			}
		}
		// -----------------------------------------------------------------------------
		// STEP #2.1: recursion
		// -----------------------------------------------------------------------------
		if (!loStack.isEmpty()) {
			boolean[] repeated = convertToBooleanArray(loStack, hiStack, c);
			int[] r = new int[c + 1];
			int maxValue = calculateR(r, m, sort1, c, repeated);

			int[] SApr = getSuffix(r, round(highestOneBitPos(maxValue)), sort1, sort2, auxSort, contSort);

			// using the Suffix array of r finishes order sort1
			for (int j = 1; j < SApr.length; j++) {
				int idx = SApr[j];
				sort1[j - 1] = (idx < m) ? idx * 3 + 1 : (idx - m) * 3 + 2;
			}
		}
		int[] partialSA = getPartialSA(sort1, c, data.length);
		int maxBitSA = round(highestOneBitPos(sort1.length));
		// -----------------------------------------------------------------------------
		// STEP #3: Sorting the non sample Suffices b0
		// -----------------------------------------------------------------------------
		int b = moduleThree(sort2, 0, 0, data.length);

		loStack = new Stack(0);
		hiStack = new Stack(b - 1);

		bit = MaxBit;
		while (loStack != null && bit >= 0) {
			radixSort(sort2, loStack, hiStack, data, 0, bit, auxSort, contSort);
			bit -= Cte_Radix;
		}
		bit = maxBitSA;
		while (loStack != null && bit >= 0) {
			radixSort(sort2, loStack, hiStack, partialSA, 1, bit, auxSort, contSort);
			bit -= Cte_Radix;
		}
		// -----------------------------------------------------------------------------
		// STEP #4: Merge
		// -----------------------------------------------------------------------------
		return merge(sort1, c, sort2, b, partialSA, data);
	}

	// -------------------------------------------------------------------------------
	// PRIVATE METHODS
	// -------------------------------------------------------------------------------

	/**
	 * 
	 * @param c
	 *            suffix array of c
	 * @param b
	 *            suffix array of b0
	 * @param partialSA
	 *            partial suffix matrix with the c indexes
	 * @param data
	 *            the original data
	 * @return the SA
	 */
	private static int[] merge(final int[] c, final int sizeC, final int[] b, final int sizeB, final int[] partialSA,
			final int[] data) {
		int[] ans = new int[data.length];
		int indexC = 0, indexB = 0, indexAns = 0;
		while (indexC != sizeC && indexB != sizeB) {
			if (compare(b[indexB], c[indexC], data, partialSA) < 0)
				ans[indexAns++] = b[indexB++];
			else
				ans[indexAns++] = c[indexC++];
		}
		while (indexC != sizeC)
			ans[indexAns++] = c[indexC++];
		while (indexB != sizeB)
			ans[indexAns++] = b[indexB++];
		return ans;
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
	private static int compare(int valueB, int valueC, final int[] data, final int[] partialSA) {
		int ans;
		if (valueC % 3 == 1) {
			ans = data[valueB] - data[valueC];
			if (ans == 0)
				ans = partialSA[valueB + 1] - partialSA[valueC + 1];
		} else {
			ans = data[valueB] - data[valueC];
			if (ans == 0)
				ans = data[valueB + 1] - data[valueC + 1];
			if (ans == 0)
				ans = partialSA[valueB + 2] - partialSA[valueC + 2];
		}
		return ans;
	}

	/**
	 * @param sort
	 *            c sorted
	 * @param size
	 *            size of the data
	 * @return a partial suffix matrix with the c indexes
	 */
	private static int[] getPartialSA(final int[] sort, int c, final int size) {
		int[] SA = new int[size];
		for (int j = 1, i = 0; i < c; i++, j++)
			SA[sort[i]] = j;
		return SA;

	}

	/**
	 * 
	 * @param r
	 *            for every element it contains the ranking obtained when ordering
	 *            the corresponding c(original values ​​of sort1) element
	 * @param m
	 *            in sort[m] starts the module 2 values
	 * @param sort
	 *            the sorted indexes
	 * @param repeated
	 *            array that represents the repeated elements
	 * @return the maximum value in r
	 */
	private static int calculateR(int[] r, final int m, final int[] sort, final int c, final boolean[] repeated) {
		int d = 0;
		for (int i = 0; i < c; ++i) {
			int idx = sort[i];
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
	private static boolean[] convertToBooleanArray(Stack loStack, Stack hiStack, final int size) {
		boolean[] repeated = new boolean[size];

		while (!loStack.isEmpty()) {
			int lo = loStack.pop();
			int hi = hiStack.pop();
			for (int i = lo + 1; i <= hi; i++)
				repeated[i] = true;
		}
		return repeated;
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
	private static byte[] getMap(final CharSequence charSequence) {
		byte[] map = new byte[AlfabetEnd];

		// mark the appearances whit 1
		for (int i = 0; i < charSequence.length(); i++) {
			char c = charSequence.charAt(i);
			if (map[c] == 0)
				map[c] = 1;
		}

		prefix(map);
		return map;
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
	private static int[] transform(CharSequence charSequence, final byte[] map) {
		int[] data = new int[charSequence.length() + 1];

		for (int i = 0; i < charSequence.length(); i++)
			data[i] = map[charSequence.charAt(i)];
		return data;
	}

	private static void radixSort(int[] array, Stack loStack, Stack hiStack, final int[] data, final int d,
			final int bit, final int[] auxSort, final int[] contSort) {
		// radixSKA(array, loStack, hiStack, data, d, bit);
		radixCount(array, loStack, hiStack, data, d, bit, auxSort, contSort);
	}

	/**
	 * https://www.youtube.com/watch?v=zqs87a_7zxw
	 * 
	 * @param array
	 *            the array of indexes to sort
	 * @param loStack,
	 *            hiStack contains the intervals to be sorted, ends with the
	 *            intervals where there are repeated elements
	 * @param data,
	 *            d, bit every element i of array is ordered according to the
	 *            function: (data[i + d] >>> bit) & Cte_Radix_Bit
	 * 
	 *            this extracts a hexadecimal number according to the bit shift of
	 *            the d-th data element after i.
	 */
	private static void radixSKA(int[] array, Stack loStack, Stack hiStack, final int[] data, final int d,
			final int bit) {
		final Stack loStackAux = new Stack(), hiStackAux = new Stack();
		int[][] c = new int[2][Cte_Radix_Bit + 1];

		while (!loStack.isEmpty()) {
			int lo = loStack.pop();
			int hi = hiStack.pop();

			Arrays.fill(c[0], 0);
			Arrays.fill(c[1], 0);

			count(c, lo, hi, array, data, d, bit);
			prefix(c);
			locateUsingRotations(c, lo, array, data, d, bit);
			extractIntervals(loStackAux, hiStackAux, c, lo, hi);
		}

		Stack.exchange(loStack, loStackAux);
		Stack.exchange(hiStack, hiStackAux);
	}

	/**
	 * using the range of the hexadecimal numbers, and the previous interval
	 * calculates the intervals where there are repeated elements in the array
	 * 
	 * @param loStack,
	 *            hiStack are used to store the new intervals
	 * @param c
	 *            the range of the hexadecimal numbers
	 * @param lo
	 *            the previous interval
	 */
	private static void extractIntervals(Stack loStack, Stack hiStack, int[][] c, final int lo, final int hi) {
		for (int act, prev = 0, i = 0; prev != hi - lo + 1; prev = act, i++) {
			act = c[1][i];
			int loAux = prev + lo;
			int hiAux = act - 1 + lo;
			if (loAux < hiAux) {
				loStack.add(loAux);
				hiStack.add(hiAux);
			}
		}
	}

	/**
	 * 
	 * sorts the array elements based on the function, using the intervals of each
	 * value.
	 * 
	 * @param c
	 *            the intervals of each value.
	 * @param lo
	 *            the initial position in the array
	 * @param array
	 *            the array
	 * @param data,
	 *            d, bit function: (data[i + d] >>> bit) & Cte_Radix_Bit
	 */
	private static void locateUsingRotations(int[][] c, final int lo, int[] array, final int[] data, final int d,
			final int bit) {
		boolean fin = false;
		while (!fin) {
			fin = true;
			for (int i = 0; i < c[0].length; i++) {
				fin &= c[0][i] == c[1][i];
				for (int j = c[0][i]; j < c[1][i]; j++)
					exchange(array, lo + j, lo + c[0][(data[array[j + lo] + d] >>> bit) & Cte_Radix_Bit]++);
			}
		}
	}

	/**
	 * 
	 * save in c[1] the number of hexadecimal elements as a result of the function
	 * on the elements of array from lo to hi
	 * 
	 * @param c
	 *            to save
	 * @param lo,
	 *            hi array the array and the interval
	 * @param data,
	 *            d, bit function: (data[i + d] >>> bit) & Cte_Radix_Bit
	 */
	private static void count(int[][] c, final int lo, final int hi, final int[] array, final int[] data, final int d,
			final int bit) {
		for (int i = lo; i <= hi; i++)
			c[1][(data[array[i] + d] >>> bit) & Cte_Radix_Bit]++;
	}

	/**
	 * 
	 * calculates the interval of each hexadecimal number based on the count
	 * 
	 * @param c
	 *            the count.
	 */
	private static void prefix(int[][] c) {
		for (int i = 1; i < Cte_Radix_Bit; i++)
			c[0][i + 1] = c[1][i] = c[1][i] + c[1][i - 1];
		c[0][1] = c[1][0];
		c[1][Cte_Radix_Bit] += c[1][Cte_Radix_Bit - 1];
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
	 * @param array
	 *            the array to calculate the prefix
	 */
	private static void prefix(byte[] array) {
		for (int i = 1; i < array.length; i++)
			array[i] += array[i - 1];
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

	/**
	 * swap the values of the indexes in array
	 * 
	 * @param array
	 *            the array
	 * @param i
	 *            first index
	 * @param j
	 *            second index
	 */
	private static void exchange(int[] array, final int i, final int j) {
		int t = array[i];
		array[i] = array[j];
		array[j] = t;
	}

	/**
	 * 
	 * @param array
	 *            the array of indexes to sort
	 * @param loStack,
	 *            hiStack contains the intervals to be sorted, ends with the
	 *            intervals where there are repeated elements
	 * @param data,
	 *            d, bit every element i of array is ordered according to the
	 *            function: (data[i + d] >>> bit) & Cte_Radix_Bit
	 * 
	 *            this extracts a hexadecimal number according to the bit shift of
	 *            the d-th data element after i.
	 */
	public static void radixCount(int[] array, Stack loStack, Stack hiStack, final int[] data, final int d,
			final int bit, final int[] auxSort, final int[] contSort) {
		final Stack loStackAux = new Stack(), hiStackAux = new Stack();

		while (!loStack.isEmpty()) {
			int lo = loStack.pop();
			int hi = hiStack.pop();

			Arrays.fill(contSort, 0);

			for (int i = lo; i <= hi; i++)
				contSort[((data[array[i] + d] >>> bit) & Cte_Radix_Bit) + 2]++;

			for (int i = 3; i < contSort.length; i++)
				contSort[i] += contSort[i - 1];

			for (int i = lo; i <= hi; i++)
				auxSort[contSort[((data[array[i] + d] >>> bit) & Cte_Radix_Bit) + 1]++] = array[i];

			for (int i = lo; i <= hi; i++)
				array[i] = auxSort[i - lo];

			for (int act, prev = contSort[0], i = 1; prev != hi - lo + 1; prev = act, i++) {
				act = contSort[i];
				int loAux = prev + lo;
				int hiAux = act - 1 + lo;
				if (loAux < hiAux) {
					loStackAux.add(loAux);
					hiStackAux.add(hiAux);
				}
			}
		}

		Stack.exchange(loStack, loStackAux);
		Stack.exchange(hiStack, hiStackAux);
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

		public Stack(int i) {
			apunt = null;
			Element temp = new Element(i);
			temp.sig = null;
			apunt = temp;
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
