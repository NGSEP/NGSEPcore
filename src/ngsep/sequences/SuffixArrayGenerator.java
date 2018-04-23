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
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * @author jc.bojaca
 * @version 1
 */
interface Constants {
	/**
	 * character to BWT
	 */
	static final char SPECIAL_CHARACTER = '$';
	/**
	 * the ASCII interval that it handles
	 */
	static final int ALFABETEND = 127;
	/**
	 * models the radix factor (bits used to sort)
	 */
	final static int Cte_Radix = 8;
	/**
	 * represents the maximum achievable value with CteRadix bites
	 */
	final static int Cte_Radix_Bit = (int) Math.pow(2, Cte_Radix) - 1;
	/**
	 * DC#3
	 */
	static final int CONS = 3;
	/**
	 * number format
	 */
	static final DecimalFormat format = new DecimalFormat("0.000000000");

	/**
	 * this method uses the principle of the binary search
	 * 
	 * @param i
	 *            the number
	 * @return the position of the most significant bit
	 */
	default int highestOneBitPos(final int i) {
		int lo = 0;
		int hi = (Integer.BYTES << CONS) - 1;
		int med = 0;
		while (lo + 1 != hi) {
			med = (hi + lo) >> 1;
			switch (i >>> med) {
			case 0:
				hi = med;
				break;
			default:
				lo = med;
				break;
			}
		}
		return hi;
	}

	/**
	 * 
	 * @param i
	 *            the number to round
	 * @return the smallest multiple of nearest Cte_Radix Cte_Radix
	 */
	default int round(final int i) {
		return (i / Cte_Radix) * Cte_Radix;
	}
}

/**
 * @author jc.bojaca
 * @version
 */
public class SuffixArrayGenerator implements Constants {
	/**
	 * the char sequence
	 */
	final private CharSequence sequence;
	/**
	 * the suffix Array
	 */
	final private int[] suffixArray, ranksSA;
	/**
	 * the translator of the sequence
	 */
	final private Translator translator;

	/**
	 * 
	 * @param charSequence
	 *            the sequence to which the suffix array is calculated
	 */
	public SuffixArrayGenerator(CharSequence charSequence) {
		sequence = charSequence;
		translator = new Translator(charSequence);
		int[] data = translator.getData();
		DC3 dc3 = new DC3(data);
		suffixArray = dc3.getSuffixArray();
		ranksSA = dc3.getRanksSA();
	}

	public int[] getSA() {
		final int[] answer = new int[suffixArray.length - 1];
		System.arraycopy(suffixArray, 1, answer, 0, answer.length);
		return answer;
	}

	public char[] getBWT() {
		final char[] bwt = new char[suffixArray.length];
		bwt[ranksSA[0]] = SPECIAL_CHARACTER;
		for (int i = 1; i < suffixArray.length; i++)
			bwt[ranksSA[i]] = sequence.charAt(i - 1);
		return bwt;
	}

	public Map<Integer, Integer> getPartialSuffixArray(final int suffixFraction) {
		final Map<Integer, Integer> partialSuffixArray = new HashMap<Integer, Integer>();
		for (int i = 0; i < ranksSA.length; i += suffixFraction)
			partialSuffixArray.put(ranksSA[i], i);
		return partialSuffixArray;
	}

	public int[][] getTallyIndexes(int tallyDistance, char[] bwt) {
		final int alfabetLength = translator.getAlphabet().length();
		final int tallyRows = (bwt.length + (tallyDistance - 1)) / tallyDistance;

		final byte[] map = translator.getMap();
		final int[] arr = new int[alfabetLength + 1];
		final int[][] tallyIndexes = new int[tallyRows][alfabetLength];

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
		return translator.getAlphabet().toString();
	}

	public Map<Character, Integer> getFirstRowsInMatrix() {
		return translator.getFirstRowsInMatrix();
	}

	public Map<Character, Integer> getLastRowsInMatrix() {
		return translator.getLastRowsInMatrix();
	}
}

/**
 * 
 * this class provides support for the ranks chain
 * 
 * @author jc.bojaca
 * @version 1
 */
class DC3 implements Constants {
	/**
	 * contains the number of occurrences of each hexadecimal digit
	 */
	private final static int[] contSort = new int[Cte_Radix_Bit + 2];
	/**
	 * the suffix Array
	 */
	final private int[] suffixArray, ranksSA;

	/**
	 * is used to sort
	 */
	private final boolean[] repeated;
	private final byte[] auxSort2;
	private final int[] auxSort, sort1, sort2;
	private final Stack loStack = new Stack(), hiStack = new Stack(), loStackAux = new Stack(),
			hiStackAux = new Stack();

	DC3(int[] data) {
		final int size = data.length;
		final int size2_of_3 = (2 * size) / 3;
		final int size1_of_3 = ((size - 1) / 3) + 1;

		suffixArray = new int[size];
		ranksSA = new int[size];

		sort2 = new int[size1_of_3];

		sort1 = new int[size2_of_3];
		auxSort = new int[size2_of_3];
		auxSort2 = new byte[size2_of_3];
		repeated = new boolean[size2_of_3];

		getSuffix(data, 0);

		for (int i = 0; i < suffixArray.length; i++)
			ranksSA[suffixArray[i]] = i;
	}

	/**
	 * @return the suffixArray
	 */
	int[] getSuffixArray() {
		return suffixArray;
	}

	/**
	 * @return the suffixArray
	 */
	int[] getRanksSA() {
		return ranksSA;
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
	private void getSuffix(final int[] data, final int MaxBit) {
		// -----------------------------------------------------------------------------
		// STEP #1: create C (sort1)
		// -----------------------------------------------------------------------------
		final int m = moduleThree(sort1, 0, 1, data.length);
		final int sort1Size = moduleThree(sort1, m, 2, data.length);
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
			final int[] r = new int[sort1Size + 1];
			final int maxValueR = calculateR(r, m, sort1Size);

			getSuffix(r, round(highestOneBitPos(maxValueR)));

			int idx = 0;
			// using the Suffix array of r, finishes order sort1
			for (int j = 1; j <= sort1Size; j++) {
				idx = suffixArray[j];
				if (idx < m) {
					sort1[j - 1] = idx * 3 + 1;
				} else {
					sort1[j - 1] = (idx - m) * 3 + 2;
				}
			}
		}
		calculatePartialSA(sort1Size, data.length);
		// -----------------------------------------------------------------------------
		// STEP #3: Sorting the non sample Suffices b0
		// -----------------------------------------------------------------------------
		final int sort2Size = moduleThree(sort2, 0, 0, data.length);
		changeSortInterval(0, sort2Size - 1);

		for (int bit = MaxBit; !loStack.isEmpty() && bit >= 0; bit -= Cte_Radix)
			sort(sort2, data, 0, bit);
		for (int bit = round(highestOneBitPos(sort1Size)); !loStack.isEmpty() && bit >= 0; bit -= Cte_Radix)
			sort(sort2, ranksSA, 0, bit);
		// -----------------------------------------------------------------------------
		// STEP #4: Merge
		// -----------------------------------------------------------------------------
		long ini = System.nanoTime();
		merge(sort1Size, sort2Size, data);
		System.out.println(format.format((System.nanoTime() - ini) / ((double) 1000 * 1000 * 1000)));

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
		int ind = 0;
		int hi = 0;
		int lo = 0;
		while (!loStack.isEmpty()) {
			lo = loStack.pop();
			hi = hiStack.pop();

			Arrays.fill(contSort, 0);

			for (int i = lo; i <= hi; i++) {
				ind = data[array[i] + d] >>> bit;
				auxSort2[i] = (byte) ind;
				contSort[Cte_Radix_Bit & ind]++;
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
	 * @param data
	 *            the original data
	 * @return the SA
	 */
	private int[] merge(final int sort1Size, final int sort2Size, final int[] data) {
		int index1 = 0, index2 = 0, indexAns = 0;
		int valueB = 0;
		int ans = 0;
		int valueC = sort1[0];
		int valueCMod = valueC % CONS;

		General: while (index2 != sort2Size) {
			valueB = sort2[index2];
			while (true) {
				switch (valueCMod) {
				case 1:
					ans = data[valueB] - data[valueC];
					if (0 == ans)
						ans = ranksSA[valueB] - ranksSA[valueC];
					break;

				default:
					ans = data[valueB] - data[valueC];
					if (0 == ans)
						ans = data[valueB + 1] - data[valueC + 1];
					if (0 == ans)
						ans = ranksSA[valueB + 1] - ranksSA[valueC + 1];
					break;
				}
				if (ans > 0) {
					suffixArray[indexAns++] = valueC;
					index1++;
					if (index1 == sort1Size)
						break General;
					valueC = sort1[index1];
					valueCMod = valueC % CONS;
				} else
					break;
			}

			suffixArray[indexAns++] = valueB;
			index2++;
		}

		for (; index1 != sort1Size; indexAns++, index1++)
			suffixArray[indexAns] = sort1[index1];
		for (; index2 != sort2Size; indexAns++, index2++)
			suffixArray[indexAns] = sort2[index2];

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
	 * @return the maximum value in r
	 */
	private int calculateR(int[] r, final int m, final int sort1Size) {
		int d = 0;
		int idx = 0;
		for (int i = 0; i < sort1Size; ++i) {
			idx = sort1[i];
			if (!repeated[i]) {
				++d;
			}
			switch (idx % CONS) {
			case 1:
				r[idx / CONS] = d;
				break;

			default:
				r[idx / CONS + m] = d;
				break;
			}
		}
		return d;
	}

	/**
	 * 
	 * change the format of repeated intervals to a Boolean array
	 * 
	 * @param loStack
	 *            repeated intervals
	 * @param hiStack
	 *            repeated intervals
	 * @param size
	 *            of the array
	 */
	private void calculateRepeatedIndexes(Stack loStack, Stack hiStack, final int size) {
		int lo = 0;
		int hi = 0;

		Arrays.fill(repeated, 0, size, false);

		while (!loStack.isEmpty()) {
			lo = loStack.pop();
			hi = hiStack.pop();
			for (int i = lo + 1; i <= hi; i++) {
				repeated[i] = true;
			}
		}
	}

	/**
	 * @param sort1Size
	 *            c sorted
	 * @param dataSize
	 *            size of the data
	 */
	private void calculatePartialSA(final int sort1Size, final int dataSize) {
		int j = 1;
		for (int i = 0; i < sort1Size; i++) {
			ranksSA[sort1[i] - 1] = j;
			j++;
		}
	}

	/**
	 * fill the array from startPoint with the values ​​whose module equals
	 * startValue that are less than maxValue
	 * 
	 * @param array
	 *            the array where the values ​​go
	 * @param lo
	 *            the point from where the array begins to fill
	 * @param factor
	 *            the value to which the module has to be equal
	 * @param maxValue
	 *            the maximum possible value
	 * @return the maximum position that was filled +1
	 */
	private static int moduleThree(int[] array, final int lo, final int factor, final int maxValue) {
		int index = lo;
		for (int value = factor; value < maxValue; value += CONS)
			array[index++] = value;
		return index;
	}
}

/**
 * 
 * this class provides support for the ranks chain
 * 
 * @author jc.bojaca
 * @version 1
 */
class Translator implements Constants {
	/**
	 * the array with the transformation
	 */
	private final int[] data;
	/**
	 * the letters in the sequence
	 */
	private final StringBuilder alphabet;
	/**
	 * map of the letters
	 */
	private final byte[] map = new byte[ALFABETEND];
	/**
	 * number of each letter
	 */
	private final static int[] counts = new int[ALFABETEND];
	/**
	 * the position of the final index in the count
	 */
	private final Map<Character, Integer> lastRowsInMatrix = new TreeMap<>();
	/**
	 * the position of the initial index in the count
	 */
	private final Map<Character, Integer> firstRowsInMatrix = new TreeMap<>();

	/**
	 * constructor
	 * 
	 * @param charSequence
	 *            the original sequence
	 */
	Translator(CharSequence charSequence) {
		final int size = charSequence.length();
		int c = 0;

		Arrays.fill(counts, 0);

		firstRowsInMatrix.put(SPECIAL_CHARACTER, 0);
		lastRowsInMatrix.put(SPECIAL_CHARACTER, 0);
		alphabet = new StringBuilder(ALFABETEND);
		data = new int[size + 1];

		for (int i = 0; i < size; i++) {
			c = (int) charSequence.charAt(i);
			counts[c]++;
			data[i] = c;
		}

		for (int acum = 1, cont = 1, i = 0; i < ALFABETEND; i++)
			if (0 != counts[i]) {
				firstRowsInMatrix.put((char) i, acum);
				lastRowsInMatrix.put((char) i, (acum += counts[i]) - 1);
				map[i] = (byte) cont++;
				alphabet.append((char) i);
			}

		for (int i = 0; i < size; i++)
			data[i] = map[data[i]];
	}

	/**
	 * @return the data
	 */
	int[] getData() {
		return data;
	}

	/**
	 * @return the data
	 */
	byte[] getMap() {
		return map;
	}

	/**
	 * @return the alphabet
	 */
	StringBuilder getAlphabet() {
		return alphabet;
	}

	/**
	 * @return the firstRowsInMatrix
	 */
	Map<Character, Integer> getFirstRowsInMatrix() {
		return firstRowsInMatrix;
	}

	/**
	 * @return the lastRowsInMatrix
	 */
	Map<Character, Integer> getLastRowsInMatrix() {
		return lastRowsInMatrix;
	}

	@Override
	public String toString() {
		final StringBuilder str = new StringBuilder();
		str.append(alphabet.toString());
		str.append(System.getProperty("line.separator"));
		str.append(firstRowsInMatrix.toString());
		str.append(System.getProperty("line.separator"));
		str.append(lastRowsInMatrix.toString());
		return str.toString();
	}
}

/**
 * This class represents an int stack
 * 
 * @author jc.bojaca
 * @version 1
 */
class Stack {
	/**
	 * the firts element of the stack
	 */
	private Element apunt;

	/**
	 * the constructor
	 */
	Stack() {
		apunt = null;
	}

	/**
	 * 
	 * @return if the stack not have values.
	 */
	boolean isEmpty() {
		return null == apunt;
	}

	/**
	 * get the first element into stack
	 * 
	 * @return the first element
	 */
	int pop() {
		final int value = apunt.val;
		apunt = apunt.sig;
		return value;
	}

	/**
	 * add one element into stack
	 * 
	 * @param i
	 *            the element
	 */
	void add(int i) {
		final Element toAdd = new Element(i);
		toAdd.sig = apunt;
		apunt = toAdd;
	}

	/**
	 * to exchange to stacks
	 * 
	 * @param a
	 *            stack 1
	 * @param b
	 *            stack 2
	 */
	static void exchange(Stack a, Stack b) {
		final Element temp = a.apunt;
		a.apunt = b.apunt;
		b.apunt = temp;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		final StringBuilder str = new StringBuilder("[");
		Element element = null;
		for (element = apunt; null != element.sig; element = element.sig) {
			str.append(element.toString());
			str.append(',');
		}
		str.append(element.toString());
		str.append(']');
		return str.toString();
	}

	/**
	 * this class represents one element of stack
	 * 
	 * @author jc.bojaca
	 *
	 */
	private static class Element {
		/**
		 * the next element
		 */
		private Element sig;
		/**
		 * the value
		 */
		private final int val;

		/**
		 * constructor
		 */
		private Element() {
			val = 0;
		}

		/**
		 * constructor
		 * 
		 * @param val
		 *            the value
		 */
		private Element(int val) {
			this.val = val;
		}

		/**
		 * to string
		 * 
		 * @return the string of the translator
		 */
		@Override
		public String toString() {
			return String.valueOf(val);
		}
	}
}