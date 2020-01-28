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

import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * 
 * @author jc.bojaca
 *
 */
public class DC3SuffixArrayGenerator implements SuffixArrayGenerator {
	
	/** map of the letters */
	private final byte[] map = new byte[Byte.MAX_VALUE + 2];
	/** contains the number of occurrences of each hexadecimal digit */
	private final int[] contSort = new int[Byte.MAX_VALUE + 2];
	/** range to sort */
	private Deque<int[]> index;
	/** auxiliary range to sort */
	private Deque<int[]> indexAux;
	/** the suffix Array */
	private final int[] suffixArray;
	/** the ranks array and the auxiliary array for the dc3 */
	private int[] suffixArray2;
	/** the array to save the bytes values for data */
	private byte[] auxSort;
	/** the ranks array: Each position contains the position in the suffix array */
	private byte[] ranks;
	/** the letters in the sequence */
	private String alphabet;
	/** map of the letters */
	private Map<Character, Integer> countsA;

	/**
	 * Instantiates a new suffix array generator.
	 *
	 * @param charSequence the sequence to which the suffix array is calculated
	 */
	public DC3SuffixArrayGenerator(CharSequence sequence) {
		long ini = System.nanoTime();
		byte[] data = transform(sequence);
		System.err.println("Creating int arrays of length: "+data.length);
		index = new LinkedList<>();
		indexAux = new LinkedList<>();

		suffixArray = new int[data.length];
		suffixArray2 = new int[data.length];
		System.err.println("Created int arrays");
		auxSort = new byte[data.length];
		ranks = new byte[data.length * numerOfBytesOf((data.length << 1) / 3)];

		getSuffixArray(data, 1);
		//Dispose second array
		suffixArray2 = null;
		auxSort = null;
		ranks = null;

		System.err.println(
				"			Time to create the suffix array: " + (System.nanoTime() - ini) / ((double) 1000 * 1000 * 1000) + " s");
	}

	/**
	 * @return the suffix array. The first position is sequence.length
	 */
	public int [] getSuffixArray() {
		return suffixArray;
	}

	/**
	 * build the alphabet, the map, and a byte array from sequence
	 * 
	 * @param sequence the sequence
	 * @return an array with a number instead of each letter, the number keep the
	 *         lexicographical order
	 */
	private byte[] transform(final CharSequence sequence) {
		final int size = sequence.length();

		Arrays.fill(contSort, 0);
		StringBuilder alphabetSB = new StringBuilder();
		byte[] data = new byte[size + 1];

		// Calculates the counts per character and translates the characters to integers
		for (int i = 0; i < size; i++) {
			byte c = (byte) sequence.charAt(i);
			contSort[c]++;
			data[i] = c;
		}

		// Goes over the counts array to fill the alphabet and the map to transform the
		// raw data
		for (int i = 0; i < contSort.length; i++) {
			if (contSort[i] != 0) {
				alphabetSB.append((char) i);
				map[i] = (byte) alphabetSB.length();
			}
		}
		alphabet = alphabetSB.toString();
		for (int i = 0; i < size; i++)
			data[i] = map[data[i]];

		countsA = new HashMap<>(alphabet.length());
		for (int i = 0; i < alphabet.length(); i++) {
			char c = alphabet.charAt(i);
			countsA.put(c, contSort[c]);
		}
		return data;
	}

	/**
	 * 
	 * calculate the suffix array of data;
	 * 
	 * @param data  the chain to get the suffix array
	 * @param bytes the number of bytes for each value in data.
	 */
	private void getSuffixArray(byte[] data, int bytes) {
		final int m = arrayFill(0, 1 * bytes, data.length, 3 * bytes);
		final int sort1Size = arrayFill(m, 2 * bytes, data.length, 3 * bytes);
		index.push(new int[] { 0, sort1Size - 1 });
		for (int i = 0; i < 3 * bytes; i++)
			sort(data, i);
		div(0, sort1Size - 1, bytes);

		if (!index.isEmpty()) {
			final int a = sort1Size - markRepeated(sort1Size);
			final int numBytes = numerOfBytesOf(a);
			byte[] r = getR(numBytes, sort1Size, m);
			getSuffixArray(r, numBytes);
			returnOfRecurtion(sort1Size, m);
		}

		int byteRank = calculateRanks(sort1Size);

		final int sort2Size = arrayFill(sort1Size, 0, data.length, 3 * bytes);
		index.push(new int[] { sort1Size, sort2Size - 1 });
		for (int i = 0; i < bytes; i++)
			sort(data, i);
		divMul(sort1Size, sort2Size - 1, bytes, byteRank);
		for (int i = 0; i < byteRank; i++)
			sort(ranks, i);
		div(sort1Size, sort2Size - 1, byteRank);

		merge(data, bytes, byteRank, sort1Size, sort2Size);
	}

	/**
	 * divide and multiply the range in suffixArray2.
	 * 
	 * @param lo     start point
	 * @param hi     final point
	 * @param value  divide
	 * @param value2 multiply
	 */
	private void divMul(int lo, int hi, int value, int value2) {
		if (value == value2)
			return;
		for (int i = lo; i <= hi; i++) {
			suffixArray2[i] /= value;
			suffixArray2[i] *= value2;
		}
	}

	/**
	 * divide the range in suffixArray2.
	 * 
	 * @param lo    start point
	 * @param hi    final point
	 * @param value the factor
	 */
	private void div(int lo, int hi, int value) {
		if (value == 1)
			return;
		for (int i = lo; i <= hi; i++)
			suffixArray2[i] /= value;
	}

	/**
	 * transform the suffixArray of r in non sample suffixes
	 * 
	 * @param size the size of the suffix array
	 * @param m    the start point of module 2 values.
	 */
	private void returnOfRecurtion(int size, int m) {
		int idx = 0;
		for (int j = 1; j <= size; j++) {
			idx = suffixArray[j];
			if (idx < m)
				suffixArray2[j - 1] = idx * 3 + 1;
			else
				suffixArray2[j - 1] = (idx - m) * 3 + 2;
		}
	}

	/**
	 * merge the sample and non sample suffixes
	 * 
	 * @param data      the data to make the suffix array
	 * @param byteData  the number of bytes for each element in data
	 * @param byteRank  the number of bytes for each element in ranks
	 * @param sort1Size the size of non sample suffixes
	 * @param sort2Size the size of non suffixes
	 */
	private void merge(final byte[] data, int byteData, int byteRank, int sort1Size, int sort2Size) {
		int index2 = sort1Size;
		int index1 = 0;
		int indexAns = 0;
		int valueB = 0;
		int ans = 0;
		int valueC = suffixArray2[0];
		int valueCMod = valueC % 3;

		General: while (index2 != sort2Size) {
			valueB = suffixArray2[index2];

			while (true) {
				ans = 0;
				for (int x = 0, i = valueB * byteData, j = valueC * byteData; ans == 0 && x < byteData; x++, j++, i++)
					ans = data[i] - data[j];

				if (ans == 0) {
					if (valueCMod == 1) {
						for (int x = 0, i = valueB * byteRank, j = valueC * byteRank; ans == 0
								&& x < byteRank; x++, j++, i++)
							ans = ranks[i] - ranks[j];
					} else {
						for (int x = 0, i = (valueB + 1) * byteData, j = (valueC + 1) * byteData; ans == 0
								&& x < byteData; x++, j++, i++)
							ans = data[i] - data[j];
						for (int x = 0, i = (valueB + 1) * byteRank, j = (valueC + 1) * byteRank; ans == 0
								&& x < byteRank; x++, j++, i++)
							ans = ranks[i] - ranks[j];
					}
				}

				if (ans > 0) {
					suffixArray[indexAns++] = valueC;
					index1++;
					if (index1 == sort1Size)
						break General;
					valueC = suffixArray2[index1];
					valueCMod = valueC % 3;
				} else {
					break;
				}
			}

			suffixArray[indexAns++] = valueB;
			index2++;
		}

		for (; index1 != sort1Size; indexAns++, index1++)
			suffixArray[indexAns] = suffixArray2[index1];
		for (; index2 != sort2Size; indexAns++, index2++)
			suffixArray[indexAns] = suffixArray2[index2];
	}

	/**
	 * the number of bytes necessary to represent the value.
	 * 
	 * @param value the value
	 * @return the number of bytes.
	 */
	private static int numerOfBytesOf(int value) {
		int v = 1;
		int corr = 7;
		while (value >> corr != 0) {
			v++;
			corr += 7;
		}
		return v;
	}

	/**
	 * transform the order in the non sample suffixes in a ranks
	 * 
	 * @param size size of non sample suffixes
	 * @return the the number of bytes necessary to represent an element.
	 */
	private int calculateRanks(int size) {
		int bytes = numerOfBytesOf(size);
		byte[] cont = new byte[bytes];

		int j;
		for (int i = 0; i < size; i++) {
			pluss(cont, bytes - 1);
			j = (suffixArray2[i] - 1) * bytes;

			for (int x = 0; x < bytes; x++)
				ranks[j++] = cont[x];
		}
		return bytes;
	}

	/**
	 * 
	 * create the recursion array
	 * 
	 * @param m        in suffixArray2[m] starts the module 2 values
	 * @param numBytes the number of bytes for each element.
	 * @return recursion array
	 */
	private byte[] getR(int numBytes, int size, int m) {
		byte[] r = new byte[(size + 1) * numBytes];
		byte[] cont = new byte[numBytes];
		int j;
		int idx = 0;
		for (int i = 0; i < size; ++i) {
			idx = suffixArray2[i];
			if (auxSort[i] == 0)
				pluss(cont, numBytes - 1);

			j = idx / 3;
			if (idx % 3 == 2)
				j += m;

			j *= numBytes;
			for (int k = 0; k < numBytes; k++) {
				r[j++] = cont[k];
			}
		}
		return r;
	}

	/**
	 * change the format of repeated intervals to a byte array
	 * 
	 * @param size the size of non sample suffixes
	 * @return the number of repeated elements
	 */
	private int markRepeated(int size) {
		Arrays.fill(auxSort, 0, size, (byte) 0);
		int total = 0;

		int lo;
		int hi;
		int[] range;
		while (!index.isEmpty()) {
			range = index.pop();
			lo = range[0];
			hi = range[1];
			total += hi - lo;
			for (int i = lo + 1; i <= hi; i++)
				auxSort[i] = 1;
		}
		return total;
	}

	/**
	 * sort the suffixes of suffixArra2 in the range in index based in the byte of
	 * data + mov.
	 * 
	 * @param data the data indexed in suffixArra2
	 * @param mov  the plus (to sort more bytes than one)
	 */
	private void sort(byte[] data, int mov) {
		int prev;
		int cont;
		int act;
		byte value;
		int lo;
		int hi;
		int[] range;
		while (!index.isEmpty()) {
			range = index.pop();
			lo = range[0];
			hi = range[1];

			Arrays.fill(contSort, 0);
			for (int i = lo; i <= hi; i++) {
				value = data[suffixArray2[i] + mov];
				auxSort[i] = value;
				contSort[value]++;
			}

			cont = contSort[0] += lo - 1;
			prev = 0;
			for (int i = 1; prev != hi; i++) {
				prev = cont;
				contSort[i] = cont += contSort[i];
			}

			for (int i = lo; i <= hi; i++)
				suffixArray[contSort[auxSort[i]]--] = suffixArray2[i];

			System.arraycopy(suffixArray, lo, suffixArray2, lo, (hi - lo + 1));

			prev = contSort[0];
			for (int i = 1; prev != hi; prev = act, i++) {
				act = contSort[i];
				if (prev + 1 < act)
					indexAux.push(new int[] { prev + 1, act });
			}
		}
		Deque<int[]> temp = index;
		index = indexAux;
		indexAux = temp;
	}

	/**
	 * fill the array from startPoint with the values less than maxValue
	 * 
	 * @param pos      the point from where the array begins to fill
	 * @param lo       the value to which the module has to be equal
	 * @param maxValue the max value
	 * @param inc      the increment.
	 * @return
	 */
	private int arrayFill(final int pos, final int lo, final int maxValue, int inc) {
		int i = pos;
		for (int value = lo; value < maxValue; value += inc)
			suffixArray2[i++] = value;
		return i;
	}

	/**
	 * default sum of array like a number.
	 * 
	 * @param array the array to sum
	 * @param i     the index to sum
	 */
	private static final void pluss(byte[] array, int i) {
		array[i]++;
		if (array[i] < 0 && i != 0) {
			array[i] = 0;
			pluss(array, i - 1);
		}
	}
}
