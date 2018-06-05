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
import java.util.Map;
import java.util.TreeMap;
import java.util.function.IntUnaryOperator;

/**
 * The Class SuffixArrayGenerator.
 *
 * @author jc.bojaca
 * @version 1
 */
public class SuffixArrayGenerator {
    /** the number of bytes to sort each character */
    public static final int BYTESTOSORT = 16;

    /** Character to BWT */
    public static final char SPECIAL_CHARACTER = '$';

    /** The ASCII interval that it handles */
    public static final byte LAST_ASCII_CHARACTER = 127;

    /** Models the radix factor (bits used to sort) */
    private static final int CTE_RADIX = 8;

    /** Represents the maximum achievable value with CteRadix bites */
    private static final int CTE_RADIX_BIT = (int) Math.pow(2, CTE_RADIX) - 1;

    /** Difference cover CONS modulo */
    private static final int CONS = 3;

    /** Number format to print */
    private static final DecimalFormat FORMAT = new DecimalFormat("0.000000000");

    /** The position of the most significant bit round to Cte_Radix Cte_Radix */
    public static final IntUnaryOperator HighestMostSignificatBitRound = (final int i) -> {
	int hi = (Integer.BYTES << CONS) - 1;
	int med = 0;
	int lo = 0;
	while (lo + 1 != hi) {
	    med = (hi + lo) >>> 1;
	    if (0 == i >>> med)
		hi = med;
	    else
		lo = med;
	}
	return (hi / CTE_RADIX) * CTE_RADIX;
    };

    /** the letters in the sequence */
    private String alphabet;

    /** map of the letters */
    private final byte[] map = new byte[SuffixArrayGenerator.LAST_ASCII_CHARACTER];

    /** number of each letter */
    private final int[] counts = new int[SuffixArrayGenerator.LAST_ASCII_CHARACTER];

    /** the suffix Array */
    private final int[] suffixArray;
    /** the ranks array: Each position contains the position in the suffix array */
    private final int[] ranksSA;

    /** contains the number of occurrences of each hexadecimal digit */
    private final int[] contSort = new int[CTE_RADIX_BIT + 2];
    /** the elements repeated in the pseudo ordering */
    private final boolean[] repeated;
    /** the array to save the function values for data */
    private final byte[] auxSort2;
    /** the array to save and copy the sort */
    private final int[] auxSort;
    /** the sample Suffices */
    private final int[] sort1;
    /** the non sample Suffices */
    private final int[] sort2;
    /** start points to sort */
    private final Stack loStack = new Stack();
    /** final points to sort */
    private final Stack hiStack = new Stack();
    /** start auxiliary points to sort */
    private final Stack loStackAux = new Stack();
    /** final auxiliary points to sort */
    private final Stack hiStackAux = new Stack();

    /**
     * Instantiates a new suffix array generator.
     *
     * @param charSequence
     *            the sequence to which the suffix array is calculated
     */
    public SuffixArrayGenerator(CharSequence sequence) {
	long ini = System.nanoTime();
	int[] data = preprocess(sequence);
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

	// TODO: check if it has to be recalculated
	for (int i = 0; i < suffixArray.length; i++)
	    ranksSA[suffixArray[i]] = i;

	System.out.println("Time to create a suffix array: "
		+ FORMAT.format((System.nanoTime() - ini) / ((double) 1000 * 1000 * 1000)) + " s");
    }

    private int[] preprocess(CharSequence sequence) {
	final int size = sequence.length();

	Arrays.fill(counts, 0);
	StringBuilder alphabetSB = new StringBuilder();
	int[] data = new int[size + 1];

	// Calculates the counts per character and translates the characters to integers
	for (int i = 0; i < size; i++) {
	    int c = sequence.charAt(i);
	    counts[c]++;
	    data[i] = c;
	}

	// Goes over the counts array to fill the alphabet and the map to transform the
	// raw data
	for (int i = 0; i < counts.length; i++) {
	    if (counts[i] != 0) {
		alphabetSB.append((char) i);
		map[i] = (byte) alphabetSB.length();
	    }
	}
	alphabet = alphabetSB.toString();
	for (int i = 0; i < size; i++)
	    data[i] = map[data[i]];
	return data;
    }

    /**
     * @return the suffix array. The first position is sequence.length
     */
    public int[] getSuffixArray() {
	return suffixArray;
    }

    /**
     * @return the reverse of the suffix array
     */
    public int[] getReverseSuffixArray() {
	return ranksSA;
    }

    /**
     * @return the alphabet
     */
    public String getAlphabet() {
	return alphabet;
    }

    /**
     * @return The number of times that each character in the alphabet appears
     */
    public Map<Character, Integer> getCharacterCounts() {
	Map<Character, Integer> countsA = new TreeMap<>();
	for (int i = 0; i < alphabet.length(); i++) {
	    char c = alphabet.charAt(i);
	    countsA.put(c, counts[c]);
	}
	return countsA;
    }

    /**
     * calculate the suffix array of data;
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
	    for (int bit = MaxBit; !loStack.isEmpty() && bit >= 0; bit -= CTE_RADIX)
		sort(sort1, data, d, bit);
	// -----------------------------------------------------------------------------
	// STEP #2.1: recursion
	// -----------------------------------------------------------------------------
	if (!loStack.isEmpty()) {
	    calculateRepeatedIndexes(loStack, hiStack, sort1Size);
	    final int[] r = new int[sort1Size + 1];
	    final int maxValueR = calculateR(r, m, sort1Size);

	    getSuffix(r, HighestMostSignificatBitRound.applyAsInt(maxValueR));

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
	generetareTheRanksArray(sort1Size);
	// -----------------------------------------------------------------------------
	// STEP #3: Sorting the non sample Suffices b0
	// -----------------------------------------------------------------------------
	final int sort2Size = moduleThree(sort2, 0, 0, data.length);
	changeSortInterval(0, sort2Size - 1);

	for (int bit = MaxBit; !loStack.isEmpty() && bit >= 0; bit -= CTE_RADIX)
	    sort(sort2, data, 0, bit);
	for (int bit = HighestMostSignificatBitRound.applyAsInt(sort1Size); !loStack.isEmpty()
		&& bit >= 0; bit -= CTE_RADIX)
	    sort(sort2, ranksSA, 0, bit);
	// -----------------------------------------------------------------------------
	// STEP #4: Merge
	// -----------------------------------------------------------------------------
	merge(sort1Size, sort2Size, data);
    }

    /**
     * define the interval to sort
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
     * sort values in array based in the specific function: (data[i + d] >>> bit) &
     * Cte_Radix_Bit
     * 
     * @param array
     *            the array of indexes to sort
     * @param loStack
     *            , hiStack contains the intervals to be sorted, ends with the
     *            intervals where there are repeated elements
     * @param data
     *            the original array
     * @param d
     *            the byte reference
     * @param bit
     *            the bite reference
     */
    private void sort(int[] array, final int[] data, final int d, final int bit) {
	radixCount(array, data, d, bit);
    }

    /**
     * sort values in array based in the specific function: (data[i + d] >>> bit) &
     * Cte_Radix_Bit
     * 
     * @param array
     *            the array of indexes to sort
     * @param loStack
     *            , hiStack contains the intervals to be sorted, ends with the
     *            intervals where there are repeated elements
     * @param data
     *            the original array
     * @param d
     *            the byte reference
     * @param bit
     *            the bite reference
     */
    private void radixCount(int[] array, final int[] data, final int d, final int bit) {
	int ind = 0;
	int hi = 0;
	int lo = 0;
	int prev;
	int cont;
	int act;
	while (!loStack.isEmpty()) {
	    lo = loStack.pop();
	    hi = hiStack.pop();

	    Arrays.fill(contSort, 0);

	    for (int i = lo; i <= hi; i++) {
		ind = data[array[i] + d] >>> bit;
		auxSort2[i] = (byte) ind;
		contSort[CTE_RADIX_BIT & ind]++;
	    }

	    contSort[0] += lo - 1;
	    prev = 0;
	    cont = contSort[0];
	    for (int i = 1; prev != hi; i++) {
		prev = cont;
		contSort[i] = cont += contSort[i];
	    }

	    for (int i = lo; i <= hi; i++) {
		auxSort[contSort[Byte.toUnsignedInt(auxSort2[i])]--] = array[i];
	    }

	    System.arraycopy(auxSort, lo, array, lo, (hi - lo + 1));

	    prev = contSort[0];
	    for (int i = 1; prev != hi; prev = act, i++) {
		act = contSort[i];
		if (prev + 1 < act) {
		    loStackAux.add(prev + 1);
		    hiStackAux.add(act);
		}
	    }
	}
	loStack.exchange(loStackAux);
	hiStack.exchange(hiStackAux);
    }

    /**
     * merge sort1 and sort2
     * 
     * @param sort1Size
     *            suffix array of c
     * @param sort2Size
     *            suffix array of b0
     * @param data
     *            the original data (contains the values referenced in sort1 and
     *            sort2)
     * @return the SA
     */
    private int[] merge(final int sort1Size, final int sort2Size, final int[] data) {
	int index1 = 0;
	int index2 = 0;
	int indexAns = 0;
	int valueB = 0;
	int ans = 0;
	int valueC = sort1[0];
	int valueCMod = valueC % CONS;

	General: while (index2 != sort2Size) {
	    valueB = sort2[index2];
	    while (true) {
		if (valueCMod == 1) {
		    ans = data[valueB] - data[valueC];
		    if (0 == ans)
			ans = ranksSA[valueB] - ranksSA[valueC];
		} else {
		    ans = data[valueB] - data[valueC];
		    if (0 == ans)
			ans = data[valueB + 1] - data[valueC + 1];
		    if (0 == ans)
			ans = ranksSA[valueB + 1] - ranksSA[valueC + 1];
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
     * create the recursion array
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
	    if (idx % CONS == 1) {
		r[idx / CONS] = d;
	    } else {
		r[idx / CONS + m] = d;
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
     * transform the order in the sort1 in a ranks in RA
     * 
     * @param sort1Size
     *            c sorted
     * @param dataSize
     *            size of the data
     */
    private void generetareTheRanksArray(final int sort1Size) {
	for (int j = 1, i = 0; i < sort1Size; i++, j++)
	    ranksSA[sort1[i] - 1] = j;
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
     * @return the maximum position that was filled + 1
     */
    private static final int moduleThree(int[] array, final int lo, final int factor, final int maxValue) {
	int index = lo;
	for (int value = factor; value < maxValue; value += CONS)
	    array[index++] = value;
	return index;
    }
}

/**
 * This class represents an int stack
 * 
 * @author jc.bojaca
 * @version 1
 */
class Stack {
    /** The first element of the stack */
    private Element apunt;

    /** The constructor */
    Stack() {
	apunt = null;
    }

    /** @return if the stack not have values. */
    boolean isEmpty() {
	return null == apunt;
    }

    /** @return the first element */
    int pop() {
	final int value = apunt.val;
	apunt = apunt.sig;
	return value;
    }

    /** Add one element into stack @param i add the element */
    void add(int i) {
	final Element toAdd = new Element(i);
	toAdd.sig = apunt;
	apunt = toAdd;
    }

    /**
     * Exchange this stack with the given stack
     * 
     * @param b
     *            the stack to exchange
     */
    void exchange(Stack b) {
	final Element temp = this.apunt;
	this.apunt = b.apunt;
	b.apunt = temp;
    }

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
	/** the next element */
	private Element sig;

	/** the value */
	private final int val;

	/** constructor */
	private Element() {
	    val = 0;
	}

	/** constructor @param val the value */
	private Element(int val) {
	    this.val = val;
	}

	@Override
	public String toString() {
	    return String.valueOf(val);
	}
    }
}
