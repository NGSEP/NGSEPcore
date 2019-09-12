package ngsep.sequences;

import static java.lang.Math.ceil;
import static java.lang.Math.log;
import static java.lang.Math.max;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class FMIndexSingleSequence2 {
	// -----------Constants----------------------------------------------
	private final static double LN_2 = log(2);
	private static final int DEFAULT_SUFFIX_FRACTION = 16;
	private static final int DEFAULT_TALLY_DISTANCE = 2;// *64

	// -----------About the sequence-------------------------------------
	private String alphabet;
	private int size;

	// -----------About the Index-----------------------------------------
	private Map<Integer, Integer> sa = new TreeMap<>();
	private int[][] tally;
	private long[][] bwt;
	private int[] count;
	private byte[] map;
	private int[] rsa;

	// ------------About the wrapper--------------------------------------
	private int bitsPerLetter;
	private int tallyFraction;
	private int suffixFraction;

	// ------------To avoid the extra letter into BWT --------------------
	private int bwtExceptIndex, bwtExceptLongIndex;

	public FMIndexSingleSequence2(final CharSequence sequence) {
		this(sequence, DEFAULT_TALLY_DISTANCE, DEFAULT_SUFFIX_FRACTION);
	}

	public FMIndexSingleSequence2(CharSequence sequence, int tallyFraction, int suffixFraction) {
		this.tallyFraction = tallyFraction;
		this.suffixFraction = suffixFraction;
		calculate(sequence);
	}

	/**
	 * Searches the given sequence in this FMIndex
	 * 
	 * @param searchSequence Sequence to search
	 * @return Set<Integer> Set of start positions for the given sequence
	 */
	public Set<Integer> search(CharSequence seq) {
		char lastChar = seq.charAt(seq.length() - 1);
		int lo = count[lastChar];
		int hi = count[lastChar + 1] - 1;

		for (int k = seq.length() - 2; k >= 0 && lo <= hi; k--) {
			char c = seq.charAt(k);
			int cInd = count[c];
			lo = cInd + cumletters(lo - 1, map[c]);
			hi = cInd + cumletters(hi, map[c]) - 1;
		}

		Set<Integer> ans = new HashSet<>();
		for (int i = lo; i <= hi; i++) {
			int ind = i;
			int c = 0;
			Integer val = sa.get(ind);
			while (val == null) {
				byte cha = bwtCharAt(ind);
				ind = count[alphabet.charAt(cha)] + cumletters(ind, cha) - 1;
				val = sa.get(ind);
				c++;
			}
			ans.add(val + c);
		}
		return ans;
	}

	/**
	 * Return the subsequence of the indexed sequence between the given coordinates
	 * 
	 * @param start position of the sequence (0-based, included)
	 * @param end   position of the sequence (0-based, excluded)
	 * @return CharSequence segment of the indexed sequence between the given
	 *         coordinates
	 */
	public CharSequence getSequence(int start, int end) {
		if (start >= end)
			throw new StringIndexOutOfBoundsException(
					"Start position " + start + " should be smaller than end position: " + end);

		int ind = (end + suffixFraction - 1) / suffixFraction;
		int posSA = rsa[ind];
		if (ind == rsa.length - 1) {
			ind = posSA;
			posSA = 0;
		} else
			ind *= suffixFraction;

		while (ind > end) {
			byte cha = bwtCharAt(posSA);
			posSA = count[alphabet.charAt(cha)] + cumletters(posSA, cha) - 1;
			ind--;
		}

		char[] seq = new char[end - start];
		for (int i = end - start - 1; i >= 0; i--) {
			byte cha = bwtCharAt(posSA);
			char character = alphabet.charAt(cha);
			seq[i] = character;
			posSA = count[character] + cumletters(posSA, cha) - 1;
		}
		return new String(seq);
	}

	private void calculate(CharSequence sequence) {
		size = sequence.length();
		count = countLetters(sequence);

		alphabet = getAlphabet();
		map = getMap();
		byte[] data = transform(sequence);

		countsCumsum();
		reMap();

		bitsPerLetter = max(1, (int) log2(alphabet.length()));
		SuffixArrayGenerator2 suffixArrayGenerator = new SuffixArrayGenerator2(data);
		int[] suffixArray = suffixArrayGenerator.getSuffixArray();
		int[] reverseSA = suffixArrayGenerator.getReverseSuffixArray();

		bwt = getBWT(suffixArray, data, bitsPerLetter);
		tally = buildTally();

		partialSufixes(suffixArray, reverseSA);
	}

	private void partialSufixes(int[] suffixArray, int[] reverseSA) {
		sa = new HashMap<>();
		rsa = new int[1 + (reverseSA.length + suffixFraction - 1) / suffixFraction];
		sa.put(0, suffixArray[0]);
		for (int i = 0, j = 0; i < reverseSA.length; i += suffixFraction, j++) {
			sa.put(reverseSA[i], i);
			rsa[j] = reverseSA[i];
		}
		rsa[rsa.length - 1] = reverseSA.length - 1;
	}

	private void countsCumsum() {
		count[0] = 1;
		int cum = 1;
		int p = 0;
		for (int i = 1; i < count.length; i++) {
			p = count[i];
			count[i] = cum;
			cum += p;
		}
	}

	private int cumletters(int ind, byte let) {
		int longInd = ind >>> 6;
		int letterInd = ind & 0x3f;

		int tallyInd = longInd / tallyFraction;
		int countByte = tally[tallyInd][let];

		for (int i = (longInd / tallyFraction) * tallyFraction; i < longInd; i++)
			countByte += Long.bitCount(trans(bwt[i], let));

		int a = 64 - (letterInd + 1);
		countByte += Long.bitCount((trans(bwt[longInd], let) << a) >>> a);

		// remove the extra zero letter of Exception $
		if (let == 0 && ind >= bwtExceptIndex && longInd >= bwtExceptLongIndex
				&& (longInd / tallyFraction) * tallyFraction <= bwtExceptLongIndex)
			countByte--;

		return countByte;
	}

	private byte bwtCharAt(int i) {
		long[] longs = bwt[i >>> 6];
		int j = i & 0x3f;
		byte ind = 0;
		for (int bit = 0; bit < bitsPerLetter; bit++)
			ind |= ((longs[bit] >>> j) & 1) << bit;
		return ind;
	}

	private int[][] buildTally() {
		int tallyRows = 1 + (bwt.length + tallyFraction - 1) / tallyFraction;
		final int[] arr = new int[alphabet.length()];

		int[][] tally = new int[tallyRows][alphabet.length()];

		int i = 0, j = 1;
		while (i < ((bwt.length - 1) / tallyFraction) * tallyFraction) {
			for (int lll = 0; lll < tallyFraction; lll++) {
				for (int k = 0; k < arr.length; k++)
					arr[k] += Long.bitCount(trans(bwt[i], (byte) k));
				if (i == bwtExceptLongIndex)
					arr[0]--;
				i++;
			}
			System.arraycopy(arr, 0, tally[j++], 0, arr.length);
		}

		while (i < bwt.length - 1) {
			for (int k = 0; k < arr.length; k++)
				arr[k] += Long.bitCount(trans(bwt[i], (byte) k));
			if (i == bwtExceptLongIndex)
				arr[0]--;
			i++;
		}
		int a = 64 - (size + 1 - (i << 6));
		for (int k = 0; k < arr.length; k++)
			arr[(int) k] += Long.bitCount((trans(bwt[i], (byte) k) << a) >>> a);
		if (i == bwtExceptLongIndex)
			arr[0]--;

		System.arraycopy(arr, 0, tally[j++], 0, arr.length);
		return tally;
	}

	private long[][] getBWT(int[] suffixArray, byte[] data, int bitsPerLetter) {
		long[][] BWT = new long[(suffixArray.length + 63) >>> 6][bitsPerLetter];
		int roundLength = suffixArray.length & ~0x3f;

		int j = 0, i = 0;
		long a;
		while (i < roundLength) {
			long[] longs = BWT[j];
			for (int k = 0; k < 64; k++) {
				if (suffixArray[i] > 0)
					a = data[suffixArray[i] - 1] - 1;
				else {
					bwtExceptIndex = i;
					a = 0;
				}
				for (int bit = 0; bit < bitsPerLetter; bit++)
					longs[bit] |= (1 & (a >>> bit)) << k;
				i++;
			}
			j++;
		}
		long[] longs = BWT[j];
		for (int k = 0; i < suffixArray.length; k++) {
			if (suffixArray[i] > 0)
				a = data[suffixArray[i] - 1] - 1;
			else {
				bwtExceptIndex = i;
				a = 0;
			}
			for (int bit = 0; bit < bitsPerLetter; bit++)
				longs[bit] |= (1 & (a >>> bit)) << k;
			i++;
		}
		bwtExceptLongIndex = bwtExceptIndex >>> 6;
		return BWT;
	}

	private byte[] transform(CharSequence sequence) {
		final int size = sequence.length();
		final byte[] data = new byte[size + 1];
		for (int i = 0; i < size; i++)
			data[i] = map[sequence.charAt(i)];
		return data;
	}

	private byte[] getMap() {
		final byte[] map = new byte[Byte.MAX_VALUE + 2];
		byte k = 1;
		for (int i = 0; i < alphabet.length(); i++)
			map[alphabet.charAt(i)] = k++;
		return map;
	}

	private void reMap() {
		for (int i = 0; i < map.length; i++)
			map[i]--;
	}

	private static int[] countLetters(CharSequence sequence) {
		final int[] count = new int[Byte.MAX_VALUE + 2];
		final int size = sequence.length();
		for (int i = 0; i < size; i++) {
			byte c = (byte) sequence.charAt(i);
			count[c]++;
		}
		return count;
	}

	private String getAlphabet() {
		final StringBuilder alphabetBuilder = new StringBuilder();
		for (int i = 0; i < count.length; i++) {
			if (count[i] != 0)
				alphabetBuilder.append((char) i);
		}
		return alphabetBuilder.toString();
	}

	private long trans(long[] ind, byte m) {
		long ans = -1;
		for (int i = 0; i < bitsPerLetter; i++)
			ans &= (~-(1 & m >>> i)) ^ ind[i];

		return ans;
	}

	private final static double log2(int x) {
		return ceil(log(x) / LN_2);
	}
}