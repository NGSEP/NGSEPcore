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

/**
 * 
 * this is a test of generating the SA with the DC3 algorithm
 * 
 * @author Juan Camilo Bojaca
 *
 */
public class SuffixArrayGenerator {
	private final static int Cte_Radix = 4; // == log2(16)
	private final static int Cte_Radix_Bit = (int) Math.pow(2, Cte_Radix) - 1;
	private final static int AlfabetStart = 34;
	private final static int AlfabetEnd = 127;
	private int[] ans;

	public SuffixArrayGenerator(CharSequence charSequence) {
		int[] data = new int[charSequence.length() + 1];
		byte[] map = new byte[AlfabetEnd - AlfabetStart + 1];
		byte character;

		for (int i = 0; i < charSequence.length(); i++) {
			character = (byte) (charSequence.charAt(i) - AlfabetStart);
			if (map[character] == 0)
				map[character]++;
		}

		for (int i = 1; i < map.length; i++)
			map[i] += map[i - 1];

		data = new int[charSequence.length() + 1];
		for (int i = 0; i < charSequence.length(); i++)
			data[i] = map[charSequence.charAt(i) - AlfabetStart];

		ans = getSuffix(data, (highestOneBitPos(map[map.length - 1]) / Cte_Radix) * Cte_Radix);
	}

	private int[] getSuffix(final int[] data, final int MaxBit) {
		// variables
		int[] ans, sort, sort2, SA;
		// aux
		int[] contSort, auxSort;
		Stack loStack, hiStack, loStackAux = new Stack(), hiStackAux = new Stack();
		int hi, hiAux, lo, loAux, Ci, Bi, d, bit, maxBitSA, ib, ic, k, idx;
		boolean[] ranks;
		// -----------------------------------------------------------------------------
		// STEP #1: create C
		// -----------------------------------------------------------------------------
		sort = new int[(int) ((2 * data.length) / 3)];
		Ci = 0;
		for (int i = 1; i < data.length; i += 3)
			sort[Ci++] = i;
		int m = Ci;
		for (int i = 2; i < data.length; i += 3)
			sort[Ci++] = i;

		// -----------------------------------------------------------------------------
		// STEP #2: radix sort C
		// -----------------------------------------------------------------------------
		contSort = new int[Cte_Radix_Bit + 3];
		auxSort = new int[sort.length];

		loStack = new Stack(0);
		hiStack = new Stack(Ci - 1);

		d = 0;
		bit = MaxBit;
		while (!loStack.isEmpty() && d < 3) {

			while (!loStack.isEmpty()) {
				lo = loStack.pop();
				hi = hiStack.pop();
				Arrays.fill(contSort, 0);

				for (int i = lo; i <= hi; i++)
					contSort[((data[sort[i] + d] >>> bit) & Cte_Radix_Bit) + 2]++;

				for (int i = 3; i < contSort.length; i++)
					contSort[i] += contSort[i - 1];

				for (int i = lo; i <= hi; i++)
					auxSort[contSort[((data[sort[i] + d] >>> bit) & Cte_Radix_Bit) + 1]++] = sort[i];

				for (int i = lo; i <= hi; i++)
					sort[i] = auxSort[i - lo];

				for (int act, prev = contSort[0], i = 1; prev != hi - lo + 1; prev = act, i++) {
					act = contSort[i];
					loAux = prev + lo;
					hiAux = act - 1 + lo;
					if (loAux < hiAux) {
						loStackAux.add(loAux);
						hiStackAux.add(hiAux);
					}
				}
			}
			Stack.exchange(loStack, loStackAux);
			Stack.exchange(hiStack, hiStackAux);

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
			ranks = new boolean[sort.length];

			while (!loStack.isEmpty()) {
				lo = loStack.pop();
				hi = hiStack.pop();
				for (int i = lo + 1; i <= hi; i++)
					ranks[i] = true;
			}

			sort2 = new int[Ci + 1];
			d = 0;
			for (int i = 0; i < sort.length; ++i) {
				idx = sort[i];
				if (!ranks[i])
					++d;
				sort2[idx % 3 == 1 ? idx / 3 : idx / 3 + m] = d;
			}

			int[] SApr = getSuffix(sort2, (highestOneBitPos(d) / Cte_Radix) * Cte_Radix);

			for (int j = 1; j < SApr.length; j++) {
				idx = SApr[j];
				sort[j - 1] = (idx < m) ? idx * 3 + 1 : (idx - m) * 3 + 2;
			}
		}

		SA = new int[data.length];
		for (int j = 1, i = 0; i < sort.length; i++, j++)
			SA[sort[i]] = j;
		maxBitSA = (highestOneBitPos(sort.length) / Cte_Radix) * Cte_Radix;

		// -----------------------------------------------------------------------------
		// STEP #3: Sorting the non sample Suffices b0
		// -----------------------------------------------------------------------------
		sort2 = new int[((data.length - 1) / 3) + 1];

		Bi = 0;
		for (int i = 0; i < data.length; i += 3)
			sort2[Bi++] = i;

		loStack = new Stack(0);
		hiStack = new Stack(Bi - 1);

		bit = MaxBit;
		while (loStack != null && bit >= 0) {
			while (!loStack.isEmpty()) {
				lo = loStack.pop();
				hi = hiStack.pop();
				Arrays.fill(contSort, 0);

				for (int i = lo; i <= hi; i++)
					contSort[((data[sort2[i]] >>> bit) & Cte_Radix_Bit) + 2]++;

				for (int i = 3; i < contSort.length; i++)
					contSort[i] += contSort[i - 1];

				for (int i = lo; i <= hi; i++)
					auxSort[contSort[((data[sort2[i]] >>> bit) & Cte_Radix_Bit) + 1]++] = sort2[i];

				for (int i = lo; i <= hi; i++)
					sort2[i] = auxSort[i - lo];

				for (int i = 1; contSort[i - 1] != hi - lo + 1; i++) {
					loAux = contSort[i - 1] + lo;
					hiAux = contSort[i] - 1 + lo;
					if (loAux < hiAux) {
						loStackAux.add(loAux);
						hiStackAux.add(hiAux);
					}
				}
			}
			Stack.exchange(loStack, loStackAux);
			Stack.exchange(hiStack, hiStackAux);
			bit -= Cte_Radix;
		}
		bit = maxBitSA;
		while (loStack != null && bit >= 0) {
			while (!loStack.isEmpty()) {
				lo = loStack.pop();
				hi = hiStack.pop();
				Arrays.fill(contSort, 0);

				for (int i = lo; i <= hi; i++)
					contSort[((SA[sort2[i] + 1] >>> bit) & Cte_Radix_Bit) + 2]++;

				for (int i = 3; i < contSort.length; i++)
					contSort[i] += contSort[i - 1];

				for (int i = lo; i <= hi; i++)
					auxSort[contSort[((SA[sort2[i] + 1] >>> bit) & Cte_Radix_Bit) + 1]++] = sort2[i];

				for (int i = lo; i <= hi; i++)
					sort2[i] = auxSort[i - lo];

				for (int i = 1; contSort[i - 1] != hi - lo + 1; i++) {
					loAux = contSort[i - 1] + lo;
					hiAux = contSort[i] - 1 + lo;
					if (loAux < hiAux) {
						loStackAux.add(loAux);
						hiStackAux.add(hiAux);
					}
				}
			}

			Stack.exchange(loStack, loStackAux);
			Stack.exchange(hiStack, hiStackAux);

			bit -= Cte_Radix;
		}
		// -----------------------------------------------------------------------------
		// STEP #4: Merge
		// -----------------------------------------------------------------------------
		ans = new int[Bi + Ci];
		ib = 0;
		ic = 0;
		k = 0;
		while (ib != Ci && ic != Bi) {
			int v, oi = sort[ib], oj = sort2[ic];
			if (oi % 3 == 1) {
				v = data[oj] - data[oi];
				if (v == 0)
					v = SA[oj + 1] - SA[oi + 1];
			} else {
				v = data[oj] - data[oi];
				if (v == 0)
					v = data[oj + 1] - data[oi + 1];
				if (v == 0)
					v = SA[oj + 2] - SA[oi + 2];
			}
			if (v < 0)
				ans[k++] = sort2[ic++];
			else
				ans[k++] = sort[ib++];
		}
		while (ib != Ci)
			ans[k++] = sort[ib++];
		while (ic != Bi)
			ans[k++] = sort2[ic++];

		return ans;
	}

	private static int highestOneBitPos(int i) {
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

	public int [] getSA() {
		int [] answer = new int [ans.length-1];
		for(int i=1;i<ans.length;i++) {
			answer[i-1] = ans[i];
		}
		return answer;
	}

	/**
	 * This class represents an integer stack
	 * 
	 * @author jc.bojaca
	 */
	private static class Stack {
		private Int apunt;

		public Stack() {
			apunt = null;
		}

		public Stack(int i) {
			apunt = null;
			Int temp = new Int(i);
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
			Int temp = new Int(i);
			temp.sig = apunt;
			apunt = temp;
		}

		public static void exchange(Stack a, Stack b) {
			Int temp = a.apunt;
			a.apunt = b.apunt;
			b.apunt = temp;
		}

		private static class Int {
			Int sig;
			int val;

			public Int(int val) {
				this.val = val;
			}
		}
	}
}
