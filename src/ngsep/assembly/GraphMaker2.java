package ngsep.assembly;

import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.LimitedSequence;

public class GraphMaker2 {
	static final int SEARCH_KMER_LENGTH = 21;

	private final int sequencesSize;
	private final List<DNAMaskedSequence> sequences;
	private Map<Integer, ReadOverlap> embeddedOverlaps = new Hashtable<>();

	public GraphMaker2(List<DNAMaskedSequence> sequences) {
		this.sequences = sequences;
		this.sequencesSize = sequences.size();
		Collections.sort(sequences, (LimitedSequence l1, LimitedSequence l2) -> l1.length() - l2.length());
		emmbeddeDetection(sequences);
	}

	public void emmbeddeDetection(List<DNAMaskedSequence> sequencesToIndex) {
		/** create the FMIndex **/
		System.out.println("start");
		long i = System.currentTimeMillis();
		FMIndex index = new FMIndex();
		index.loadUnnamedSequences("", sequencesToIndex, 100, 25);
		System.out.println("index: " + (System.currentTimeMillis() - i));
		i = System.currentTimeMillis();
		/** the arrays with the overlap **/
		@SuppressWarnings("unchecked")
		TreeMap<Integer, int[]>[][] a = new TreeMap[KmerEmmbeddedIterator.NUMBER_OF_KMERS + 1][sequencesToIndex.size()];
		for (int y = 0; y <= KmerEmmbeddedIterator.NUMBER_OF_KMERS; y++)
			for (int x = 0; x < sequencesToIndex.size(); x++)
				a[y][x] = new TreeMap<>();

		for (int idSequence = 0; idSequence < sequencesSize; idSequence++) {
			System.out.println(idSequence + "->" + sequencesToIndex.size());
			for (TreeMap<Integer, int[]>[] b : a)
				for (TreeMap<Integer, int[]> c : b)
					c.clear();

			KmerEmmbeddedIterator ki = new KmerEmmbeddedIterator(sequences.get(idSequence));
			int k = 0;
			for (String string : ki.firts()) {
				for (ReadAlignment aln : index.search(string)) {
					align(aln, a, k * KmerEmmbeddedIterator.SEARCH_KMER_LENGTH);
				}
				k++;
				rotate(a);
			}
			for (String string : ki.lasts()) {
				index.search(string);
				k++;
				rotate(a);
			}

			System.out.println(k);
		}
		System.out.println("asdljkhjfgdfdkhjgjfdrse: " + (System.currentTimeMillis() - i));
	}

	private void align(ReadAlignment aln, TreeMap<Integer, int[]>[][] a, int posinRefrenece) {
		int idSequenceAligned = Integer.parseInt(aln.getSequenceName().substring(1));
		int pos = aln.getFirst();

		int minposkmer = 0, minnode = 0, diff = 0, min = 6;
		for (int dif = KmerEmmbeddedIterator.SEARCH_KMER_LENGTH, i = 1; i <= KmerEmmbeddedIterator.NUMBER_OF_KMERS; i++, dif += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH) {
			TreeMap<Integer, int[]> treeMap = a[i][idSequenceAligned];
			if (treeMap.isEmpty())
				continue;

			int t = pos - dif;
			Integer p1 = treeMap.ceilingKey(t);
			Integer p2 = treeMap.floorKey(t);
			if (p1 != null && Math.abs(p1 - t) < min) {
				diff = p1;
				min = Math.abs(p1 - t);
				minposkmer = i;
				minnode = p1;
			}
			if (p2 != null && Math.abs(p2 - t) < min) {
				diff = p2;
				min = Math.abs(p2 - t);
				minposkmer = i;
				minnode = p2;
			}
		}

		if (min == 6) {
			a[0][idSequenceAligned].put(pos, new int[] { pos, posinRefrenece, 1 });
		} else {
			int[] aux = a[minposkmer][idSequenceAligned].remove(minnode);
			aux[2]++;
			a[0][idSequenceAligned].put(diff, aux);
		}
	}

	private void rotate(TreeMap<Integer, int[]>[][] a) {
		TreeMap<Integer, int[]>[] temp = a[KmerEmmbeddedIterator.NUMBER_OF_KMERS];
		for (TreeMap<Integer, int[]> t : temp)
			t.clear();
		System.arraycopy(a, 0, a, 1, KmerEmmbeddedIterator.NUMBER_OF_KMERS);
		a[0] = temp;
	}
}
