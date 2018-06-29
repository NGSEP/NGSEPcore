package ngsep.assembly;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;

import ngsep.sequences.DNAMaskedSequence;

public class GraphMaker {
	static final int SEARCH_KMER_LENGTH = 15;
	static final int SEARCH_KMER_DISTANCE = 10;
	static final int SEARCH_KMER_OVERLAP = SEARCH_KMER_LENGTH - SEARCH_KMER_DISTANCE;
	static final int MAX_NUMBER_OF_NON_KMERS = 5;

	private final List<DNAMaskedSequence> sequences;
	private final FMIndexSequences index;

	private Map<Integer, List<ReadOverlap>> overlapsForward = new Hashtable<>();
	private Map<Integer, ReadOverlap> embeddedOverlaps = new Hashtable<>();

	private String[] kmers;
	private int lastKemrLength;
	private Queue<KmerAligment>[] procesing;
	private Queue<KmerAligment> hits;
	private Iterator<KmerAligment> iter;
	private KmerAligment kmerAligment;

	/**
	 * 
	 * @param sequences
	 * @param index
	 * 
	 *            <pre>
	 *  sequences is sorted by length
	 *            </pre>
	 */
	@SuppressWarnings("unchecked")
	public GraphMaker(List<DNAMaskedSequence> sequences, FMIndexSequences index) {
		this.sequences = sequences;
		this.index = index;
		kmers = new String[100 + ((sequences.get(0).length() - SEARCH_KMER_LENGTH + (SEARCH_KMER_DISTANCE - 1))
				/ SEARCH_KMER_DISTANCE)];
		procesing = new Queue[MAX_NUMBER_OF_NON_KMERS + 1];
		for (int i = 0; i < procesing.length; i++)
			procesing[i] = new LinkedList<>();
		hits = new PriorityQueue<>(new Comparator<KmerAligment>() {
			@Override
			public int compare(KmerAligment arg0, KmerAligment arg1) {
				return arg1.length() - arg0.length();
			}
		});
		createOverlapGraph();
	}

	/** @return the overlapsForward */
	public Map<Integer, List<ReadOverlap>> getOverlapsForward() {
		return overlapsForward;
	}

	/** @return the embeddedOverlaps */
	public Map<Integer, ReadOverlap> getEmbeddedOverlaps() {
		return embeddedOverlaps;
	}

	private void createOverlapGraph() {

		for (int idSequence = 0; idSequence < sequences.size(); idSequence++) {
			System.out.println("|| " + idSequence + "(" + sequences.size() + ") ||");
			long ini = System.nanoTime();
			if (embeddedOverlaps.containsKey(idSequence))
				continue;

			long sum = 0;
			int Kmers = kmers(sequences.get(idSequence));
			for (int idOtherSequence = idSequence + 1; idOtherSequence < sequences.size(); idOtherSequence++) {
				sum += calculateHits(Kmers, idOtherSequence);

				while (!hits.isEmpty())
					if (processOverlap(idSequence, idOtherSequence, hits.poll()))
						break;
			}
			long fin = System.nanoTime() - ini;
			System.out.println(fin + " , " + sum + " (" + (fin - sum) / (double) 1000000000 + ")");
		}
	}

	private long calculateHits(int kmersLength, int indexSequence) {
		long sum = 0;
		Set<Integer> search;
		hits.clear();
		for (int i = 0; i < kmersLength - 1; i++) {
			long ini = System.nanoTime();
			search = index.search(kmers[i], indexSequence, false);
			sum += System.nanoTime() - ini;
			Hits: for (Integer hit : search) {
				for (int j = 1; j < procesing.length; j++) {
					iter = procesing[j].iterator();
					while (iter.hasNext()) {
						kmerAligment = iter.next();
						if (kmerAligment.addAlignment(hit, SEARCH_KMER_LENGTH)) {
							addProcesing(kmerAligment);
							iter.remove();
							continue Hits;
						}
					}
				}
				addProcesing(i * SEARCH_KMER_DISTANCE, hit, SEARCH_KMER_LENGTH);
			}
			rotateProcesing();
		}
		long ini = System.nanoTime();
		search = index.search(kmers[kmersLength - 1], indexSequence, false);
		sum += System.nanoTime() - ini;
		Hits2: for (Integer hit : search) {
			for (Queue<KmerAligment> queue : procesing)
				for (KmerAligment kmerAligment2 : queue)
					if (kmerAligment2.addAlignment(hit, lastKemrLength))
						continue Hits2;
			addProcesing((kmersLength - 1) * SEARCH_KMER_DISTANCE, hit, lastKemrLength);
		}
		finalizeProcesing();
		return sum;
	}

	private void addProcesing(int ini1, int ini2, int lentgh) {
		addProcesing(new KmerAligment(ini1, ini2, ini2 + lentgh));
	}

	private void addProcesing(KmerAligment kmerAligment) {
		procesing[0].add(kmerAligment);
	}

	private void finalizeProcesing() {
		for (int i = 0; i < procesing.length; i++) {
			hits.addAll(procesing[i]);
			procesing[i].clear();
		}
	}

	private void rotateProcesing() {
		for (int i = 1; i < procesing.length; i++)
			for (KmerAligment keAligment : procesing[i])
				keAligment.addLast2(SEARCH_KMER_LENGTH);

		Queue<KmerAligment> aux = procesing[MAX_NUMBER_OF_NON_KMERS];
		hits.addAll(aux);
		aux.clear();
		System.arraycopy(procesing, 0, procesing, 1, MAX_NUMBER_OF_NON_KMERS);
		procesing[0] = aux;
	}

	private int kmers(CharSequence sequence) {
		String sequenceS = sequence.toString();
		sequenceS = sequenceS.toUpperCase();
		int j = 0;
		for (int i = 0; i <= sequence.length() - SEARCH_KMER_LENGTH; i += SEARCH_KMER_DISTANCE)
			kmers[j++] = sequenceS.substring(i, i + SEARCH_KMER_LENGTH);
		if (sequence.length() % SEARCH_KMER_LENGTH != 0) {
			kmers[j] = sequenceS.substring(
					((sequence.length() - SEARCH_KMER_OVERLAP) / SEARCH_KMER_DISTANCE) * SEARCH_KMER_DISTANCE);
			lastKemrLength = kmers[j].length();
			j++;
		} else
			lastKemrLength = SEARCH_KMER_LENGTH;
		return j;
	}

	private boolean processOverlap(int idSequence1, int idSequence2, KmerAligment overlap) {
		ReadOverlap overlapp = new ReadOverlap(idSequence1, overlap.getFirst1(), overlap.getLast1(), idSequence2,
				overlap.getFirst2(), overlap.getLast2(), false);
		int length1 = sequences.get(idSequence1).length();
		int length2 = sequences.get(idSequence2).length();

		int dl1 = overlap.getFirst1();
		int dl2 = overlap.getFirst2();
		int dr1 = length1 - overlap.getLast1();
		int dr2 = length2 - overlap.getLast2();
		int t = 10;
		int diffBorder = 100;

		// if (overlap.isNegativeStrand()) {
		// if (dl1 + t >= dr2 && dr1 + t >= dl2) {
		// if (dl2 > diffBorder || dr2 > diffBorder) {
		// // Not a true overlap
		// return;
		// }
		// if (!embeddedOverlaps.containsKey(overlap.getIndexSequence2())) {
		// ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(),
		// overlap.getFirst2(),
		// overlap.getLast2(), overlap.getIndexSequence1(), overlap.getFirst1(),
		// overlap.getLast1(),
		// true);
		// embeddedOverlaps.put(overlap.getIndexSequence2(), ov2);
		// }
		// return;
		// }
		// if (dl1 <= dr2 + t && dr1 <= dl2 + t) {
		// if (dl1 > diffBorder || dr1 > diffBorder) {
		// // Not a true overlap
		// return;
		// }
		// if (!embeddedOverlaps.containsKey(overlap.getIndexSequence1()))
		// embeddedOverlaps.put(overlap.getIndexSequence1(), overlap);
		// return;
		// }
		// if (dl1 < dr2 && dr1 > dl2) {
		// // 2 left of 1 negative strand
		// if (dl1 > diffBorder || dl2 > diffBorder) {
		// // Not a true overlap
		// return;
		// }
		// addOverlap(overlapsBackward, overlap);
		// ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(),
		// overlap.getFirst2(), overlap.getLast2(),
		// overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(),
		// true);
		// addOverlap(overlapsForward, ov2);
		// return;
		// }
		// if (dl1 > dr2 && dr1 < dl2) {
		// // 1 left of 2
		// if (dr2 > diffBorder || dr1 > diffBorder) {
		// // Not a true overlap
		// return;
		// }
		// addOverlap(overlapsForward, overlap);
		// ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(),
		// overlap.getFirst2(), overlap.getLast2(),
		// overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(),
		// true);
		// addOverlap(overlapsBackward, ov2);
		// return;
		// }
		// } else {
		if (dl1 + t >= dl2 && dr1 + t >= dr2) {
			if (dl2 > diffBorder || dr2 > diffBorder) {
				return false;
			}
			if (!embeddedOverlaps.containsKey(idSequence2)) {
				ReadOverlap ov2 = new ReadOverlap(idSequence2, overlap.getFirst2(), overlap.getLast2(), idSequence1,
						overlap.getFirst1(), overlap.getLast1(), false);
				embeddedOverlaps.put(idSequence2, ov2);
				return true;
			}
			return false;
		}
		if (dl1 <= dl2 + t && dr1 <= dr2 + t) {
			if (dl1 > diffBorder || dr1 > diffBorder) {
				// Not a true overlap
				return false;
			}
			if (!embeddedOverlaps.containsKey(idSequence1))
				embeddedOverlaps.put(idSequence1, overlapp);
			return true;
		}
		if (dl1 < dl2 && dr1 > dr2) {
			// 2 left of 1
			if (dl1 > diffBorder || dr2 > diffBorder) {
				// Not a true overlap
				return false;
			}
			ReadOverlap ov2 = new ReadOverlap(idSequence2, overlap.getFirst2(), overlap.getLast2(), idSequence1,
					overlap.getFirst1(), overlap.getLast1(), false);
			overlapsForward.computeIfAbsent(idSequence2, key -> new ArrayList<>()).add(ov2);
			return true;
		}
		if (dl1 > dl2 && dr1 < dr2) {
			// 1 left of 2
			if (dl2 > diffBorder || dr1 > diffBorder) {
				// Not a true overlap
				return false;
			}
			overlapsForward.computeIfAbsent(idSequence1, key -> new ArrayList<>()).add(overlapp);
			return true;
		}
		// }
		return false;
	}

	class KmerAligment {
		private final static int ERRORRATE = 5;
		private final int first1;
		private final int first2;
		private int last2;

		public KmerAligment(int first1, int first2, int last2) {
			super();
			this.first1 = first1;
			this.first2 = first2;
			this.last2 = last2;
		}

		public boolean addAlignment(int hit, int kmerSize) {
			int dif = Math.abs(last2 - OverlapGraph.SEARCH_KMER_OVERLAP - hit);
			if (dif < ERRORRATE) {
				last2 = hit + kmerSize;
				return true;
			}
			return false;
		}

		public void addLast2(int value) {
			last2 += value;
		}

		public int length() {
			return last2 - first2;
		}

		/** @return the first1 */
		public int getFirst1() {
			return first1;
		}

		/** @return the last2 */
		public int getLast1() {
			return first1 + length();
		}

		/** @return the first2 */
		public int getFirst2() {
			return first2;
		}

		/** @return the last2 */
		public int getLast2() {
			return last2;
		}
	}
}
