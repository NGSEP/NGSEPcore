package ngsep.assembly;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;

public class OverlapingDetector1 implements OverlappingDetector {
	static final int SEARCH_KMER_LENGTH = 15;
	static final int SEARCH_KMER_DISTANCE = 10;
	static final int SEARCH_KMER_OVERLAP = SEARCH_KMER_LENGTH
			- SEARCH_KMER_DISTANCE;
	private FMIndex index;
	private Map<Integer, Edge> ignored;
	private List<DNAMaskedSequence> sequences;
	Map<Integer, List<Edge>> edges = new HashMap<>();

	@Override
	public Map<Integer, List<Edge>> getEdges(FMIndex index,
			List<DNAMaskedSequence> sequences, Map<Integer, Edge> ignored) {
		this.index = index;
		this.sequences = sequences;
		this.ignored = ignored;
		return null;
	}

	@SuppressWarnings("unchecked")
	public void findOverlaps2() {
		List<ReadOverlap2>[] procesingHits = new List[sequences.size()];
		Queue<ReadOverlap2>[] hits = new Queue[sequences.size()];

		for (int idSequence = 0; idSequence < sequences.size(); idSequence++) {
			System.out
					.println("  --> " + idSequence + " - " + sequences.size());
			long ini = System.nanoTime();
			if (ignored.containsKey(idSequence))
				continue;
			long time = findHits(idSequence, procesingHits, hits);

			for (int i = 0; i < hits.length; i++) {
				Queue<ReadOverlap2> queue = hits[i];
				while (!queue.isEmpty())
					if (processOverlap(idSequence, i, queue.poll()))
						break;

			}
			long fin = System.nanoTime() - ini;
			System.out.println(fin + " " + time + " (" + (fin - time)
					/ (double) 1000000000 + ")");
		}
	}

	private long findHits(int idSequence, List<ReadOverlap2>[] procesingHits,
			Queue<ReadOverlap2>[] hits) {
		for (int i = 0; i < procesingHits.length; i++)
			// TODO optimizar (array), limpiar colas, cambio en las estructuras
			// para
			// soportar eso
			procesingHits[i] = new LinkedList<>();
		for (int i = 0; i < hits.length; i++)
			hits[i] = new PriorityQueue<>(new Comparator<ReadOverlap2>() {

				@Override
				public int compare(ReadOverlap2 arg0, ReadOverlap2 arg1) {
					return arg1.length() - arg0.length();
				}
			});

		long sum = 0;
		String sequence = sequences.get(idSequence).toString();
		for (int i = 0; i <= sequence.length() - SEARCH_KMER_LENGTH; i += SEARCH_KMER_DISTANCE) {
			InitializeHits(procesingHits);
			sum += processKmer(idSequence, i,
					sequence.substring(i, i + SEARCH_KMER_LENGTH),
					procesingHits);
			finalizeHits(procesingHits, hits);
		}
		if (sequence.length() % SEARCH_KMER_LENGTH != 0) {
			int KmerStarPosition = ((sequence.length() - SEARCH_KMER_OVERLAP) / SEARCH_KMER_DISTANCE)
					* SEARCH_KMER_DISTANCE;
			InitializeHits(procesingHits);
			sum += processKmer(idSequence, KmerStarPosition,
					sequence.substring(KmerStarPosition), procesingHits);
			finalizeHits(procesingHits, hits);
		}
		for (int i = 0; i < procesingHits.length; i++)
			for (ReadOverlap2 read : procesingHits[i])
				hits[i].add(read);
		return sum;

	}

	private void InitializeHits(List<ReadOverlap2>[] procesingHits) {
		for (List<ReadOverlap2> set : procesingHits)
			for (ReadOverlap2 read : set)
				read.nonAdded();
	}

	private long processKmer(int idSequence, int KmerStarPosition, String kmer,
			List<ReadOverlap2>[] procesingHits) {
		long ini = System.nanoTime();
		Iterable<ReadAlignment> itre = index.search(kmer);
		ini = System.nanoTime() - ini;
		for (ReadAlignment aln : itre) {
			int idSequenceAligned = Integer.parseInt(aln.getSequenceName());
			if (idSequenceAligned > idSequence) {
				if (!aln.isNegativeStrand()) {
					boolean alienado = false;
					for (ReadOverlap2 read : procesingHits[idSequenceAligned]) {
						if (read.addAlignment(aln, kmer.length())) {
							alienado = true;
							break;
						}
					}
					if (!alienado)
						procesingHits[idSequenceAligned].add(new ReadOverlap2(
								KmerStarPosition, aln.getFirst() - 1, aln
										.getLast()));
				} else {
					// TODO caso negativo
					System.out.println("ouchh");
					continue;
				}
			}
		}
		return ini;
	}

	private void finalizeHits(List<ReadOverlap2>[] procesingHits,
			Queue<ReadOverlap2>[] hits) {
		for (int i = 0; i < procesingHits.length; i++) {
			Iterator<ReadOverlap2> iter = procesingHits[i].iterator();
			while (iter.hasNext()) {
				ReadOverlap2 read = iter.next();
				if (!read.added()) {
					read.errorCountPlusPlus();
					if (read.getErrorCount() >= 5) {
						iter.remove();
						if (read.length() / (double) sequences.get(i).length() >= 0.15)
							hits[i].add(read);
					}
				}
			}
		}
	}

	private boolean processOverlap(int idSequence1, int idSequence2,
			ReadOverlap2 overlap) {
		Edge overlapp = new Edge(idSequence1, idSequence2, overlap.getFirst1(),
				overlap.getLast1(), overlap.getFirst2(), overlap.getLast2());
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
			if (!ignored.containsKey(idSequence2)) {
				Edge ov2 = new Edge(idSequence2, idSequence1,
						overlap.getFirst2(), overlap.getLast2(),
						overlap.getFirst1(), overlap.getLast1());
				ignored.put(idSequence2, ov2);
				return true;
			}
			return false;
		}
		if (dl1 <= dl2 + t && dr1 <= dr2 + t) {
			if (dl1 > diffBorder || dr1 > diffBorder) {
				// Not a true overlap
				return false;
			}
			if (!ignored.containsKey(idSequence1))
				ignored.put(idSequence1, overlapp);
			return true;
		}
		if (dl1 < dl2 && dr1 > dr2) {
			// 2 left of 1
			if (dl1 > diffBorder || dr2 > diffBorder) {
				// Not a true overlap
				return false;
			}
			Edge ov2 = new Edge(idSequence2, idSequence1, overlap.getFirst2(),
					overlap.getLast2(), overlap.getFirst1(), overlap.getLast1());
			addOverlap(ov2);
			return true;
		}
		if (dl1 > dl2 && dr1 < dr2) {
			// 1 left of 2
			if (dl2 > diffBorder || dr1 > diffBorder) {
				// Not a true overlap
				return false;
			}
			addOverlap(overlapp);
			return true;
		}
		// }
		return false;
	}

	private void addOverlap(Edge overlap) {
		List<Edge> overlapsSeq1 = edges.get(overlap.getIndexSequence1());
		if (overlapsSeq1 == null) {
			overlapsSeq1 = new ArrayList<>();
			edges.put(overlap.getIndexSequence1(), overlapsSeq1);
		}
		overlapsSeq1.add(overlap);
	}
}

class ReadOverlap2 {
	private final static int ERRORRATE = 5;

	private final int first1;
	private final int first2;
	private int last1;
	private int last2;
	private int errorCount;
	private boolean added;

	public ReadOverlap2(int first1, int first2, int last2) {
		super();
		this.first1 = first1;
		this.first2 = first2;
		this.last2 = last2;
		this.last1 = first1 + OverlapingDetector1.SEARCH_KMER_LENGTH;
		errorCount = 0;
		added = true;
	}

	public boolean addAlignment(ReadAlignment aln, int kmerSize) {
		if (added())
			return false;
		int dif = Difenrece(last2 - OverlapingDetector1.SEARCH_KMER_OVERLAP,
				aln.getFirst());
		if (dif < ERRORRATE) {
			added = true;
			last2 = aln.getLast();
			last1 += OverlapingDetector1.SEARCH_KMER_DISTANCE * errorCount;
			last1 += kmerSize - OverlapingDetector1.SEARCH_KMER_OVERLAP;
			errorCount = 0;
			return true;
		} else
			return false;
	}

	public boolean added() {
		return added;
	}

	private final static int Difenrece(int a, int b) {
		return Math.abs(a - b);
	}

	/**
	 * @return the errorCount
	 */
	public int getErrorCount() {
		return errorCount;
	}

	/**
	 * @param errorCount
	 *            the errorCount to set
	 */
	public void errorCountPlusPlus() {
		this.errorCount++;
	}

	/**
	 * @param diference
	 *            the diference to set
	 */
	public void nonAdded() {
		this.added = false;
	}

	/** @return the first1 */
	public int getFirst1() {
		return first1;
	}

	/** @return the first2 */
	public int getFirst2() {
		return first2;
	}

	/**
	 * @return the last1
	 */
	public int getLast1() {
		return last1;
	}

	/**
	 * @return the last2
	 */
	public int getLast2() {
		return last2;
	}

	public int length() {
		return last1 - first1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "ReadOverlap2 [first1=" + first1 + ", first2=" + first2
				+ ", last1=" + last1 + ", last2=" + last2 + "]";
	}

}
