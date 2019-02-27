package ngsep.assembly;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.FMIndex;

public class GraphBuilder_OverlapFinder_Quee implements
		GraphBuilder_OverlapFinder {
	static final int SEARCH_KMER_LENGTH = 10;
	static final int SEARCH_KMER_DISTANCE = 10;
	static final int SEARCH_KMER_OVERLAP = SEARCH_KMER_LENGTH - SEARCH_KMER_DISTANCE;

	private List<CharSequence> sequences;
	private AssemblyGraph assemblyGraph;
	private FMIndex index;

	private Map<Integer, List<ReadOverlap>> overlapsForward = new Hashtable<>();
	private Map<Integer, List<ReadOverlap>> overlapsBackward = new Hashtable<>();
	private Map<Integer, ReadOverlap> embeddedOverlaps = new Hashtable<>();
	
	@Override
	public void calculate(List<CharSequence> seq, FMIndex index) {
		this.index = index;
		this.sequences = seq;
		findOverlaps();
		assemblyGraph = buidAssemblyGraph();
	}

	@Override
	public AssemblyGraph getGrap() {
		return assemblyGraph;
	}

	
	public void findOverlaps() {
		for (int s = 0; s < sequences.size(); s++) {
			if (embeddedOverlaps.containsKey(s))
				continue;

			Map<Integer, List<ReadAlignment>> seqHits = findHitsRead(s);
			buildOverlapsFromKmerAlignments(s, seqHits);
		}
	}
	/**
	 * find the hints for CharSequence
	 * 
	 * @param idRead
	 * @param read
	 * @return a map whit all the hits per sequence.
	 */
	private Map<Integer, List<ReadAlignment>> findHitsRead(int idRead) {
		Map<Integer, List<ReadAlignment>> seqHits = new HashMap<>(sequences.size());
		String read = sequences.get(idRead).toString();

		for (int i = 0; i < read.length() - SEARCH_KMER_LENGTH; i += SEARCH_KMER_LENGTH) {
			String kmer = read.substring(i, i + SEARCH_KMER_LENGTH);

			List<ReadAlignment> regions = index.search(kmer);
			for (ReadAlignment aln : regions) {
				// Use mate start to store the kmer start site producing the hit
				aln.setMateFirst(i);
				int k = Integer.parseInt(aln.getSequenceName());
				if (k > idRead) {
					seqHits.computeIfAbsent(k, key -> new ArrayList<>()).add(aln);
				}
			}
		}
		return seqHits;
	}
	
	private void buildOverlapsFromKmerAlignments(int searchId, Map<Integer, List<ReadAlignment>> seqHits) {
		CharSequence searchSequence = sequences.get(searchId);
		for (Entry<Integer, List<ReadAlignment>> entry : seqHits.entrySet()) {
			int k = entry.getKey();
			List<ReadAlignment> alns = entry.getValue();

			ReadOverlap next = null;
			for (ReadAlignment aln : alns) {
				// System.out.println(next == null );
				if (next == null || !next.addKmerAlignment(k, aln)) {
					if (next != null)
						processOverlap(next, searchSequence.length(), sequences.get(k).length());
					next = new ReadOverlap(searchId, aln.getMateFirst(), aln.getMateFirst() + SEARCH_KMER_LENGTH - 1, k,
							aln.getFirst(), aln.getLast(), aln.isNegativeStrand());

				}
			}
			if (next != null)
				processOverlap(next, searchSequence.length(), sequences.get(k).length());
		}
	}
	
	/**
	 * PRE: length1 > length2
	 * 
	 * @param overlap
	 * @param length1
	 * @param length2
	 */
	private void processOverlap(ReadOverlap overlap, int length1, int length2) {
		int l1 = overlap.getLast1() - overlap.getFirst1();
		int l2 = overlap.getLast2() - overlap.getFirst2();
		if (l1 < 0.2 * length1 && l2 < 0.2 * length2)
			return;
		int dl1 = overlap.getFirst1();
		int dl2 = overlap.getFirst2();
		int dr1 = length1 - overlap.getLast1();
		int dr2 = length2 - overlap.getLast2();
		int t = 10;
		int diffBorder = 100;

		if (overlap.isNegativeStrand()) {
			if (dl1 + t >= dr2 && dr1 + t >= dl2) {
				if (dl2 > diffBorder || dr2 > diffBorder) {
					// Not a true overlap
					return;
				}
				if (!embeddedOverlaps.containsKey(overlap.getIndexSequence2())) {
					ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(),
							overlap.getLast2(), overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(),
							true);
					embeddedOverlaps.put(overlap.getIndexSequence2(), ov2);
				}
				return;
			}
			if (dl1 <= dr2 + t && dr1 <= dl2 + t) {
				if (dl1 > diffBorder || dr1 > diffBorder) {
					// Not a true overlap
					return;
				}
				if (!embeddedOverlaps.containsKey(overlap.getIndexSequence1()))
					embeddedOverlaps.put(overlap.getIndexSequence1(), overlap);
				return;
			}
			if (dl1 < dr2 && dr1 > dl2) {
				// 2 left of 1 negative strand
				if (dl1 > diffBorder || dl2 > diffBorder) {
					// Not a true overlap
					return;
				}
				addOverlap(overlapsBackward, overlap);
				ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(),
						overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), true);
				addOverlap(overlapsForward, ov2);
				return;
			}
			if (dl1 > dr2 && dr1 < dl2) {
				// 1 left of 2
				if (dr2 > diffBorder || dr1 > diffBorder) {
					// Not a true overlap
					return;
				}
				addOverlap(overlapsForward, overlap);
				ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(),
						overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), true);
				addOverlap(overlapsBackward, ov2);
				return;
			}
		} else {
			if (dl1 + t >= dl2 && dr1 + t >= dr2) {
				if (dl2 > diffBorder || dr2 > diffBorder) {
					// Not a true overlap
					return;
				}
				if (!embeddedOverlaps.containsKey(overlap.getIndexSequence2())) {
					ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(),
							overlap.getLast2(), overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(),
							false);
					embeddedOverlaps.put(overlap.getIndexSequence2(), ov2);
				}
				return;
			}
			if (dl1 <= dl2 + t && dr1 <= dr2 + t) {
				if (dl1 > diffBorder || dr1 > diffBorder) {
					// Not a true overlap
					return;
				}
				if (!embeddedOverlaps.containsKey(overlap.getIndexSequence1()))
					embeddedOverlaps.put(overlap.getIndexSequence1(), overlap);
				return;
			}
			if (dl1 < dl2 && dr1 > dr2) {
				// 2 left of 1
				if (dl1 > diffBorder || dr2 > diffBorder) {
					// Not a true overlap
					return;
				}
				addOverlap(overlapsBackward, overlap);
				ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(),
						overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), false);
				addOverlap(overlapsForward, ov2);
				return;
			}
			if (dl1 > dl2 && dr1 < dr2) {
				// 1 left of 2
				if (dl2 > diffBorder || dr1 > diffBorder) {
					// Not a true overlap
					return;
				}
				addOverlap(overlapsForward, overlap);
				ReadOverlap ov2 = new ReadOverlap(overlap.getIndexSequence2(), overlap.getFirst2(), overlap.getLast2(),
						overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), false);
				addOverlap(overlapsBackward, ov2);
				return;
			}
		}
	}
	
	private void addOverlap(Map<Integer, List<ReadOverlap>> overlaps, ReadOverlap overlap) {
		List<ReadOverlap> overlapsSeq1 = overlaps.get(overlap.getIndexSequence1());
		if (overlapsSeq1 == null) {
			overlapsSeq1 = new ArrayList<>();
			overlaps.put(overlap.getIndexSequence1(), overlapsSeq1);
		}
		overlapsSeq1.add(overlap);
	}

	private AssemblyGraph buidAssemblyGraph() {
//		int[] map = new int[sequences.size()];
//		Arrays.fill(map, 1);
//		Queue<Integer> queue = new PriorityQueue<Integer>(
//				(Integer a, Integer b) -> a - b);
//		for (int i : embeddedOverlaps.keySet()) {
//			queue.add(i);
//			map[i] = 0;
//		}
//
//		for (int i = 1; i < map.length; i++)
//			map[i] += map[i - 1];
//
//		List<CharSequence> list = new LinkedList<CharSequence>();
//
//		int i = 0;
//		Integer t = queue.poll();
//		while (!queue.isEmpty()) {
//			while (i < t) {
//				list.add(sequences.get(i));
//				i++;
//			}
//			t = queue.poll();
//			i++;
//		}
//		while (i < sequences.size()) {
//			list.add(sequences.get(i));
//			i++;
//		}
//
//		System.out.println("		Sequences: " + sequences.size());
//		System.out.println("		Emmbeded Sequences: " + embeddedOverlaps.size());
//
//		AssemblyGraph assemblyGraph = new AssemblyGraph(list);
//		for (Entry<Integer, Embedded> entry : embeddedOverlaps.entrySet()) {
//			Embedded emb = entry.getValue();
//			int embeddedId = entry.getKey();
//			int parentId = emb.getIdSequence();
//			int parentPos = emb.getPos();
//			boolean isReverse = emb.isReversed();
//			AssembyEmbedded embedded = new AssembyEmbedded(
//					sequences.get(embeddedId), parentPos, isReverse);
//			assemblyGraph.addEmbedded(map[parentId], embedded);
//		}
//
//		int count = 0;
//		for (Overlap overlap : overlaps) {
//			if (!embeddedOverlaps.containsKey(overlap.getIdFrom())
//					&& !embeddedOverlaps.containsKey(overlap.getIdTo())) {
//				boolean startFrom = !overlap.isFromReversed();
//				boolean startTo = startFrom ^ !overlap.isFromReversed();
//				AssemblyVertex from = assemblyGraph.getVertex(
//						map[overlap.getIdFrom()], startFrom);
//				AssemblyVertex to = assemblyGraph.getVertex(
//						map[overlap.getIdTo()], startTo);
//				assemblyGraph.addEdge(to, from, overlap.getLength());
//				count++;
//			}
//		}
//
//		System.out.println("		Overlaps: " + count+"/"+overlaps.size());
//
		return assemblyGraph;
	}
}

class ReadOverlap {
	private int indexSequence1;
	private int indexSequence2;
	private int first1;
	private int last1;
	private int first2;
	private int last2;
	private boolean negativeStrand;

	public ReadOverlap(int indexSequence1, int first1, int last1, int indexSequence2, int first2, int last2,
			boolean negativeStrand) {
		super();
		this.indexSequence1 = indexSequence1;
		this.indexSequence2 = indexSequence2;
		this.first1 = first1;
		this.last1 = last1;
		this.first2 = first2;
		this.last2 = last2;
		this.negativeStrand = negativeStrand;
	}

	public boolean addKmerAlignment(int idxSeq2, ReadAlignment aln) {
		if (aln == null)
			return false;
		if (idxSeq2 != indexSequence2)
			return false;
		if (aln.isNegativeStrand() != negativeStrand)
			return false;
		int kmerFirst = aln.getMateFirst();
		int kmerLast = kmerFirst + GraphBuilder_OverlapFinder_Quee.SEARCH_KMER_LENGTH;
		int expectedAlnFirst = last2 + (kmerFirst - last1);
		int expectedAlnLast = last2 + (kmerLast - last1);
		if (negativeStrand) {
			expectedAlnFirst = first2 - (kmerLast - last1);
			expectedAlnLast = first2 - (kmerFirst - last1);
		}
		if (Math.abs(aln.getFirst() - expectedAlnFirst) > 5)
			return false;
		if (Math.abs(aln.getLast() - expectedAlnLast) > 5)
			return false;
		last1 = kmerLast;
		if (negativeStrand) {
			first2 = aln.getFirst();
		} else {
			last2 = aln.getLast();
		}
		return true;
	}

	public int getIndexSequence1() {
		return indexSequence1;
	}

	public int getIndexSequence2() {
		return indexSequence2;
	}

	public int getFirst1() {
		return first1;
	}

	public int getLast1() {
		return last1;
	}

	public int getFirst2() {
		return first2;
	}

	public int getLast2() {
		return last2;
	}

	public int length() {
		return last1 - first2;
	}

	public boolean isNegativeStrand() {
		return negativeStrand;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		// return "" + indexSequence2 + "(" + first1 + "," + last1 + ")";
		return "ReadOverlap [indexSequence1=" + indexSequence1 + ", indexSequence2=" + indexSequence2 + ", first1="
				+ first1 + ", last1=" + last1 + ", first2=" + first2 + ", last2=" + last2 + ", negativeStrand="
				+ negativeStrand + "]";
	}
}