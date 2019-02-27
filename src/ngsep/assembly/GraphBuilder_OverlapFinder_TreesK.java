package ngsep.assembly;

import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.FMIndex;

public class GraphBuilder_OverlapFinder_TreesK implements
		GraphBuilder_OverlapFinder {
	private static final int NUMBER_INVALID_KMERS = 5;
	/** represent the maximum bleed + 1 between two kmers of the same hit */
	private static final int MAX_KMER_BLEED = 6;
	/** Possible end of Fastq file */

	private List<CharSequence> sequences;
	private AssemblyGraph assemblyGraph;

	private TreeMap<Integer, int[]>[][] hits;
	private Map<Integer, ReadOverlap> embeddedOverlaps = new Hashtable<>();

	@Override
	@SuppressWarnings("unchecked")
	public void calculate(List<CharSequence> seq, FMIndex index) {
		this.sequences = seq;
		hits = new TreeMap[NUMBER_INVALID_KMERS + 1][sequences.size()];
		for (int y = 0; y <= NUMBER_INVALID_KMERS; y++)
			for (int x = 0; x < sequences.size(); x++)
				hits[y][x] = new TreeMap<>();
		emmbeddeDetection(seq, index);
		assemblyGraph = buidAssemblyGraph();
	}

	@Override
	public AssemblyGraph getGrap() {
		return assemblyGraph;
	}

	private void emmbeddeDetection(List<CharSequence> sequencesToIndex,
			FMIndex index) {
		int n;
		int k;
		Iterator<String> iter;
		boolean[] mark = new boolean[sequencesToIndex.size()];
		Hashtable<Integer, List<int[]>> hitsPoint = new Hashtable<>();
		for (int idSequence = 0; idSequence < sequencesToIndex.size(); idSequence++) {
//			System.out.println(idSequence);
			hitsPoint.clear();
			KmerEmmbeddedIterator ki = new KmerEmmbeddedIterator(
					sequencesToIndex.get(idSequence));
			n = (int) Math.round(ki.getNumber()
					* KmerEmmbeddedIterator.PORCENTAJE_OF_TOLERANCE);
			if (n <= 0)
				n = 1;

			k = 0;
			iter = ki.firts().iterator();
			for (int i = 0; i < n; i++) {
				for (ReadAlignment aln : index.search(iter.next(),
						idSequence + 1, sequencesToIndex.size()))
					alignAdd(aln, k);
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}

			for (int i = idSequence + 1; i < mark.length; i++) {
				boolean a = false;
				for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
					if (!hits[j][i].isEmpty()) {
						a = true;
						break;
					}
				mark[i] = a;
			}
			while (iter.hasNext()) {
				for (ReadAlignment aln : index.search(iter.next(),
						idSequence + 1, sequencesToIndex.size())) {
					if (mark[Integer.parseInt(aln.getSequenceName())])
						align(aln, k);
				}
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}
			for (int i = idSequence + 1; i < mark.length; i++) {
				if (mark[i]) {
					boolean a = false;
					for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
						if (!hits[j][i].isEmpty()) {
							for (int[] array : hits[j][i].values()) {
								if (array[0]
										+ sequencesToIndex.get(idSequence)
												.length() <= sequencesToIndex
										.get(i).length()) {
									a = true;
									hitsPoint.computeIfAbsent(i,
											(o) -> new LinkedList<>()).add(
											array);
								}
							}
							hits[j][i].clear();
						}
					mark[i] = a;
				}
			}
			k = 0;
			iter = ki.firts().iterator();
			for (int i = 0; i < n; i++) {
				for (ReadAlignment aln : index.search(iter.next(),
						idSequence + 1, sequencesToIndex.size()))
					alignAdd(aln, k);
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}
			for (int i = idSequence + 1; i < mark.length; i++) {
				boolean a = false;
				for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
					if (!hits[j][i].isEmpty()) {
						a = true;
						break;
					}
				mark[i] = a;
			}
			while (iter.hasNext()) {
				for (ReadAlignment aln : index.search(iter.next(),
						idSequence + 1, sequencesToIndex.size())) {
					if (mark[Integer.parseInt(aln.getSequenceName())])
						align(aln, k);
				}
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}
			for (int i = idSequence + 1; i < mark.length; i++) {
				if (mark[i]) {
					boolean a = false;
					for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
						if (!hits[j][i].isEmpty()) {
							for (int[] array : hits[j][i].values()) {
								if (array[0]
										+ sequencesToIndex.get(idSequence)
												.length() <= sequencesToIndex
										.get(i).length()) {
									a = true;
									hitsPoint.computeIfAbsent(i,
											(o) -> new LinkedList<>()).add(
											array);
								}
							}
							hits[j][i].clear();
						}
					mark[i] = a;
				}
			}

			int falt = ki.getNumber()
					* KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
			k = sequencesToIndex.get(idSequence).length() - falt;
			iter = ki.lasts().iterator();
			for (int i = 0; i < n; i++) {
				for (ReadAlignment aln : index.search(iter.next(),
						idSequence + 1, sequencesToIndex.size()))
					if (mark[Integer.parseInt(aln.getSequenceName())]) {
						if (aln.getFirst() + falt <= sequencesToIndex.get(
								Integer.parseInt(aln.getSequenceName()))
								.length())
							alignAdd(aln, k);
					}
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				falt -= KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}
			for (int i = idSequence + 1; i < mark.length; i++) {
				if (mark[i]) {
					boolean a = false;
					for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
						if (!hits[j][i].isEmpty()) {
							a = true;
							break;
						}
					mark[i] = a;
				}
			}
			while (iter.hasNext()) {
				for (ReadAlignment aln : index.search(iter.next(),
						idSequence + 1, sequencesToIndex.size())) {
					if (mark[Integer.parseInt(aln.getSequenceName())])
						align(aln, k);
				}
				k += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH;
				rotate();
			}

			for (int i = idSequence + 1; i < mark.length; i++) {
				if (mark[i]) {
					for (int j = 1; j <= NUMBER_INVALID_KMERS; j++)
						if (!hits[j][i].isEmpty()) {
							for (int[] array : hits[j][i].values()) {
								for (int[] jj : hitsPoint.get(i))
									if (Math.abs(array[0]
											+ (ki.getNumber() * KmerEmmbeddedIterator.SEARCH_KMER_LENGTH)
											- jj[0]
											- sequencesToIndex.get(idSequence)
													.length()) < sequencesToIndex
											.get(idSequence).length()
											/ (double) 10) {
										embeddedOverlaps
												.put(idSequence,
														new ReadOverlap(
																idSequence,
																0,
																sequencesToIndex
																		.get(idSequence)
																		.length() - 1,
																i,
																jj[0],
																array[0]
																		+ (ki.getNumber() * KmerEmmbeddedIterator.SEARCH_KMER_LENGTH),
																false));
										break;
									}
							}
							hits[j][i].clear();
						}
				}
			}
		}
	}

	/**
	 * try to align the read whit each possible aligns in a
	 * 
	 * @param aln
	 *            the read (one alignment for the kmer)
	 * @param posinRefrenece
	 *            index of the kmer in the original sequence.
	 */
	private void alignAdd(ReadAlignment aln, int posinRefrenece) {
		int idSequenceAligned = Integer.parseInt(aln.getSequenceName());
		int pos = aln.getFirst();

		int minposkmer = 0;
		int minnode = 0;
		int min = MAX_KMER_BLEED;
		for (int dif = KmerEmmbeddedIterator.SEARCH_KMER_LENGTH, i = 1; i <= NUMBER_INVALID_KMERS; i++, dif += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH) {
			TreeMap<Integer, int[]> treeMap = hits[i][idSequenceAligned];
			if (treeMap.isEmpty())
				continue;

			int t = pos - dif;
			Integer p1 = treeMap.ceilingKey(t);
			Integer p2 = treeMap.floorKey(t);
			if (p1 != null && Math.abs(p1 - t) < min) {
				min = Math.abs(p1 - t);
				minposkmer = i;
				minnode = p1;
			}
			if (p2 != null && Math.abs(p2 - t) < min) {
				min = Math.abs(p2 - t);
				minposkmer = i;
				minnode = p2;
			}
		}

		if (min == MAX_KMER_BLEED) {
			hits[0][idSequenceAligned].put(pos,
					new int[] { pos, posinRefrenece });
		} else {
			int[] aux = hits[minposkmer][idSequenceAligned].remove(minnode);
			hits[0][idSequenceAligned].put(pos, aux);
		}
	}

	/**
	 * try to align the read whit each possible aligns in a
	 * 
	 * @param aln
	 *            the read (one alignment for the kmer)
	 * @param posinRefrenece
	 *            index of the kmer in the original sequence.
	 */
	private void align(ReadAlignment aln, int posinRefrenece) {
		int idSequenceAligned = Integer.parseInt(aln.getSequenceName());
		int pos = aln.getFirst();

		int minposkmer = 0;
		int minnode = 0;
		int min = MAX_KMER_BLEED;
		for (int dif = KmerEmmbeddedIterator.SEARCH_KMER_LENGTH, i = 1; i <= NUMBER_INVALID_KMERS; i++, dif += KmerEmmbeddedIterator.SEARCH_KMER_LENGTH) {
			TreeMap<Integer, int[]> treeMap = hits[i][idSequenceAligned];
			if (treeMap.isEmpty())
				continue;

			int t = pos - dif;
			Integer p1 = treeMap.ceilingKey(t);
			Integer p2 = treeMap.floorKey(t);
			if (p1 != null && Math.abs(p1 - t) < min) {
				min = Math.abs(p1 - t);
				minposkmer = i;
				minnode = p1;
			}
			if (p2 != null && Math.abs(p2 - t) < min) {
				min = Math.abs(p2 - t);
				minposkmer = i;
				minnode = p2;
			}
		}

		if (min != MAX_KMER_BLEED) {
			int[] aux = hits[minposkmer][idSequenceAligned].remove(minnode);
			hits[0][idSequenceAligned].put(pos, aux);
		}
	}

	/**
	 * rotate the columns of hits
	 */
	private void rotate() {
		TreeMap<Integer, int[]>[] temp = hits[NUMBER_INVALID_KMERS];
		for (TreeMap<Integer, int[]> t : temp)
			t.clear();
		System.arraycopy(hits, 0, hits, 1, NUMBER_INVALID_KMERS);
		hits[0] = temp;
	}

	private AssemblyGraph buidAssemblyGraph() {
		// int[] map = new int[sequences.size()];
		// Arrays.fill(map, 1);
		// Queue<Integer> queue = new PriorityQueue<Integer>(
		// (Integer a, Integer b) -> a - b);
		// for (int i : embeddedOverlaps.keySet()) {
		// queue.add(i);
		// map[i] = 0;
		// }
		//
		// for (int i = 1; i < map.length; i++)
		// map[i] += map[i - 1];
		//
		// List<CharSequence> list = new LinkedList<CharSequence>();
		//
		// int i = 0;
		// Integer t = queue.poll();
		// while (!queue.isEmpty()) {
		// while (i < t) {
		// list.add(sequences.get(i));
		// i++;
		// }
		// t = queue.poll();
		// i++;
		// }
		// while (i < sequences.size()) {
		// list.add(sequences.get(i));
		// i++;
		// }
		//
		// System.out.println("		Sequences: " + sequences.size());
		// System.out.println("		Emmbeded Sequences: " +
		// embeddedOverlaps.size());
		//
		// AssemblyGraph assemblyGraph = new AssemblyGraph(list);
		// for (Entry<Integer, Embedded> entry : embeddedOverlaps.entrySet()) {
		// Embedded emb = entry.getValue();
		// int embeddedId = entry.getKey();
		// int parentId = emb.getIdSequence();
		// int parentPos = emb.getPos();
		// boolean isReverse = emb.isReversed();
		// AssembyEmbedded embedded = new AssembyEmbedded(
		// sequences.get(embeddedId), parentPos, isReverse);
		// assemblyGraph.addEmbedded(map[parentId], embedded);
		// }
		//
		// int count = 0;
		// for (Overlap overlap : overlaps) {
		// if (!embeddedOverlaps.containsKey(overlap.getIdFrom())
		// && !embeddedOverlaps.containsKey(overlap.getIdTo())) {
		// boolean startFrom = !overlap.isFromReversed();
		// boolean startTo = startFrom ^ !overlap.isFromReversed();
		// AssemblyVertex from = assemblyGraph.getVertex(
		// map[overlap.getIdFrom()], startFrom);
		// AssemblyVertex to = assemblyGraph.getVertex(
		// map[overlap.getIdTo()], startTo);
		// assemblyGraph.addEdge(to, from, overlap.getLength());
		// count++;
		// }
		// }
		//
		// System.out.println("		Overlaps: " + count+"/"+overlaps.size());

		return assemblyGraph;
	}
}

class KmerEmmbeddedIterator {
	public static final int NUMBER_OF_KMERS = 10;
	public static final int SEARCH_KMER_LENGTH = 10;
	public static final double PORCENTAJE_OF_KMERS = 0.1;
	public static final double PORCENTAJE_OF_TOLERANCE = 0.3;
	public static final double FACTOR = PORCENTAJE_OF_KMERS
			/ (double) SEARCH_KMER_LENGTH;
	private final String s;
	private final int number;

	public KmerEmmbeddedIterator(CharSequence s) {
		this.s = s.toString();
		number = (int) Math.round(s.length() * FACTOR);
	}

	public Iterable<String> firts() {
		return new Iterable<String>() {
			@Override
			public Iterator<String> iterator() {
				return new Iterator<String>() {
					private int i = 0;
					private int c = 0;

					@Override
					public boolean hasNext() {
						return c < number;
					}

					@Override
					public String next() {
						String ans = s.substring(i, i + SEARCH_KMER_LENGTH - 1);
						i += SEARCH_KMER_LENGTH;
						c++;
						return ans;
					}
				};
			}
		};
	}

	public Iterable<String> lasts() {
		return new Iterable<String>() {
			@Override
			public Iterator<String> iterator() {
				return new Iterator<String>() {
					private int i = s.length();
					private int c = 0;

					@Override
					public boolean hasNext() {
						return c < number;
					}

					@Override
					public String next() {
						String ans = s.substring(i - SEARCH_KMER_LENGTH - 1, i);
						i -= SEARCH_KMER_LENGTH;
						c++;
						return ans;
					}
				};
			}
		};
	}

	public int getNumber() {
		return number;
	}
}
