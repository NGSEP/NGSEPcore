package ngsep.assembly;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.FMIndex;
import static ngsep.assembly.TimeUtilities.timeGroup;
import static ngsep.assembly.TimeUtilities.timeIt;
import static ngsep.assembly.TimeUtilities.progress;

public class GraphBuilderOverlapFinderQueue2 implements GraphBuilderOverlapFinder{
	private final static double Rate_of_changes = 0.02;
	private final static double Rate_of_cuts = 0.01;
	private final static double Rate_of_cover = 7;
	private final static double MinKmerCoverRate = -1;

	private GraphBuilderKmerIterator kmerIterator;
	private List<CharSequence> sequences;
	private SimplifiedAssemblyGraph assemblyGraph;
	private FMIndex index;
	private Map<Integer, List<int[]>> hits;

	@Override
	public void calculate(List<CharSequence> seq, FMIndex index) {
		timeGroup("--------- Build Graph --------", () -> {
			this.kmerIterator = new GraphBuilderKmerIterator(Rate_of_changes, Rate_of_cuts, Rate_of_cover);
			this.hits = new HashMap<Integer, List<int[]>>(seq.size());
			this.index = index;
			this.sequences = seq;
			assemblyGraph = new SimplifiedAssemblyGraph(seq);

			printRates();
			timeIt(" --- Find overlaps ", () -> findOverlapsAndEmbedded());
			assemblyGraph.removeAllEmbeddedsIntoGraph();
		});
		assemblyGraph.printInfo();
	}

	public void printRates() {
		System.out.println(" | SEARCH_KMER_LENGTH: " + kmerIterator.SEARCH_KMER_LENGTH + " |");
		System.out.println(" | SEARCH_KMER_DISTANCE: " + kmerIterator.SEARCH_KMER_DISTANCE + " |");
		System.out.println(" | MAX_KMER_DES: " + kmerIterator.MAX_KMER_DES + " |");
	}

	@Override
	public SimplifiedAssemblyGraph getGrap() {
		return assemblyGraph;
	}
	

	public void findOverlapsAndEmbedded() {
		for (int seqId = 0, excp = 0; seqId < sequences.size(); seqId++) {
			progress(" --- Find overlaps ", seqId + assemblyGraph.amuontOfEmbeddedSequences() - excp, sequences.size());
			if (assemblyGraph.isEmbedded(seqId)) {
				excp++;
				continue;
			}

			calculateHits(seqId, kmerIterator.positiveStrand(sequences.get(seqId)));
			for (Entry<Integer, List<int[]>> entry : hits.entrySet())
				Aling(seqId, entry.getKey(), false, entry.getValue());

			calculateHits(seqId, kmerIterator.negativeStrand(sequences.get(seqId)));
			for (Entry<Integer, List<int[]>> entry : hits.entrySet())
				Aling(seqId, entry.getKey(), true, entry.getValue());
		}
	}

	private void calculateHits(int id_Ref, Iterable<Entry<Integer, String>> kmerIters) {
		hits.clear();
		for (Entry<Integer, String> entry : kmerIters) {
			int pos_Ref = entry.getKey();
			for (ReadAlignment readAlignment : index.search(entry.getValue())) {
				int id_Lec = Integer.parseInt(readAlignment.getSequenceName());
				int pos_Lec = readAlignment.getFirst();

				if (id_Ref < id_Lec)
					hits.computeIfAbsent(id_Lec, x -> new LinkedList<int[]>()).add(new int[] { pos_Ref, pos_Lec });
			}
		}
	}

	private void Aling(int id_Ref, int id_Lec, boolean isReverse, List<int[]> hits) {
		if (hits.size() < 2)
			return;

		TreeMap<Integer, int[]> tree = new TreeMap<Integer, int[]>();
		Set<Integer> remove = new HashSet<>();
		Stack<Integer> keys = new Stack<>();
		Stack<int[]> values = new Stack<>();
		int[] aux;

		for (Entry<Integer, List<Integer>> entry : groupByPosRef(hits).entrySet()) {
			int pRef = entry.getKey();

			for (int posibleKey : entry.getValue()) {
				int key = closestValidKey(tree, posibleKey, kmerIterator.MAX_KMER_DES);
				if (key == -1)
					aux = new int[] { 1, pRef, posibleKey + pRef };
				else {
					remove.add(key);
					aux = tree.get(key);
					aux[0]++;
				}
				keys.add(posibleKey);
				values.add(aux);
			}
			removeAll(tree, remove);
			addAll(tree, keys, values);
		}

		for (int[] aln : tree.values())
			if (aln[0] > 1) {
				keys.add(aln[2] - aln[1]);
				values.add(aln);
			}
		tree.clear();
		addAll(tree, keys, values);
		dectect(id_Ref, id_Lec, isReverse, tree);
	}

	private void dectect(int id_Ref, int id_Lec, boolean isReverse, TreeMap<Integer, int[]> tree) {
		int lenghtRef = sequences.get(id_Ref).length();
		int lenghtLec = sequences.get(id_Lec).length();
		int embbedLimit = lenghtLec - lenghtRef;

		// lec is embedded in ref
		for (int[] aln : tree.subMap(embbedLimit, true, 0, true).values()) {
			int pos_Lec = aln[2] - aln[1];
			double rate = kmerIterator.SEARCH_KMER_LENGTH / (double) numberOfKmers(lenghtLec);
			if (rate < MinKmerCoverRate)
				continue;

			assemblyGraph.addEmbedded(id_Ref, id_Lec, isReverse ? lenghtRef + pos_Lec - lenghtLec : -pos_Lec, isReverse,
					rate);
			return;
		}

		// lec -> ref || lec -> ref'
		for (int[] aln : tree.subMap(0, true, Integer.MAX_VALUE, true).values()) {
			int pos_Lec = aln[2] - aln[1];
			double rate = kmerIterator.SEARCH_KMER_LENGTH / (double) numberOfKmers(lenghtLec - pos_Lec);
			if (rate < MinKmerCoverRate)
				continue;

			assemblyGraph.addEdge((id_Lec << 1) + 1, (id_Ref << 1) + (isReverse ? 1 : 0), lenghtLec - pos_Lec, rate);
			break;
		}

		// ref -> lec || ref' -> lec
		for (int[] aln : tree.descendingMap().subMap(embbedLimit, true, Integer.MIN_VALUE, true).values()) {
			int pos_Lec = aln[2] - aln[1];
			double rate = kmerIterator.SEARCH_KMER_LENGTH / (double) numberOfKmers(lenghtRef + pos_Lec);
			if (rate < MinKmerCoverRate)
				continue;

			assemblyGraph.addEdge((id_Ref << 1) + (isReverse ? 0 : 1), id_Lec << 1, lenghtRef + pos_Lec, rate);
			break;
		}
	}

	private static Map<Integer, List<Integer>> groupByPosRef(List<int[]> hits) {
		Map<Integer, List<Integer>> group = new HashMap<Integer, List<Integer>>();
		int prev = -1;
		List<Integer> list = null;
		for (int[] hit : hits) {
			if (prev != hit[0]) {
				prev = hit[0];
				list = new LinkedList<>();
				group.put(prev, list);
			}
			list.add(hit[1] - prev);
		}
		return group;
	}

	private static <K> void removeAll(TreeMap<K, ? extends Object> tree, Set<K> set) {
		for (K i : set)
			tree.remove(i);
		set.clear();
	}

	private static <T, K> void addAll(TreeMap<K, T> tree, Stack<K> keys, Stack<T> values) {
		while (!keys.isEmpty())
			tree.put(keys.pop(), values.pop());
	}

	private int numberOfKmers(int size) {
		return size / (kmerIterator.SEARCH_KMER_LENGTH - kmerIterator.SEARCH_KMER_DISTANCE);
	}

	private static int closestValidKey(TreeMap<Integer, ? extends Object> treeMap, int key, int diff) {
		int ans = -1;
		int min = diff;
		Integer celingKey = treeMap.ceilingKey(key), floorKey = treeMap.floorKey(key);

		if (celingKey != null && celingKey - key < min) {
			min = celingKey - key;
			ans = celingKey;
		}
		if (floorKey != null && key - floorKey < min) {
			min = key - floorKey;
			ans = floorKey;
		}
		return ans;
	}
}
