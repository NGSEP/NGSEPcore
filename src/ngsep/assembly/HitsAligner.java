package ngsep.assembly;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.Map.Entry;

@FunctionalInterface
public interface HitsAligner {
	/**
	 * Add to a graph if the set of hits is an align
	 * 
	 * @param id_Ref    id reference sequence
	 * @param id_Lec    id sample sequence
	 * @param isReverse if reference sequence is reversed
	 * @param hits      a set of {posRef,posLect}
	 */
	void Aling(int id_Ref, int id_Lec, boolean isReverse, List<int[]> hits);
}

class TreesHitAligner implements HitsAligner {
	private SimplifiedAssemblyGraph sag;
	private AssemblyConfiguration config;
	private List<CharSequence> sequences;

	/*
	 * tree with key relative position of sample sequence into reference sequence
	 * the value is {countOfHits,posRef,posLect} (the first positions)
	 */
	private TreeMap<Integer, int[]> tree = new TreeMap<Integer, int[]>();
	/**
	 * set wit keys to remove
	 */
	private Set<Integer> remove = new HashSet<>();
	/**
	 * stacks with pending values in for the tree
	 */
	private Stack<Integer> keys = new Stack<>();
	private Stack<int[]> values = new Stack<>();

	public TreesHitAligner(SimplifiedAssemblyGraph sag, AssemblyConfiguration config, List<CharSequence> sequences) {
		this.sag = sag;
		this.config = config;
		this.sequences = sequences;
	}

	@Override
	public void Aling(int id_Ref, int id_Lec, boolean isReverse, List<int[]> hits) {
		clear();
		int[] aux;
		for (Entry<Integer, List<Integer>> entry : groupByPosRef(hits).entrySet()) {
			int pRef = entry.getKey();

			for (int posibleKey : entry.getValue()) {
				int key = closestValidKey(tree, posibleKey, config.overlap().getMaxKmerDiff());
				if (key == -1)
					aux = new int[] { 1, pRef, posibleKey + pRef, 0, 0 };
				else {
					remove.add(key);
					aux = tree.get(key);
					aux[0]++;
					aux[3] = pRef;
					aux[4] = posibleKey + pRef;
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

	private void clear() {
		tree.clear();
		remove.clear();
		keys.clear();
		values.clear();
	}

	private void dectect(int id_Ref, int id_Lec, boolean isReverse, TreeMap<Integer, int[]> tree) {
		int borderRate = 100 * config.overlap().getKmerLength();

		int lenghtRef = sequences.get(id_Ref).length();
		int lenghtLec = sequences.get(id_Lec).length();
		int embbedLimit = lenghtLec - lenghtRef;

		// lec is embedded in ref
		for (int[] aln : tree.subMap(embbedLimit, true, 0, true).values()) {
			int pos_Lec = aln[2] - aln[1];
			double rate = aln[0] / (double) numberOfKmers(lenghtLec);

			if (rate < config.overlap().getMinKmerCoverRate())
				continue;

			if (aln[2] > borderRate || (aln[4] + config.overlap().getKmerLength()) < lenghtLec - borderRate)
				continue;

			sag.addEmbedded(id_Ref, id_Lec, isReverse ? lenghtRef + pos_Lec - lenghtLec : -pos_Lec, isReverse, rate);
			return;
		}

		// lec -> ref || lec -> ref'
		for (int[] aln : tree.subMap(0, true, Integer.MAX_VALUE, true).values()) {
			int pos_Lec = aln[2] - aln[1];
			int len = lenghtLec - pos_Lec;
			double rate = aln[0] / (double) numberOfKmers(len);

			if (rate < config.overlap().getMinKmerCoverRate())
				continue;

			if (aln[1] > borderRate || (aln[4] + config.overlap().getKmerLength() - pos_Lec) < len - borderRate)
				continue;

			sag.addEdge((id_Lec << 1) + 1, (id_Ref << 1) + (isReverse ? 1 : 0), lenghtLec - pos_Lec, rate);
			break;
		}

		// ref -> lec || ref' -> lec
		for (int[] aln : tree.descendingMap().subMap(embbedLimit, true, Integer.MIN_VALUE, true).values()) {
			int pos_Lec = aln[2] - aln[1];
			int len = lenghtRef + pos_Lec;
			double rate = aln[0] / (double) numberOfKmers(len);

			if (rate < config.overlap().getMinKmerCoverRate())
				continue;

			if (aln[2] > borderRate || (aln[3] + config.overlap().getKmerLength() + pos_Lec) < len - borderRate)
				continue;

			sag.addEdge((id_Ref << 1) + (isReverse ? 0 : 1), id_Lec << 1, lenghtRef + pos_Lec, rate);
			break;
		}
	}

	private int numberOfKmers(int size) {
		return size / (config.overlap().getKmerLength() + config.overlap().getKmerDistance());
	}

	private static TreeMap<Integer, List<Integer>> groupByPosRef(List<int[]> hits) {
		TreeMap<Integer, List<Integer>> group = new TreeMap<Integer, List<Integer>>();
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
