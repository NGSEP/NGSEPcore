package ngsep.assembly;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.TreeMap;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.FMIndex;

public class GraphBuilderOverlapFinderTree implements GraphBuilderOverlapFinder {
	private final static double Rate_of_changes = 0.07;
	private final static double Rate_of_cuts = 0.03;
	private final static double Rate_of_cover = 2;
	public static final double permitedBorderRate = 0.15;

	private List<CharSequence> sequences;
	private AssemblyGraph assemblyGraph;
	private FMIndex index;

	private Map<Integer, Embedded> embeddedOverlaps = new Hashtable<>();
	private Map<Integer, TreeMap<Integer, Alignment>> alignments;
	private Map<Integer, TreeMap<Integer, Alignment>> alignmentsAux;
	private Map<Integer, List<Integer>> removeAux;
	private Deque<Overlap> overlaps = new ArrayDeque<Overlap>();
	private GraphBuilderKmerIterator kmerIterator;

	@Override
	public void calculate(List<CharSequence> seq, FMIndex index) {
		this.kmerIterator = new GraphBuilderKmerIterator(Rate_of_changes, Rate_of_cuts, Rate_of_cover);
		this.index = index;
		this.sequences = seq;
		System.out.println("		Finding overlaps");
		long ini = System.currentTimeMillis();
		findOverlaps();
		System.out.println("		found overlaps: " + (System.currentTimeMillis() - ini) / (double) 1000 + " s");

		System.out.println("		Building graph");
		ini = System.currentTimeMillis();
		assemblyGraph = buidAssemblyGraph();
		System.out.println("		built graph: " + (System.currentTimeMillis() - ini) / (double) 1000 + " s");

		try (PrintWriter pr = new PrintWriter(new File("grafo.txt"))) {
			for (AssemblyEdge edge : assemblyGraph.getEdges()) {
				pr.println(edge.getVertex1().getIndex() + "-" + edge.getVertex2().getIndex());
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	@Override
	public AssemblyGraph getGrap() {
		return assemblyGraph;
	}

	public void findOverlaps() {
		alignments = new HashMap<>(sequences.size());
		alignmentsAux = new HashMap<>(sequences.size());
		removeAux = new HashMap<>(sequences.size());
		int Excp = 0;
		int diff = Math.min(1, (int) (sequences.size() * 0.0025));
		int prev = 0;
		for (int seqId = 0; seqId < sequences.size(); seqId++) {
			int p = seqId + embeddedOverlaps.size() - Excp;
			if (p >= prev) {
				printProgressBar(System.out, "		", p, sequences.size());
				prev += diff;
			}

			if (embeddedOverlaps.containsKey(seqId)) {
				Excp++;
				continue;
			}

			alignments.clear();
			findAlignments(seqId, kmerIterator.positiveStrand(sequences.get(seqId)));
			detectOverlap(seqId, false);

			alignments.clear();
			findAlignments(seqId, kmerIterator.negativeStrand(sequences.get(seqId)));
			detectOverlap(seqId, true);
		}
		printProgressBar(System.out, "		", sequences.size(), sequences.size());
		System.out.println();
	}

	private void findAlignments(int id_Ref, Iterable<Entry<Integer, String>> kmerIters) {
		for (Entry<Integer, String> entry : kmerIters) {
			int pos_Ref = entry.getKey();
			alignmentsAux.clear();
			removeAux.clear();
			for (ReadAlignment readAlignment : index.search(entry.getValue())) {
				int id_Lec = Integer.parseInt(readAlignment.getSequenceName());
				int pos_Lec = readAlignment.getFirst();

				if (id_Ref < id_Lec)
					addAlingments(id_Ref, pos_Ref, id_Lec, pos_Lec);
			}
			for (Entry<Integer, TreeMap<Integer, Alignment>> entry2 : alignments.entrySet()) {
				int id_Lec = entry2.getKey();
				TreeMap<Integer, Alignment> treeMap = entry2.getValue();

				List<Integer> removes = removeAux.get(id_Lec);
				if (removes != null)
					for (int r : removes)
						if (treeMap.containsKey(r))
							treeMap.remove(r);

				TreeMap<Integer, Alignment> treeAux = alignmentsAux.get(id_Lec);
				if (treeAux != null)
					treeMap.putAll(treeAux);
			}
		}
	}

	private void addAlingments(int id_Ref, int pos_Ref, int id_Lec, int pos_Lec) {
		if (embeddedOverlaps.containsKey(id_Lec))
			return;
		TreeMap<Integer, Alignment> treeMap = alignments.computeIfAbsent(id_Lec,
				k -> new TreeMap<Integer, Alignment>());
		TreeMap<Integer, Alignment> treeAux = alignmentsAux.computeIfAbsent(id_Lec,
				k -> new TreeMap<Integer, Alignment>());
		List<Integer> removes = removeAux.computeIfAbsent(id_Lec, k -> new LinkedList<Integer>());

		int key = pos_Lec - pos_Ref;
		Alignment aln = aling(treeMap, removes, key);

		if (aln != null) {
			aln.setLengthLec(kmerIterator.SEARCH_KMER_LENGTH + (pos_Lec - aln.getPosLec()));
			aln.setLengthRef(kmerIterator.SEARCH_KMER_LENGTH + (pos_Ref - aln.getPosRef()));
			aln.addhit();
			treeAux.put(key, aln);
		} else {
			treeAux.put(key,
					new Alignment(pos_Ref, pos_Lec, kmerIterator.SEARCH_KMER_LENGTH, kmerIterator.SEARCH_KMER_LENGTH));
		}
	}

	private Alignment aling(TreeMap<Integer, Alignment> treeMap, List<Integer> removes, int key) {
		int removeKey = key;
		int min = kmerIterator.MAX_KMER_DES;
		Integer celingKey = treeMap.ceilingKey(key), floorKey = treeMap.floorKey(key);

		if (celingKey != null && celingKey - key < min) {
			min = celingKey - key;
			removeKey = celingKey;
		}
		if (floorKey != null && key - floorKey < min) {
			min = key - floorKey;
			removeKey = floorKey;
		}

		removes.add(removeKey);
		return (min != kmerIterator.MAX_KMER_DES) ? treeMap.get(removeKey) : null;
	}

	private void detectOverlap(int id_Ref, boolean isReverse) {
		double max = 0.12;
		int lenghtRef = sequences.get(id_Ref).length();
		int Diff = kmerIterator.MAX_KMER_DES;
		OtherSequence: for (Entry<Integer, TreeMap<Integer, Alignment>> entry : alignments.entrySet()) {
			int id_Lec = entry.getKey();
			int lenghtLec = sequences.get(id_Lec).length();

			int embbedLimit = lenghtLec - lenghtRef;
			TreeMap<Integer, Alignment> tree = entry.getValue();

			for (Alignment aln : tree.subMap(embbedLimit - Diff, true, 0 + Diff, true).values()) {

				double rate = kmerIterator.SEARCH_KMER_LENGTH * aln.getHits() / (double) (lenghtLec * Rate_of_cover);
				if (aln.getHits() < 2 || rate < max)
					continue;

				int pos_Lec = aln.getPosLec() - aln.getPosRef();
				if (pos_Lec > 0 || pos_Lec < embbedLimit)
					continue;
				if (aln.getPosLec() > lenghtLec * permitedBorderRate
						|| aln.getPosLec() + aln.getLengthLec() < lenghtLec * (1 - permitedBorderRate))
					continue;

				// lec is embedded in ref.
				if (!embeddedOverlaps.containsKey(id_Lec))
					embeddedOverlaps.put(id_Lec, new Embedded(id_Ref,
							isReverse ? lenghtRef + pos_Lec - lenghtLec : -pos_Lec, isReverse, rate));
				continue OtherSequence;
			}

			for (Alignment aln : tree.subMap(0 - Diff, true, Integer.MAX_VALUE, true).values()) {
				int pos_Lec = aln.getPosLec() - aln.getPosRef();
				double rate = kmerIterator.SEARCH_KMER_LENGTH * aln.getHits()
						/ (double) (Rate_of_cover * (lenghtLec - pos_Lec));
				if (aln.getHits() < 2 || rate < max)
					continue;

				if (pos_Lec < 0)
					continue;
				if (aln.getPosRef() > lenghtRef * permitedBorderRate
						|| aln.getPosLec() + aln.getLengthLec() < lenghtLec * (1 - permitedBorderRate))
					continue;

				// lec -> ref || lec -> ref'
				overlaps.add(new Overlap(id_Lec, false, id_Ref, isReverse, lenghtLec - pos_Lec, rate));
				break;
			}

			for (Alignment aln : tree.descendingMap().subMap(embbedLimit + Diff, true, Integer.MIN_VALUE, true)
					.values()) {
				int pos_Lec = aln.getPosLec() - aln.getPosRef();
				double rate = kmerIterator.SEARCH_KMER_LENGTH * aln.getHits()
						/ (double) (Rate_of_cover * (lenghtRef + pos_Lec));
				if (aln.getHits() < 2 || rate < max)
					continue;

				if (pos_Lec > embbedLimit)
					continue;
				if (aln.getPosLec() > lenghtLec * permitedBorderRate
						|| aln.getPosRef() + aln.getLengthRef() < lenghtRef * (1 - permitedBorderRate))
					continue;

				// ref -> lec || ref' -> lec
				overlaps.add(new Overlap(id_Ref, isReverse, id_Lec, false, lenghtRef + pos_Lec, rate));
				break;
			}
		}
	}

	public static void printProgressBar(PrintStream pr, String prev, double act, double max, String post) {
		StringBuilder string = new StringBuilder();
		string.append(prev);
		string.append("[");
		int N = 70, i = 0;
		while (i < N * (act / max)) {
			string.append("#");
			i++;
		}
		while (i < N) {
			string.append(" ");
			i++;
		}
		string.append("] ");
		DecimalFormat df = new DecimalFormat("00.00");
		string.append(df.format(100 * act / max));
		string.append("%");
		string.append(post);
		string.append("\r");
		pr.print(string);
	}

	public static void printProgressBar(PrintStream pr, double act, double max, String post) {
		printProgressBar(pr, "", act, max, post);
	}

	public static void printProgressBar(PrintStream pr, String prev, double act, double max) {
		printProgressBar(pr, prev, act, max, "");
	}

	public static void printProgressBar(PrintStream pr, double act, double max) {
		printProgressBar(pr, "", act, max, "");
	}

	private AssemblyGraph buidAssemblyGraph() {
		int[] map = new int[sequences.size()];
		Arrays.fill(map, 1);
		Queue<Integer> queue = new PriorityQueue<Integer>((Integer a, Integer b) -> a - b);
		for (int i : embeddedOverlaps.keySet()) {
			queue.add(i);
			map[i] = 0;
		}

		map[0] = 0;
		for (int i = 1; i < map.length; i++)
			map[i] += map[i - 1];

		List<CharSequence> list = new LinkedList<CharSequence>();

		int i = 0;
		Integer t;
		while (!queue.isEmpty()) {
			t = queue.poll();
			while (i < t && i < sequences.size()) {
				list.add(sequences.get(i));
				i++;
			}
			i++;
		}
		while (i < sequences.size()) {
			list.add(sequences.get(i));
			i++;
		}

		System.out.println("		Sequences: " + list.size());
		System.out.println("		Emmbeded Sequences: " + embeddedOverlaps.size());

		AssemblyGraph assemblyGraph = new AssemblyGraph(list);
		for (Entry<Integer, Embedded> entry : embeddedOverlaps.entrySet()) {
			Embedded emb = entry.getValue();
			int embeddedId = entry.getKey();
			int parentId = emb.getIdSequence();
			int parentPos = emb.getPos();
			boolean isReverse = emb.isReversed();
			AssemblyEmbedded embedded = new AssemblyEmbedded(sequences.get(embeddedId), parentPos, isReverse);
			assemblyGraph.addEmbedded(map[parentId], embedded);
		}

		int count = 0;
		for (Overlap overlap : overlaps) {
			if (!embeddedOverlaps.containsKey(overlap.getIdFrom())
					&& !embeddedOverlaps.containsKey(overlap.getIdTo())) {
				boolean startFrom = overlap.isFromReversed();
				boolean startTo = !overlap.isToReversed();
				AssemblyVertex from = assemblyGraph.getVertex(map[overlap.getIdFrom()], startFrom);
				AssemblyVertex to = assemblyGraph.getVertex(map[overlap.getIdTo()], startTo);
				assemblyGraph.addEdge(to, from, overlap.getLength());
				count++;
			}
		}

		System.out.println("		Overlaps: " + count + "/" + overlaps.size());

		return assemblyGraph;
	}
}

class Alignment {
	private int posRef;
	private int posLec;
	private int lengthLec;
	private int lengthRef;
	private int hits;

	public Alignment(int posRef, int posLec, int lengthLec, int lengthRef) {
		this.posRef = posRef;
		this.posLec = posLec;
		this.lengthLec = lengthLec;
		this.lengthRef = lengthRef;
		this.hits = 1;
	}

	public void addhit() {
		hits++;
	}

	public int getHits() {
		return hits;
	}

	public void setHits(int hits) {
		this.hits = hits;
	}

	public int getPosRef() {
		return posRef;
	}

	public void setPosRef(int posRed) {
		this.posRef = posRed;
	}

	public int getPosLec() {
		return posLec;
	}

	public void setPosLec(int posLec) {
		this.posLec = posLec;
	}

	public int getLengthLec() {
		return lengthLec;
	}

	public void setLengthLec(int length) {
		this.lengthLec = length;
	}

	public int getLengthRef() {
		return lengthRef;
	}

	public void setLengthRef(int lengthRef) {
		this.lengthRef = lengthRef;
	}
}

class Overlap {
	private int idFrom;
	private boolean fromReversed;
	private int idTo;
	private boolean toReversed;
	private int length;
	private double rate;

	public Overlap(int idFrom, boolean fromReversed, int idTo, boolean toReversed, int length, double rate) {
		this.idFrom = idFrom;
		this.fromReversed = fromReversed;
		this.idTo = idTo;
		this.toReversed = toReversed;
		this.length = length;
		this.rate = rate;
	}

	public double getRate() {
		return rate;
	}

	public void setRate(double rate) {
		this.rate = rate;
	}

	public int getIdFrom() {
		return idFrom;
	}

	public void setIdFrom(int idFrom) {
		this.idFrom = idFrom;
	}

	public boolean isFromReversed() {
		return fromReversed;
	}

	public void setFromReversed(boolean fromReversed) {
		this.fromReversed = fromReversed;
	}

	public int getIdTo() {
		return idTo;
	}

	public void setIdTo(int idTo) {
		this.idTo = idTo;
	}

	public boolean isToReversed() {
		return toReversed;
	}

	public void setToReversed(boolean toReversed) {
		this.toReversed = toReversed;
	}

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

}

class Embedded {
	private int idSequence;
	private int pos;
	private boolean reversed;
	private double rate;

	public Embedded(int idSequence, int pos, boolean reversed, double rate) {
		this.idSequence = idSequence;
		this.pos = pos;
		this.reversed = reversed;
		this.rate = rate;
	}

	public double getRate() {
		return rate;
	}

	public void setRate(double rate) {
		this.rate = rate;
	}

	public int getIdSequence() {
		return idSequence;
	}

	public void setIdSequence(int idSequence) {
		this.idSequence = idSequence;
	}

	public int getPos() {
		return pos;
	}

	public void setPos(int pos) {
		this.pos = pos;
	}

	public boolean isReversed() {
		return reversed;
	}

	public void setReversed(boolean reversed) {
		this.reversed = reversed;
	}

}
