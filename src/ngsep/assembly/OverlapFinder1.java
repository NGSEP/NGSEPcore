package ngsep.assembly;

import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.FMIndex;

public class OverlapFinder1 implements OverlapFinder {
	private Map<Integer, TreeMap<Integer, Alignment>> seqHits;
	private List<CharSequence> sequences;
	private AssemblyGraph assemblyGraph;
	private FMIndex index;

	private Map<Integer, List<AssemblyEdge>> overlapsForward = new Hashtable<>();
	private Map<Integer, Integer> embeddedOverlaps = new Hashtable<>();

	@Override
	public void calculate(List<CharSequence> seq, FMIndex index) {
		assemblyGraph = new AssemblyGraph(seq);
		this.index = index;
		sequences = seq;
		findOverlaps();
	}

	@Override
	public AssemblyGraph getGrap() {
		return assemblyGraph;
	}

	public void findOverlaps() {
		seqHits = new HashMap<>(sequences.size());
		for (int seqId = 0; seqId < sequences.size(); seqId++) {
			if (embeddedOverlaps.containsKey(seqId))
				continue;
			
			seqHits.clear();
			findAlignments(seqId, KmerIterator.NORMAL.iter(sequences.get(seqId)));
			detectOverlap(seqId, false);
			
			seqHits.clear();
			findAlignments(seqId, KmerIterator.REVERSE.iter(sequences.get(seqId)));
			detectOverlap(seqId, true);
		}
	}

	/**
	 * find the hints for CharSequence
	 * 
	 * @param refSeqId
	 * @param read
	 * @return a map whit all the hits per sequence.
	 */
	private void findAlignments(int refSeqId, Iterable<Entry<Integer, String>> kmerIter) {
		for (Entry<Integer, String> entry : kmerIter) {
			int posRef = entry.getKey();
			for (ReadAlignment aln : index.search(entry.getValue())) {
				int lectSeqId = Integer.parseInt(aln.getSequenceName());
				if (lectSeqId > refSeqId || embeddedOverlaps.containsKey(lectSeqId)) {
					TreeMap<Integer, Alignment> treeMap = seqHits.computeIfAbsent(lectSeqId,
							key -> new TreeMap<Integer, Alignment>());

					int posLect = aln.getFirst();
					int key = posRef - posLect, removeKey = key;
					int min = KmerIterator.MAX_KMER_DES;
					Integer celingKey = treeMap.ceilingKey(posRef - posLect),
							floorKey = treeMap.floorKey(posRef - posLect);

					if (celingKey != null && celingKey - key < min) {
						min = celingKey - key;
						removeKey = celingKey;
					}
					if (floorKey != null && key - floorKey < min) {
						min = key - floorKey;
						removeKey = floorKey;
					}

					int length = KmerIterator.SEARCH_KMER_LENGTH
							+ (min != KmerIterator.MAX_KMER_DES ? posLect - treeMap.remove(removeKey).getPosLec() : 0);
					treeMap.put(key, new Alignment(posRef, posLect, length));
				}
			}
		}
	}

	private void detectOverlap(int refSeqId, boolean isReverse) {
		for (Entry<Integer, TreeMap<Integer, Alignment>> entry : seqHits.entrySet()) {
			int lectSeqId = entry.getKey();
			for (Alignment aln : entry.getValue().values()) {
				if (aln.getLength() > KmerIterator.SEARCH_KMER_LENGTH) {
					
				}
			}
		}
	}
}

class Alignment {
	private int posRef;
	private int posLec;
	private int length;

	public Alignment(int posRef, int posLec, int length) {
		this.posRef = posRef;
		this.posLec = posLec;
		this.length = length;
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

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	@Override
	public String toString() {
		return "[posRef=" + posRef + ", posLec=" + posLec + ", length=" + length + "]";
	}
}
