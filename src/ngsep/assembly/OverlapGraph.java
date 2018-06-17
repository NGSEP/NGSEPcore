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
package ngsep.assembly;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.stream.Stream;
import ngsep.alignments.ReadAlignment;
import ngsep.math.Distribution;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.LimitedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.sequences.io.FastqFileReader;

/**
 * @author Jorge Duitama && jc.bojaca
 *
 */
public class OverlapGraph {
    static final int SEARCH_KMER_LENGTH = 15;
    static final int SEARCH_KMER_DISTANCE = 10;
    static final int SEARCH_KMER_OVERLAP = SEARCH_KMER_LENGTH - SEARCH_KMER_DISTANCE;
    private static final String[] fastq = { ".fastq", ".fastq.gz" };
    private static final String[] fasta = { ".fasta" };

    private final List<LimitedSequence> sequences = new ArrayList<>();
    private FMIndex index;

    private Map<Integer, List<ReadOverlap>> overlapsForward = new Hashtable<>();
    private Map<Integer, List<ReadOverlap>> overlapsBackward = new Hashtable<>();
    private Map<Integer, ReadOverlap> embeddedOverlaps = new Hashtable<>();

    public OverlapGraph(String filename) throws Exception {
	System.out.println("Built the Overlap Grap for:  " + filename + " ...\n");
	timeWithException("total", () -> {
	    timeWithException("Loaded the sequences from " + filename + ".", () -> {
		load(filename);
		System.out.println("Loaded " + sequences.size() + " sequences");
	    });

	    time("Sort the sequences.", () -> {
		Collections.sort(sequences, (LimitedSequence l1, LimitedSequence l2) -> l2.length() - l1.length());
	    });

	    time("Built the FMIndex.", () -> {
		index = new FMIndex();
		index.loadUnnamedSequences("", sequences);
	    });

	    time("Built overlaps.", () -> findOverlaps2());

	    time("Simplify the grap.", () -> {
		searchSequences();
	    });
	});

	// printOverlapsDistribution(System.out);
    }

    public void load(String filename) throws IOException {
	if (Stream.of(fastq)
		.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase()))) {
	    loadFastq(filename);
	} else if (Stream.of(fasta)
		.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase()))) {
	    loadFasta(filename);
	} else
	    throw new IOException("the file not is a fasta or fastq file: " + filename);
    }

    private void loadFasta(String filename) throws IOException {
	FastaSequencesHandler handler = new FastaSequencesHandler();
	QualifiedSequenceList seqsQl = handler.loadSequences(filename);
	for (QualifiedSequence seq : seqsQl) {
	    DNAMaskedSequence characters = (DNAMaskedSequence) seq.getCharacters();
	    sequences.add(characters);
	}
    }

    private void loadFastq(String filename) throws IOException {
	try (FastqFileReader reader = new FastqFileReader(filename)) {
	    Iterator<RawRead> it = reader.iterator();
	    while (it.hasNext()) {
		RawRead read = it.next();
		DNAMaskedSequence characters = (DNAMaskedSequence) read.getCharacters();
		sequences.add(characters);
	    }
	}
    }

    @SuppressWarnings("unchecked")
    public void findOverlaps2() {
	List<ReadOverlap2>[] procesingHits = new List[sequences.size()];
	Queue<ReadOverlap2>[] hits = new Queue[sequences.size()];

	for (int idSequence = 0; idSequence < sequences.size(); idSequence++) {
	    System.out.println("  --> "+idSequence + " - " + sequences.size());
	    long ini = System.nanoTime();
	    if (embeddedOverlaps.containsKey(idSequence))
		continue;
	    findHits(idSequence, procesingHits, hits);

	    for (int i = 0; i < hits.length; i++) {
		Queue<ReadOverlap2> queue = hits[i];
		while (!queue.isEmpty())
		    if (processOverlap(idSequence, i, queue.poll()))
			break;
	    }
	    System.out.println(System.nanoTime() - ini);
	}
    }

    private boolean processOverlap(int idSequence1, int idSequence2, ReadOverlap2 overlap) {
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
	// overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), true);
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
	// overlap.getIndexSequence1(), overlap.getFirst1(), overlap.getLast1(), true);
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
	    addOverlap(overlapsBackward, overlapp);
	    ReadOverlap ov2 = new ReadOverlap(idSequence2, overlap.getFirst2(), overlap.getLast2(), idSequence1,
		    overlap.getFirst1(), overlap.getLast1(), false);
	    addOverlap(overlapsForward, ov2);
	    return true;
	}
	if (dl1 > dl2 && dr1 < dr2) {
	    // 1 left of 2
	    if (dl2 > diffBorder || dr1 > diffBorder) {
		// Not a true overlap
		return false;
	    }
	    addOverlap(overlapsForward, overlapp);
	    ReadOverlap ov2 = new ReadOverlap(idSequence2, overlap.getFirst2(), overlap.getLast2(), idSequence1,
		    overlap.getFirst1(), overlap.getLast1(), false);
	    addOverlap(overlapsBackward, ov2);
	    return true;
	}
	// }
	return false;
    }

    private void findHits(int idSequence, List<ReadOverlap2>[] procesingHits, Queue<ReadOverlap2>[] hits) {
	for (int i = 0; i < procesingHits.length; i++)
	    // TODO optimizar (array), limpiar colas, cambio en las estructuras para
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
	    sum += processKmer(idSequence, i, sequence.substring(i, i + SEARCH_KMER_LENGTH), procesingHits);
	    finalizeHits(procesingHits, hits);
	}
	if (sequence.length() % SEARCH_KMER_LENGTH != 0) {
	    int KmerStarPosition = ((sequence.length() - SEARCH_KMER_OVERLAP) / SEARCH_KMER_DISTANCE)
		    * SEARCH_KMER_DISTANCE;
	    InitializeHits(procesingHits);
	    sum += processKmer(idSequence, KmerStarPosition, sequence.substring(KmerStarPosition), procesingHits);
	    finalizeHits(procesingHits, hits);
	}
	for (int i = 0; i < procesingHits.length; i++)
	    for (ReadOverlap2 read : procesingHits[i])
		hits[i].add(read);
	System.out.println(sum);

    }

    private void InitializeHits(List<ReadOverlap2>[] procesingHits) {
	for (List<ReadOverlap2> set : procesingHits)
	    for (ReadOverlap2 read : set)
		read.nonAdded();
    }

    private long processKmer(int idSequence, int KmerStarPosition, String kmer, List<ReadOverlap2>[] procesingHits) {
	long ini = System.nanoTime();
	Iterable<ReadAlignment> itre = index.search(kmer);
	ini = System.nanoTime() - ini;
	for (ReadAlignment aln : itre) {
	    int idSequenceAligned = Integer.parseInt(aln.getSequenceName().substring(1));
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
			procesingHits[idSequenceAligned]
				.add(new ReadOverlap2(KmerStarPosition, aln.getFirst() - 1, aln.getLast()));
		} else {
		    // TODO caso negativo
		    System.out.println("ouchh");
		    continue;
		}
	    }
	}
	return ini;
    }

    private void finalizeHits(List<ReadOverlap2>[] procesingHits, Queue<ReadOverlap2>[] hits) {
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

    public void findOverlaps() {
	for (int s = 0; s < sequences.size(); s++) {
	    if (embeddedOverlaps.containsKey(s))
		continue;

	    Map<Integer, List<ReadAlignment>> seqHits = findHitsRead(s);
	    // System.out.println(seqHits.size());
	    // System.out.println(seqHits.keySet());
	    buildOverlapsFromKmerAlignments(s, seqHits);
	    // System.out.println("overlapsBackward");
	    // for(ReadOverlap read:overlapsBackward.get(s)) {
	    // System.out.println(read.getIndexSequence1()+","+read.getIndexSequence2());
	    // }
	    //
	    // System.out.println("overlapsForward");
	    // for(ReadOverlap read:overlapsForward.get(s)) {
	    // System.out.println(read.getIndexSequence1()+","+read.getIndexSequence2());
	    // }
	    // System.out.println(embeddedOverlaps);
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
		int k = Integer.parseInt(aln.getSequenceName().substring(1));
		if (k > idRead) {
		    seqHits.computeIfAbsent(k, key -> new ArrayList<>()).add(aln);
		}
	    }
	}
	return seqHits;
    }

    private void buildOverlapsFromKmerAlignments(int searchId, Map<Integer, List<ReadAlignment>> seqHits) {
	LimitedSequence searchSequence = sequences.get(searchId);
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

    public List<List<Integer>> searchSequences() {
	int N = sequences.size();
	Map<Integer, Map<Integer, Integer>> Edges = getEdges();
	GraphSimplificator gs = new GraphSimplificator(N, Edges);
	Edges = gs.Edges;
	MaxPath mx = new MaxPath(N, Edges);
	time("sdfasd", () -> {
	    System.out.println(mx.maxLentghPath());
	});
	return null;
    }

    private Map<Integer, Map<Integer, Integer>> getEdges() {
	Map<Integer, Map<Integer, Integer>> ans = new HashMap<>(sequences.size() - embeddedOverlaps.size());
	Map<Integer, Integer> list;
	for (int i = 0; i < sequences.size(); i++)
	    ans.put(i, new HashMap<>());

	for (Entry<Integer, List<ReadOverlap>> entry : overlapsForward.entrySet()) {
	    List<ReadOverlap> overlaps = entry.getValue();
	    Integer vert = entry.getKey();
	    if (embeddedOverlaps.containsKey(vert))
		vert = embeddedOverlaps.get(vert).getIndexSequence2();
	    list = ans.get(vert);
	    for (ReadOverlap over : overlaps) {
		if (!embeddedOverlaps.containsKey(over.getIndexSequence2()))
		    list.put(over.getIndexSequence2(), over.getLast2() - over.getFirst2());
	    }
	}
	return ans;
    }

    public void printOverlapsDistribution(PrintStream out) {
	Distribution distForward = new Distribution(0, 20, 1);
	Distribution distBackward = new Distribution(0, 20, 1);
	int numEmbedded = embeddedOverlaps.size();
	for (int i = 0; i < sequences.size(); i++) {
	    List<ReadOverlap> ovs = overlapsForward.get(i);
	    if (ovs == null)
		distForward.processDatapoint(0);
	    else
		distForward.processDatapoint(ovs.size());
	    ovs = overlapsBackward.get(i);
	    if (ovs == null)
		distBackward.processDatapoint(0);
	    else
		distBackward.processDatapoint(ovs.size());
	}
	double[] vF = distForward.getDistribution();
	double[] vB = distBackward.getDistribution();
	System.out.println("Overlaps distribution: ");
	for (int i = 0; i < vF.length; i++) {
	    out.println("" + i + "\t" + vF[i] + "\t" + vB[i]);
	}
	out.println("More:\t" + distForward.getOutliers().size() + "\t" + distBackward.getOutliers().size());
	out.println();
	out.println("Embedded alignments: " + numEmbedded);

    }

    public static final void time(String string, Runnable runnable) {
	long time = System.currentTimeMillis();
	runnable.run();
	System.out.println(string + " time: " + ((System.currentTimeMillis() - time) / (double) 1000) + "s");
    }

    public static final void timeWithException(String string, RunnableWithException runnable) throws Exception {
	long time = System.currentTimeMillis();
	runnable.run();
	System.out.println(string + " time: " + ((System.currentTimeMillis() - time) / (double) 1000) + "s");
    }

    public static void main(String[] args) throws Exception {
	new OverlapGraph(args[0]);
    }
}

@FunctionalInterface
interface RunnableWithException {
    void run() throws Exception;
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
	int kmerLast = kmerFirst + OverlapGraph.SEARCH_KMER_LENGTH;
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
