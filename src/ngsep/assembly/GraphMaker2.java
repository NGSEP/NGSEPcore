package ngsep.assembly;

import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

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
	Collections.sort(sequences, (LimitedSequence l1, LimitedSequence l2) -> l2.length() - l1.length());
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
	TreeMap<Integer, int[]>[][] a = new TreeMap[sequencesToIndex.size()][KmerEmmbeddedIterator.NUMBER_OF_KMERS];

	for (int idSequence = 0; idSequence < sequencesSize; idSequence++) {
	    if (embeddedOverlaps.containsKey(idSequence))
		continue;
	    System.out.println(idSequence + "->" + sequencesToIndex.size());

	    KmerEmmbeddedIterator ki = new KmerEmmbeddedIterator(sequences.get(idSequence));
	    int k = 0;
	    for (String string : ki.firts()) {
		index.search(string);
		k++;
	    }
	    for (String string : ki.lasts()) {
		index.search(string);
		k++;
	    }

	    System.out.println(k);
	}
	System.out.println("asdljkhjfgdfdkhjgjfdrse: " + (System.currentTimeMillis() - i));
    }

}
