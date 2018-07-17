package ngsep.assembly;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndexSingleSequence;

/**
 * this class contains the fmIndex of each sequence
 * 
 * @author jc.bojaca
 *
 */
public class FMIndexSequences implements Serializable {
    /** serialVersionUID */
    private static final long serialVersionUID = 1L;
    /** FMIndex array */
    private final FMIndexSingleSequence[] sequences;
    /** FMIndex array for the complementary sequences */
    private final FMIndexSingleSequence[] complementarySequences;
    /** number of sequences */
    private final int N;

    public FMIndexSequences(List<? extends DNAMaskedSequence> sequences) {
	this(sequences, true);
    }

    public FMIndexSequences(List<? extends DNAMaskedSequence> sequences, boolean parallel) {
	this.N = sequences.size();
	this.sequences = new FMIndexSingleSequence[N];
	this.complementarySequences = new FMIndexSingleSequence[N];

	if (parallel) {
	    try {
		loadParallel(Runtime.getRuntime().availableProcessors(), sequences);
	    } catch (InterruptedException e) {
		e.printStackTrace();
	    }
	} else {
	    loadSerial(sequences);
	}
    }

    private void loadSerial(List<? extends DNAMaskedSequence> sequences) {
	Thread act = new FMIndexSequence(sequences, this.sequences, this.complementarySequences, 0, 1, N);
	act.run();
    }

    private void loadParallel(int nThreads, List<? extends DNAMaskedSequence> sequences) throws InterruptedException {
	System.out.println("# of threads: " + nThreads);
	Thread act;
	Queue<Thread> threads = new LinkedList<>();

	for (int i = 0; i < nThreads; i++) {
	    act = new FMIndexSequence(sequences, this.sequences, this.complementarySequences, i, nThreads, N);
	    act.start();
	    threads.add(act);
	}
	while (!threads.isEmpty())
	    threads.poll().join();
    }

    public Set<Integer> search(String kmer, int posSequnce, boolean complementary) {
	FMIndexSingleSequence[] array = complementary ? complementarySequences : sequences;
	return array[posSequnce].search(kmer);
    }

    class FMIndexSequence extends Thread {
	/** original sequences */
	private final List<? extends DNAMaskedSequence> sequences;
	/** FMIndex array (only write) */
	private final FMIndexSingleSequence[] index;
	/** FMIndex array for the complementary sequences */
	private final FMIndexSingleSequence[] complementarySequences;
	/** Id */
	private final int id;
	/** # of threads */
	private final int Nthreads;
	/** number of sequences */
	private final int N;

	public FMIndexSequence(List<? extends DNAMaskedSequence> sequences, FMIndexSingleSequence[] index,
		FMIndexSingleSequence[] complementarySequences, int id, int nthreads, int n) {
	    this.sequences = sequences;
	    this.complementarySequences = complementarySequences;
	    this.index = index;
	    this.id = id;
	    Nthreads = nthreads;
	    N = n;
	}

	@Override
	public void run() {
	    DNAMaskedSequence aux;
	    for (int i = id; i < N; i += Nthreads) {
		aux = sequences.get(i);
		index[i] = new FMIndexSingleSequence(aux);
		complementarySequences[i] = new FMIndexSingleSequence(aux.getReverseComplement());
	    }
	}
    }
}
