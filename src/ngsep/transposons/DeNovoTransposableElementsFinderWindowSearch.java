package ngsep.transposons;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.ProgressNotifier;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class DeNovoTransposableElementsFinderWindowSearch implements DeNovoTransposableElementsFinder {

	// Logging and progress
	private Logger log = Logger.getAnonymousLogger();
	private ProgressNotifier progressNotifier=null;
	
	private int windowLength = 20000;
	private int step = 5000;
	private int kmerLength = 11;
	private int numThreads = 1;
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}

	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}

	
	public int getWindowLength() {
		return windowLength;
	}
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}
	public int getStep() {
		return step;
	}
	public void setStep(int step) {
		this.step = step;
	}
	
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
	}
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	@Override
	public List<TransposableElementAnnotation> findTransposons(ReferenceGenome genome) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		QualifiedSequenceList seqs = genome.getSequencesList();
		for(QualifiedSequence seq:seqs) {
			List<TransposableElementAnnotation> answerSeq = findTransposons(seq);
			answer.addAll(answerSeq);
		}
			
		return answer;
	}
	public List<TransposableElementAnnotation> findTransposons(QualifiedSequence seq) {
		ThreadPoolManager pool = new ThreadPoolManager(numThreads, 1000);
		pool.setSecondsPerTask(10);
		List<List<TransposableElementAnnotation>> answP = new ArrayList<List<TransposableElementAnnotation>>();
		int nL = seq.getLength();
		for(int i=0;i<nL;i+=step) {
			final int start = i;
			final int end = Math.min(nL, i+windowLength);
			final int n = answP.size();
			answP.add(null);
			try {
				pool.queueTask(()->findTransposonsProcess(seq, start, end, answP, n));
			} catch (InterruptedException e) {
				throw new RuntimeException("Process interrupted", e);
			}
		}
		try {
			pool.terminatePool();
		} catch (InterruptedException e) {
			throw new RuntimeException("Process interrupted", e);
		}
		List<TransposableElementAnnotation> allEvents = new ArrayList<TransposableElementAnnotation>();
		for(List<TransposableElementAnnotation> answP1:answP) {
			if(answP1!=null) allEvents.addAll(answP1);
		}
		Collections.sort(allEvents,GenomicRegionPositionComparator.getInstance());
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		TransposableElementAnnotation lastAnn = null;
		for(TransposableElementAnnotation ann:allEvents) {
			if(lastAnn==null || lastAnn.getLast()<ann.getFirst()) {
				answer.add(ann);
				lastAnn = ann;
			}
		}
		return answer;
	}
	private void findTransposonsProcess(QualifiedSequence seq, int start, int end, List<List<TransposableElementAnnotation>> answP, int n) {
		answP.set(n,findTransposons(seq, start, end));
	}
	protected List<TransposableElementAnnotation> findTransposons(QualifiedSequence seq, int start, int end) {
		return new ArrayList<TransposableElementAnnotation>();
	}
	protected Map<Long, List<Integer>> getReverseMap(Map<Integer, Long> kmersMap) {
		Map<Long,List<Integer>> reverseMapF = new HashMap<Long, List<Integer>>();
		for(Map.Entry<Integer, Long> entry:kmersMap.entrySet()) {
			List<Integer> posKmer = reverseMapF.computeIfAbsent(entry.getValue(), v-> new ArrayList<Integer>());
			posKmer.add(entry.getKey());
		}
		return reverseMapF;
	}
}
