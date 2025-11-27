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
			log.info("Processing sequence: "+seq.getName());
			List<TransposableElementAnnotation> answerSeq = findTransposons(seq);
			answer.addAll(answerSeq);
		}
			
		return answer;
	}
	public List<TransposableElementAnnotation> findTransposons(QualifiedSequence seq) {
		ThreadPoolManager pool = new ThreadPoolManager(numThreads, 1000);
		pool.setSecondsPerTask(10);
		//Hashmap indexed by start position
		Map<Integer,List<TransposableElementAnnotation>> answP = new HashMap<Integer,List<TransposableElementAnnotation>>();
		int nL = seq.getLength();
		int overlappingSegments = windowLength/step;
		if(windowLength%step !=0) overlappingSegments++;
		log.info("Window length: "+windowLength+" step: "+step+" overlapping segments: "+overlappingSegments);
		for(int s=0;s<overlappingSegments;s++) {
			int startSteps = step*s;
			for(int i=startSteps;i<nL;i+=windowLength) {
				final int start = i;
				final int end = Math.min(nL, i+windowLength);
				try {
					pool.queueTask(()->findTransposonsProcess(seq, start, end, answP));
				} catch (InterruptedException e) {
					throw new RuntimeException("Process interrupted", e);
				}
			}
		}
		try {
			pool.terminatePool();
		} catch (InterruptedException e) {
			throw new RuntimeException("Process interrupted", e);
		}
		List<TransposableElementAnnotation> allEvents = new ArrayList<TransposableElementAnnotation>();
		for(List<TransposableElementAnnotation> answP1:answP.values()) {
			if(answP1!=null) allEvents.addAll(answP1);
		}
		Collections.sort(allEvents,GenomicRegionPositionComparator.getInstance());
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		TransposableElementAnnotation lastAnn = null;
		for(TransposableElementAnnotation ann:allEvents) {
			if(!passFilters(ann)) continue;
			if(lastAnn==null || lastAnn.getLast()<ann.getFirst()) {
				answer.add(ann);
				lastAnn = ann;
			}
		}
		return answer;
	}
	private void findTransposonsProcess(QualifiedSequence seq, int start, int end, Map<Integer,List<TransposableElementAnnotation>> answP) {
		List<TransposableElementAnnotation> overlappingEvents = findOverlappingEvents(seq, start,end, answP);
		answP.put(start,findTransposons(seq, start, end, overlappingEvents));
	}
	private List<TransposableElementAnnotation> findOverlappingEvents(QualifiedSequence seq, int start, int end, Map<Integer, List<TransposableElementAnnotation>> answP) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		int overlappingSegments = windowLength/step;
		int s = Math.max(0, start-overlappingSegments*step);
		for(int i=s;i<end;i+=step) {
			List<TransposableElementAnnotation> events = answP.get(i);
			if(events==null) continue;
			for(TransposableElementAnnotation ann:events) {
				if(start <=ann.getLast() && ann.getFirst()<=end) answer.add(ann);
			}
		}
		return answer;
	}
	protected List<TransposableElementAnnotation> findTransposons(QualifiedSequence seq, int start, int end, List<TransposableElementAnnotation> overlappingEvents) {
		return new ArrayList<TransposableElementAnnotation>();
	}
	protected boolean passFilters(TransposableElementAnnotation ann) {
		return true;
	}
}
