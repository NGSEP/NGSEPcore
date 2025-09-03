package ngsep.transposons;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.PairwiseAligner;
import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.alignments.PairwiseAlignment;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.ProgressNotifier;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.AbstractLimitedSequence;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class DeNovoTransposableElementsFinderShortIR implements DeNovoTransposableElementsFinder {
	// Logging and progress
	private Logger log = Logger.getAnonymousLogger();
	private ProgressNotifier progressNotifier=null;
	
	private int debugPos = -1;

	private int kmerLength = 10;
	private int windowLength = 2000;
	private int minElementLength = 100;
	private int minTIRLength = 12;
	private int step = 1000;
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
	
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
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
		log.info("Processing sequence "+seq.getName());
		ThreadPoolManager pool = new ThreadPoolManager(numThreads, 10*numThreads);
		pool.setSecondsPerTask(1);
		List<List<TransposableElementAnnotation>> answP = new ArrayList<List<TransposableElementAnnotation>>();
		int nL = seq.getLength();
		for(int i=0;i<nL;i+=step) {
			final int start = i;
			final int end = Math.min(nL, i+windowLength);
			final int n = answP.size();
			answP.add(null);
			try {
				pool.queueTask(()->findTIRsProcess(seq, start, end, answP, n));
			} catch (InterruptedException e) {
				throw new RuntimeException("Process interrupted", e);
			}
		}
	
		try {
			pool.terminatePool();
		} catch (InterruptedException e) {
			throw new RuntimeException("Process interrupted", e);
		}
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		TransposableElementAnnotation lastAnn = null;
		for(List<TransposableElementAnnotation> answP1:answP) {
			if(answP1==null) continue;
			for(TransposableElementAnnotation ann:answP1) {
				if(lastAnn==null || lastAnn.getLast()<ann.getFirst()) {
					answer.add(ann);
					lastAnn = ann;
				}
			}
		}
		Collections.sort(answer,GenomicRegionPositionComparator.getInstance());
		log.info("Processed sequence "+seq.getName()+". Number of TIRs: "+answer.size());
		return answer;
	}
	private void findTIRsProcess(QualifiedSequence seq, int start, int end, List<List<TransposableElementAnnotation>> answP, int answPIndex) {
		HMMTransposonDomainsFinder domainsFinder = new HMMTransposonDomainsFinder();
		domainsFinder.loadHMMsFromClasspath();
		CharSequence seqDNA = seq.getCharacters().subSequence(start, end);
		int n = seqDNA.length();
		CharSequence rc = DNAMaskedSequence.getReverseComplement(seqDNA);
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		Map<Integer,Long> kmersMapForward = KmersExtractor.extractDNAKmerCodesAsMap(seqDNA, kmerLength, 0, seqDNA.length(),true);
		Map<Long, List<Integer>> reverseMapF = getReverseMap(kmersMapForward);
		Map<Integer,Long> kmersMapReverse = KmersExtractor.extractDNAKmerCodesAsMap(rc, kmerLength, 0, rc.length(),true);
		Map<Long, List<Integer>> reverseMapR = getReverseMap(kmersMapReverse);
		PairwiseAlignerSimpleGap pwa = new PairwiseAlignerSimpleGap(minElementLength);
		pwa.setForceEnd1(false);
		pwa.setForceEnd2(false);
		int lastCandidateStart = seqDNA.length();
		for(Map.Entry<Integer, Long> entry:kmersMapReverse.entrySet()) {
			if(n-entry.getKey()>lastCandidateStart) continue;
			List<Integer> posKmerForward = reverseMapF.get(entry.getValue());
			if(posKmerForward==null) continue;
			int startTIR = -1;
			for(int i = posKmerForward.size()-1;i>=0;i--) {
				startTIR = posKmerForward.get(i);
				if(startTIR<lastCandidateStart) break;
			}
			if(start == debugPos) System.out.println("Next candidate from: "+startTIR+" to: "+(n-entry.getKey())+" lastCandStart: "+lastCandidateStart);
			if(startTIR>=lastCandidateStart) continue;
			List<Integer> posRevList = reverseMapR.get(entry.getValue()); 
			int endTIR = n - posRevList.get(posRevList.size()-1);
			if(endTIR<startTIR+minElementLength) continue;
			if(start == debugPos) System.out.println("SeqLen: "+seqDNA.length()+" Kmer: "+new String( AbstractLimitedSequence.getSequence(entry.getValue(), kmerLength, new DNASequence()))+" posForward: "+posKmerForward+" posReverse: "+posRevList+" limits TIR: "+startTIR+" - "+endTIR);
			DNAMaskedSequence candidateTIR = (DNAMaskedSequence)seqDNA.subSequence(startTIR, endTIR);
			int [] tirInfo = validateTIR(candidateTIR.toString(),pwa, start == debugPos);
			if(start == debugPos) System.out.println("Seq: "+candidateTIR+" info: "+tirInfo[0]+" "+tirInfo[1]+" "+tirInfo[2]);
			if(tirInfo[0]==1) {
				TransposableElementAnnotation tirCandidate = new TransposableElementAnnotation(seq.getName(), start + startTIR+1, start+endTIR);
				tirCandidate.setRepeatLimits(start+startTIR+tirInfo[1], start+startTIR+tirInfo[2], (byte)0);
				answer.add(tirCandidate);
				lastCandidateStart = startTIR;
				domainsFinder.assignFamily(tirCandidate, candidateTIR);
			}
		}
		//if(start > 28000 && start < 33000) System.out.println("Start window: "+start+" Final candidate regions: "+answer.size());
		answP.set(answPIndex, answer);
	}
	private int[] validateTIR(String candidateTIR, PairwiseAligner pwa, boolean debug) {
		int[] info = new int [3];
		int n = candidateTIR.length();
		info[0]=0;
		if(n>=minElementLength) 
		{
			if(debug) System.out.println("Element length: "+n);
			// TODO: Validate other filters
			int end = Math.max(minElementLength,n/2)/2;
			String leftSegment = candidateTIR.substring(0, end);
			String rightSegment = DNAMaskedSequence.getReverseComplement(candidateTIR.substring(n-end, n)).toString();
			if(debug) {
				System.out.println("L: "+leftSegment);
				System.out.println("R: "+rightSegment);
				
			}
			PairwiseAlignment alnAfterHit = pwa.calculateAlignment(leftSegment, rightSegment);
			if(debug) {
				//pwa.printAlignmentMatrix(pwa.getMatchScores(), leftSegment, rightSegment);
				System.out.println("AlnL: "+alnAfterHit.getAlignedSequence1());
				System.out.println("AlnR: "+alnAfterHit.getAlignedSequence2());
				
			}
			int internalLeft = alnAfterHit.getEnd1();
			int internalRight = n-alnAfterHit.getEnd2();
			int mismatches = alnAfterHit.getMismatches();
			int limit = Math.min(internalLeft, internalRight)/5;
			limit = Math.max(limit, 2);
			if(debug) System.out.println("Limits: "+internalLeft+" "+internalRight+" mismatches: "+mismatches+" limit: "+limit);
			if(internalLeft+minElementLength<internalRight && Math.min(internalLeft, internalRight) >= minTIRLength &&  mismatches<=limit) {
				info[0]=1;
				info[1] = internalLeft;
				info[2] = internalRight;
			}
		}
		return info;
	}
	private Map<Long, List<Integer>> getReverseMap(Map<Integer, Long> kmersMapForward) {
		Map<Long,List<Integer>> reverseMapF = new HashMap<Long, List<Integer>>();
		for(Map.Entry<Integer, Long> entry:kmersMapForward.entrySet()) {
			List<Integer> posKmer = reverseMapF.computeIfAbsent(entry.getValue(), v-> new ArrayList<Integer>());
			posKmer.add(entry.getKey());
		}
		return reverseMapF;
	}
	public static void main(String[] args) throws Exception {
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		DeNovoTransposableElementsFinderShortIR instance = new DeNovoTransposableElementsFinderShortIR();
		instance.numThreads = 1;
		List<TransposableElementAnnotation> anns = instance.findTransposons(genome);
		try (PrintStream out=new PrintStream(args[1])) {
			for(TransposableElementAnnotation ann:anns) {
				out.print(""+ann.getSequenceName()+"\t"+ann.getFirst()+"\t"+ann.getLast()+"\t"+ann.isPositiveStrand()+"\t"+ann.getInferredFamily()+"\t"+ann.getLeftEndRepeat()+"\t"+ann.getRightStartRepeat()+"\t"+ann.getOrientation());
				if(ann.getDomainAlignments()!=null && ann.getDomainAlignments().size()>0) {
					out.print("\t");
					for(TransposonDomainAlignment daln:ann.getDomainAlignments()) {
						out.print(daln.getDomainCode()+":"+daln.getStart()+":"+daln.getEvalue()+";");
					}
				}
				out.println();
			}
		}	
	}
}
