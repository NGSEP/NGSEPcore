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
	private int minLength = 100;
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
		for(List<TransposableElementAnnotation> answP1:answP) {
			if(answP1!=null) answer.addAll(answP1);
		}
		Collections.sort(answer,GenomicRegionPositionComparator.getInstance());
		log.info("Processed sequence "+seq.getName()+". Number of TIRs: "+answer.size());
		return answer;
	}
	private void findTIRsProcess(QualifiedSequence seq, int start, int end, List<List<TransposableElementAnnotation>> answP, int n) {
		CharSequence seqDNA = seq.getCharacters().subSequence(start, end);
		CharSequence rc = DNAMaskedSequence.getReverseComplement(seqDNA);
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		Map<Integer,Long> kmersMapForward = KmersExtractor.extractDNAKmerCodesAsMap(seqDNA, kmerLength, 0, seqDNA.length(),true);
		Map<Long, List<Integer>> reverseMapF = getReverseMap(kmersMapForward);
		Map<Integer,Long> kmersMapReverse = KmersExtractor.extractDNAKmerCodesAsMap(rc, kmerLength, 0, rc.length(),true);
		Map<Long, List<Integer>> reverseMapR = getReverseMap(kmersMapReverse);
		PairwiseAligner pwa = new PairwiseAlignerSimpleGap(minLength);
		for(Map.Entry<Integer, Long> entry:kmersMapReverse.entrySet()) {
			List<Integer> posKmerForward = reverseMapF.get(entry.getValue());
			if(posKmerForward==null) continue;
			int startTIR = posKmerForward.get(posKmerForward.size()-1);
			List<Integer> posRevList = reverseMapR.get(entry.getValue()); 
			int endTIR = seqDNA.length() - posRevList.get(posRevList.size()-1);
			if(endTIR<startTIR+minLength) continue;
			if(startTIR == 591) System.out.println("SeqLen: "+seqDNA.length()+" Kmer: "+new String( AbstractLimitedSequence.getSequence(entry.getValue(), kmerLength, new DNASequence()))+" posForward: "+posKmerForward+" posReverse: "+posRevList+" limits TIR: "+startTIR+" - "+endTIR);
			String candidateTIR = seqDNA.subSequence(startTIR, endTIR).toString();
			int [] tirInfo = validateTIR(candidateTIR,pwa);
			if(startTIR == 591) System.out.println("Seq: "+candidateTIR+" info: "+tirInfo[0]+" "+tirInfo[1]+" "+tirInfo[2]);
			if(tirInfo[0]==1) {
				TransposableElementAnnotation tirCandidate = new TransposableElementAnnotation(seq.getName(), start + startTIR+1, start+endTIR);
				tirCandidate.setRepeatLimits(start+tirInfo[1], start+tirInfo[2], (byte)0);
				answer.add(tirCandidate);
			}
		}
		answP.set(n, answer);
	}
	private int[] validateTIR(String candidateTIR, PairwiseAligner pwa) {
		int[] info = new int [3];
		int n = candidateTIR.length();
		info[0]=0;
		if(n>=minLength) 
		{
			// TODO: Expand and validate other filters
			int end = minLength/2;
			String leftSegment = candidateTIR.substring(0, end);
			String rightSegment = DNAMaskedSequence.getReverseComplement(candidateTIR.substring(n-end, n)).toString();
			PairwiseAlignment alnAfterHit = pwa.calculateAlignment(leftSegment, rightSegment);
			int internalLeft = alnAfterHit.getEnd1();
			int internalRight = n-alnAfterHit.getEnd2();
			if(internalLeft+minLength<internalRight) {
				info[0]=1;
				info[1] = internalLeft;
				info[2] = candidateTIR.length()-internalRight;
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
