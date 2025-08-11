package ngsep.transposons;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.alignments.PairwiseAlignment;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.ReadAlignment.Platform;
import ngsep.alignments.ReadsAligner;
import ngsep.genome.ReferenceGenome;
import ngsep.main.ProgressNotifier;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class DeNovoTransposableElementsFinderConservedEnds implements DeNovoTransposableElementsFinder {

	// Logging and progress
	private Logger log = Logger.getAnonymousLogger();
	private ProgressNotifier progressNotifier=null;
	
	private int debugPos = -1;
	
	private TransposableElementFamily filterOrder = TransposableElementFamily.LTR_UNKNOWN;
	private int subseqLength = 500;
	private int step = 250;
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

	public int getSubseqLength() {
		return subseqLength;
	}
	public void setSubseqLength(int subseqLength) {
		this.subseqLength = subseqLength;
	}

	public int getStep() {
		return step;
	}
	public void setStep(int step) {
		this.step = step;
	}

	@Override
	public List<TransposableElementAnnotation> findTransposons(ReferenceGenome genome) {
		ReadsAligner aligner = new ReadsAligner(genome,Platform.PACBIO);
		aligner.setMaxAlnsPerRead(100);
		ThreadPoolManager pool = new ThreadPoolManager(numThreads, numThreads);
		int sequenceLengthProcess = 2000000;
		pool.setSecondsPerTask(sequenceLengthProcess/1000);
		QualifiedSequenceList seqs = genome.getSequencesList();
		List<List<TransposableElementAnnotation>> answP = new ArrayList<List<TransposableElementAnnotation>>();
		for(QualifiedSequence seq:seqs) {
			int nL = seq.getLength();
			for(int i=0;i<nL;i+=sequenceLengthProcess) {
				final int start = Math.max(0, i-step);
				final int end = Math.min(nL, i+sequenceLengthProcess);
				final int n = answP.size();
				answP.add(null);
				try {
					pool.queueTask(()->findTransposonsProcess(seq, start, end, aligner, answP, n));
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
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		for(List<TransposableElementAnnotation> answP1:answP) {
			if(answP1!=null) answer.addAll(answP1);
		}
		Collections.sort(answer,genome.getComparator());
		return answer;
	}
	private void findTransposonsProcess(QualifiedSequence seq, int start, int end, ReadsAligner aligner, List<List<TransposableElementAnnotation>> answP, int n) {
		answP.set(n,findTransposons(seq, start, end, aligner));
	}
	private List<TransposableElementAnnotation> findTransposons(QualifiedSequence seq, int start, int end, ReadsAligner aligner) {
		log.info("Processing sequence. "+seq.getName()+" from "+start+" to "+end);
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		PairwiseAlignerSimpleGap pwa = new PairwiseAlignerSimpleGap(2000);
		pwa.setForceEnd1(false);
		pwa.setForceEnd2(false);
		HMMTransposonDomainsFinder domainsFinder = new HMMTransposonDomainsFinder();
		domainsFinder.loadHMMsFromClasspath();
		DNAMaskedSequence dna = (DNAMaskedSequence)seq.getCharacters();
		
		for(int i=start;i<end-subseqLength;i+=step) {
			DNAMaskedSequence read = (DNAMaskedSequence)dna.subSequence(i,i+subseqLength);
			List<ReadAlignment> alns = aligner.alignRead(new QualifiedSequence("",read));
			List<ReadAlignment> alnsSeq = selectAlnsInSequence(seq.getName(),alns);
			ReadAlignment aln = selectCloseAlignment(i,subseqLength,alnsSeq);
			if(aln==null) continue;
			//Discard possible hits due to tandem repeats
			int countB = countAlnsInBetween(i+subseqLength,aln.getFirst(),alnsSeq);
			if(countB>1) continue;
			TransposableElementAnnotation ltrAnn = inferTEFromEndAlignment(seq,i,aln, pwa, domainsFinder);
			if(ltrAnn!=null) {
				answer.add(ltrAnn);
				i=ltrAnn.getLast()-step;
			}
		}
		return answer;
	}

	private List<ReadAlignment> selectAlnsInSequence(String name, List<ReadAlignment> alns) {
		List<ReadAlignment> alnsSeq = new ArrayList<ReadAlignment>();
		for(ReadAlignment aln:alns) {
			if(name == aln.getSequenceName()) alnsSeq.add(aln);
		}
		return alnsSeq;
	}

	private ReadAlignment selectCloseAlignment(int start, int length,  List<ReadAlignment> alns) {
		ReadAlignment bestAln = null;
		for(ReadAlignment aln:alns) {
			if(aln.getFirst()<start+length) continue;
			if(aln.getFirst()>start+30000) continue;
			if(improvedAln(start, bestAln, aln)) bestAln = aln;
		}
		return bestAln;
	}

	private boolean improvedAln(int start, ReadAlignment bestAln, ReadAlignment aln) {
		if(bestAln==null) return true;
		int distanceBest = bestAln.getFirst()-start;
		int distanceAln = aln.getFirst()-start;
		int diff5000Best = Math.abs(distanceBest-5000);
		int diff5000Aln = Math.abs(distanceAln-5000);
		return diff5000Aln<diff5000Best;
	}
	
	private int countAlnsInBetween(int start, int end, List<ReadAlignment> alns) {
		int count = 0;
		for(ReadAlignment aln:alns) {
			if(start<=aln.getFirst() && aln.getLast()<=end) count++; 
		}
		return count;
	}

	private TransposableElementAnnotation inferTEFromEndAlignment(QualifiedSequence seq, int start, ReadAlignment aln, PairwiseAlignerSimpleGap pwa, HMMTransposonDomainsFinder domainsFinder) {
		System.err.println("Checking hit of sequence. Start: "+start+" aln: "+aln);
		//Find start alignment before hit
		DNAMaskedSequence dna = (DNAMaskedSequence) seq.getCharacters();
		DNAMaskedSequence leftSegment = ((DNAMaskedSequence)dna.subSequence(Math.max(0, start-300), start)).getReverseComplement();
		DNAMaskedSequence rightSegment;
		if(aln.isPositiveStrand()) rightSegment = ((DNAMaskedSequence)dna.subSequence(Math.max(0, aln.getFirst()-300), aln.getFirst())).getReverseComplement();
		else rightSegment = (DNAMaskedSequence)dna.subSequence(aln.getLast(), Math.min(seq.getLength(), aln.getLast()+300));
		
		
		PairwiseAlignment alnBeforeHit = pwa.calculateAlignment(leftSegment, rightSegment);
		int start2 = Math.max(0, start - alnBeforeHit.getEnd1());
		//Updated later if same strand
		int end2 = Math.min(seq.getLength(), aln.getLast()+alnBeforeHit.getEnd2());
		//Updated later if opposite strand
		int internalRight = Math.max(0, aln.getFirst()-alnBeforeHit.getEnd2());
		
		//Find end alignment after hit
		leftSegment = (DNAMaskedSequence)dna.subSequence(start, Math.min(seq.getLength(), start+2000));
		if(aln.isPositiveStrand()) rightSegment = (DNAMaskedSequence)dna.subSequence(aln.getFirst()-1, Math.min(seq.getLength(), aln.getFirst()+1999));
		else rightSegment = ((DNAMaskedSequence)dna.subSequence(Math.max(0, aln.getLast()-2000), aln.getLast())).getReverseComplement();
		PairwiseAlignment alnAfterHit = pwa.calculateAlignment(leftSegment, rightSegment);
		int internalLeft = start+alnAfterHit.getEnd1();
		if(aln.isPositiveStrand()) end2 = Math.min(seq.getLength(), aln.getFirst()+alnAfterHit.getEnd2());
		else internalRight = Math.max(0, aln.getLast()-alnAfterHit.getEnd2());
		
		System.err.println("Checked borders. Start: "+start2+" internal left "+internalLeft+" internal right: "+internalRight+" end "+end2);
		if(internalLeft<start2+aln.getReadLength() || internalLeft>end2) return null;
		if(internalRight<internalLeft || internalRight>end2-aln.getReadLength()) return null;
		
		//Create annotation
		TransposableElementAnnotation ann = new TransposableElementAnnotation(seq.getName(), start2+1, end2);
		ann.setRepeatLimits(internalLeft+1, internalRight+1, aln.isPositiveStrand()?TransposableElementFamily.REPEAT_ORIENTATION_FF:TransposableElementFamily.REPEAT_ORIENTATION_FR);
		
		//Assign family
		DNAMaskedSequence dnaTE = (DNAMaskedSequence) seq.getCharacters().subSequence(internalLeft, internalRight);
		domainsFinder.assignFamily(ann,dnaTE);
		if(ann.getInferredFamily()==null) {
			domainsFinder.assignFamily(ann,dnaTE.getReverseComplement());
			if(ann.getInferredFamily()!=null) ann.setNegativeStrand(true);
		}
		if(!passFilters(ann)) return null;
		return ann;
	}
	
	
	private boolean passFilters(TransposableElementAnnotation ann) {
		if(filterOrder!=null) {
			TransposableElementFamily family = ann.getInferredFamily();
			if(family==null || !family.getOrder().equals(filterOrder.getOrder())) return false;
		}
		return true;
	}

	public static void main(String[] args) throws Exception {
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		DeNovoTransposableElementsFinderConservedEnds instance = new DeNovoTransposableElementsFinderConservedEnds();
		instance.numThreads = 4;
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
