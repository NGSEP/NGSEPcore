package ngsep.transposons;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.alignments.PairwiseAligner;
import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.alignments.PairwiseAlignment;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.AbstractLimitedSequence;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;

public class DeNovoTransposableElementsFinderShortIR extends DeNovoTransposableElementsFinderWindowSearch implements DeNovoTransposableElementsFinder {
	
	private int debugPos = -1;

	private int minElementLength = 100;
	private int minTIRLength = 12;
	private boolean keepTIRsWithoutDomains = false;
	
	private HMMTransposonDomainsFinder baseFinder;
	
	public DeNovoTransposableElementsFinderShortIR () {
		super();
		setWindowLength(3000);
		setStep(1000);
		setKmerLength(10);
		baseFinder = new HMMTransposonDomainsFinder();
		Set<String> domainsTIR = new HashSet<>();
		domainsTIR.add("HTH");
		domainsTIR.add("MULE");
		domainsTIR.add("MUTATOR");
		domainsTIR.add("HAT");
		baseFinder.loadHMMsFromClasspath(domainsTIR);
	}
	
	protected List<TransposableElementAnnotation> findTransposons(QualifiedSequence seq, int start, int end) {
		CharSequence seqDNA = seq.getCharacters().subSequence(start, end);
		int n = seqDNA.length();
		int kmerLength = getKmerLength();
		CharSequence rc = DNAMaskedSequence.getReverseComplement(seqDNA);
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		Map<Integer,Long> kmersMapForward = KmersExtractor.extractDNAKmerCodesAsMap(seqDNA, kmerLength, 0, seqDNA.length(),true);
		Map<Long, List<Integer>> reverseMapF = KmersExtractor.getReverseMap(kmersMapForward);
		Map<Integer,Long> kmersMapReverse = KmersExtractor.extractDNAKmerCodesAsMap(rc, kmerLength, 0, rc.length(),true);
		Map<Long, List<Integer>> reverseMapR = KmersExtractor.getReverseMap(kmersMapReverse);
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
				lastCandidateStart = startTIR;
				int minInternalLength = baseFinder.getMinLengthProfile()*3;
				int internalLength = tirCandidate.getRightStartRepeat()-tirCandidate.getLeftEndRepeat(); 
				if(internalLength>minInternalLength) {
					//System.out.println("Finding domains for candidate between "+tirCandidate.getFirst()+" and "+tirCandidate.getLast()+" internal length: "+internalLength);
					HMMTransposonDomainsFinder domainsFinder = baseFinder.clone();
					domainsFinder.assignFamily(tirCandidate, candidateTIR);
				}
				if(keepTIRsWithoutDomains || tirCandidate.getInferredFamily()!=null) answer.add(tirCandidate);  
			}
		}
		//if(start > 28000 && start < 33000) System.out.println("Start window: "+start+" Final candidate regions: "+answer.size());
		return answer;
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
	
	public static void main(String[] args) throws Exception {
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		DeNovoTransposableElementsFinderShortIR instance = new DeNovoTransposableElementsFinderShortIR();
		
		instance.setNumThreads(1);
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
