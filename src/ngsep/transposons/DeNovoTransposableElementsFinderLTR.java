package ngsep.transposons;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.alignments.PairwiseAlignment;
import ngsep.alignments.UngappedSearchHitsCluster;
import ngsep.alignments.UngappedSearchHitsClusterBuilder;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.io.FastaSequencesHandler;

public class DeNovoTransposableElementsFinderLTR extends DeNovoTransposableElementsFinderWindowSearch {
	private TransposableElementFamily filterOrder = TransposableElementFamily.LTR_UNKNOWN;
	//private TransposableElementFamily filterOrder = null;
	private String filterDomainCode = "GAG";
	private HMMTransposonDomainsFinder baseFinder;
	
	private int debugPos = -1;
	
	public DeNovoTransposableElementsFinderLTR () {
		baseFinder = new HMMTransposonDomainsFinder();
		Set<String> domainsLTR = new HashSet<>();
		domainsLTR.add("GAG");
		domainsLTR.add("INT");
		domainsLTR.add("AP");
		domainsLTR.add("ENV");
		domainsLTR.add("RNASEH");
		domainsLTR.add("RT");
		domainsLTR.add("END");
		domainsLTR.add("RH");
		baseFinder.loadHMMsFromClasspath(domainsLTR);
	}
	
	protected List<TransposableElementAnnotation> findTransposons(QualifiedSequence seq, int start, int end, List<TransposableElementAnnotation> overlappingEvents) {
		Logger log = getLog();
		int kmerLength = getKmerLength();
		int windowLength = getWindowLength();
		int step = getStep();
		log.info("Processing sequence. "+seq.getName()+" from "+start+" to "+end);
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		
		String segmentDNA = seq.getCharacters().subSequence(start, end).toString();
		if(overlappingEvents.size()>0) segmentDNA = maskSegment(segmentDNA,start,overlappingEvents);
		int n = segmentDNA.length();
		Map<Integer,Long> kmersMapForward = KmersExtractor.extractDNAKmerCodesAsMap(segmentDNA, kmerLength, 0, n,true);
		Map<Long, List<Integer>> reverseMapF = KmersExtractor.getReverseMap(kmersMapForward);
		List<UngappedSearchHit> hits = new ArrayList<UngappedSearchHit>();
		int totalMultihitKmers = 0;
		int totalHitsMultihitKmers = 0;
		for(Map.Entry<Integer,Long> entry:kmersMapForward.entrySet()) {
			int i1 = entry.getKey();
			List<Integer> hitPos = reverseMapF.get(entry.getValue());
			if(hitPos.size()<2) continue;
			totalHitsMultihitKmers+=hitPos.size();
			totalMultihitKmers++;
			for(int i2:hitPos) {
				if(i2<i1+500) continue;
				UngappedSearchHit hit = new UngappedSearchHit(0, i2);
				hit.setQueryStart(i1);
				hit.setHitLength((short) kmerLength);
				hits.add(hit);
			}
		}
		//Filter for tandem repeat regions
		if(totalHitsMultihitKmers>5*totalMultihitKmers) return answer;
		UngappedSearchHitsClusterBuilder builder = new UngappedSearchHitsClusterBuilder();
		List<UngappedSearchHitsCluster> clusters = builder.clusterRegionKmerAlns(step,0, windowLength, hits);
		if(clusters.size()==0) return answer;
		log.info("Processing sequence. "+seq.getName()+" from "+start+" to "+end+" total hits: "+hits.size()+" total clusters: "+clusters.size());
		PairwiseAlignerSimpleGap pwa = new PairwiseAlignerSimpleGap(3000);
		pwa.setForceEnd1(false);
		pwa.setForceEnd2(false);
		
		
		HMMTransposonDomainsFinder domainsFinder = null;
		Collections.sort(clusters,(c1,c2)->c2.getCountKmerHitsCluster()-c1.getCountKmerHitsCluster());
		for(UngappedSearchHitsCluster cluster:clusters) {
			TransposableElementAnnotation ltrAnn = inferTEFromEndAlignment(seq, start, segmentDNA, cluster, pwa);
			if(ltrAnn==null) continue;
			//Assign family
			if(domainsFinder==null) domainsFinder = baseFinder.clone();
			assignFamily(ltrAnn, seq, domainsFinder);
			findTargetSideDuplication(ltrAnn, seq);
			log.info("Assigned family for annotation at "+ltrAnn.getFirst()+" "+ltrAnn.getLast()+" inferred family: "+ltrAnn.getInferredFamily());
			answer.add(ltrAnn);
		}
		return answer;
	}
	

	

	private String maskSegment(String segmentDNA, int start, List<TransposableElementAnnotation> overlappingEvents) {
		StringBuilder segmentMod = new StringBuilder(segmentDNA);
		for(TransposableElementAnnotation ann:overlappingEvents) {
			int startMask = Math.max(0, ann.getFirst()-1-start);
			int endMask = Math.min(segmentDNA.length(), ann.getLast()-start);
			for(int i=startMask;i<endMask;i++) segmentMod.setCharAt(i, 'N');
		}
		return segmentMod.toString();
	}

	private void findTargetSideDuplication(TransposableElementAnnotation ltrAnn, QualifiedSequence seq) {
		CharSequence dna = seq.getCharacters();
		int refStartLeft = Math.max(1,ltrAnn.getFirst()-20);
		int refStartRight = ltrAnn.getLast()-5;
		int refEndRight = Math.min(dna.length(),ltrAnn.getLast()+20);
		String leftSegment = dna.subSequence(refStartLeft-1, ltrAnn.getFirst()+5).toString();
		String rightSegment = dna.subSequence(refStartRight-1, refEndRight).toString();
		int tsdMinLength = 6;
		int k =4;
		int tsdStartLeft = -1;
		int tsdStartRight = -1;
		String tsd= null;
		for(int i=0;i<=rightSegment.length()-k;i++) {
			String kmer = rightSegment.substring(i,i+k);
			int idxL = leftSegment.indexOf(kmer);
			if(idxL<0) continue;
			int l = k;
			int j1 = idxL+k;
			int j2 = i+k;
			while(j1<leftSegment.length() && j2<rightSegment.length()) {
				if(leftSegment.charAt(j1)!=rightSegment.charAt(j2)) break;
				j1++;
				j2++;
				l++;
			}
			if(l>=tsdMinLength) {
				tsdStartLeft = idxL;
				tsdStartRight = i;
				tsd = rightSegment.substring(i,i+l);
				break;
			}
		}
		if(tsd!=null) {
			ltrAnn.setFirst(refStartLeft+tsdStartLeft+tsd.length());
			ltrAnn.setLast(refStartRight+tsdStartRight-1);
			ltrAnn.setTsd(tsd);
		}
	}

	private TransposableElementAnnotation inferTEFromEndAlignment(QualifiedSequence seq, int start, String segmentDNA, UngappedSearchHitsCluster cluster, PairwiseAlignerSimpleGap pwa) {
		if (start==debugPos) System.err.println("Checking hit of sequence. Start: "+start);
		//Find start alignment before hit
		int evidenceStart1 = start + cluster.getQueryEvidenceStart();
		int evidenceEnd1 = start + cluster.getQueryEvidenceEnd();
		int evidenceStart2 = start + cluster.getSubjectEvidenceStart();
		int evidenceEnd2 = start + cluster.getSubjectEvidenceEnd();
		int evidenceLength1 = evidenceEnd1-evidenceStart1;
		int evidenceLength2 = evidenceEnd2-evidenceStart2;
		
		int ltrEstimatedLength = Math.min(evidenceLength1, evidenceLength2);
		ltrEstimatedLength = Math.min(ltrEstimatedLength, 5000);
		
		if (start==debugPos) System.err.println("Checking hit of sequence. Start: "+start+" LTR length: "+ltrEstimatedLength+ " ev1: "+evidenceStart1+" "+ evidenceEnd1+" ev2: "+evidenceStart2+" "+evidenceEnd2);
		if(evidenceStart2 - evidenceEnd1 <500) return null;
		DNAMaskedSequence dna = (DNAMaskedSequence) seq.getCharacters();
		DNAMaskedSequence leftSegment = ((DNAMaskedSequence)dna.subSequence(Math.max(0, evidenceStart1-300), evidenceStart1)).getReverseComplement();
		DNAMaskedSequence rightSegment = ((DNAMaskedSequence)dna.subSequence(Math.max(0, evidenceStart2-300), evidenceStart2)).getReverseComplement();
		
		
		PairwiseAlignment alnBeforeHit = pwa.calculateAlignment(leftSegment, rightSegment);
		//TODO: Alignment quality
		
		int start2 = Math.max(0, evidenceStart1 - alnBeforeHit.getEnd1());
		int internalRight = Math.max(0, evidenceStart2-alnBeforeHit.getEnd2());
		
		if (start==debugPos) System.err.println("Checking hit of sequence. Start: "+start+" endsAlnBefore: "+alnBeforeHit.getEnd1()+" "+alnBeforeHit.getEnd2()+" mismatches: "+alnBeforeHit.getMismatches()+" Left aln limits: "+start2+ " "+internalRight);
		
		//Find end alignment after hit
		leftSegment = (DNAMaskedSequence)dna.subSequence(evidenceStart1, Math.min(seq.getLength(), evidenceStart1+ltrEstimatedLength));
		rightSegment = (DNAMaskedSequence)dna.subSequence(evidenceStart2, Math.min(seq.getLength(), evidenceStart2+ltrEstimatedLength));
		PairwiseAlignment alnAfterHit = pwa.calculateAlignment(leftSegment, rightSegment);
		int internalLeft = evidenceStart1+alnAfterHit.getEnd1();
		int end2 = Math.min(seq.getLength(), evidenceStart2+alnAfterHit.getEnd2());
		if (start==debugPos) System.err.println("Checking hit of sequence. Start: "+start+" endsAlnAfter: "+alnAfterHit.getEnd1()+" "+alnAfterHit.getEnd2()+" mismatches: "+alnAfterHit.getMismatches()+" Right aln limits: "+internalLeft+ " "+end2);
		if(alnAfterHit.getMismatches()>0.25*ltrEstimatedLength) return null;
		if(internalLeft-start2<50) return null;
		
		if (start==debugPos) System.err.println("Checked borders. Start: "+start2+" internal left "+internalLeft+" internal right: "+internalRight+" end "+end2+" Mismatches: "+alnBeforeHit.getMismatches()+" "+alnAfterHit.getMismatches());
		if(internalRight-internalLeft<500) return null;
		//Create annotation
		TransposableElementAnnotation ann = new TransposableElementAnnotation(seq.getName(), start2+1, end2);
		ann.setRepeatLimits(internalLeft+1, internalRight+1, TransposableElementFamily.REPEAT_ORIENTATION_FF);
		if (start==debugPos) System.err.println("Created annotation at "+ann.getFirst()+" "+ann.getLast());
		return ann;
	}
	private void assignFamily(TransposableElementAnnotation ann, QualifiedSequence seq, HMMTransposonDomainsFinder domainsFinder) {
		DNAMaskedSequence dnaTE = (DNAMaskedSequence) seq.getCharacters().subSequence(ann.getLeftEndRepeat()-1, ann.getRightStartRepeat());
		getLog().info("Verifying family for TE at "+ann.getFirst()+" "+ann.getLast()+" internal: "+ann.getLeftEndRepeat()+" "+ann.getRightStartRepeat()+" seq len: "+dnaTE.length());
		domainsFinder.assignFamily(ann,dnaTE);
	}
	
	
	protected boolean passFilters(TransposableElementAnnotation ann) {
		if(ann.getInferredFamily()==null) return false;
		if(filterOrder!=null) {
			TransposableElementFamily family = ann.getInferredFamily();
			if(family==null || !family.getOrder().equals(filterOrder.getOrder())) return false;
		}
		if(filterDomainCode!=null) {
			boolean found = false;
			for(TransposonDomainAlignment aln: ann.getDomainAlignments()) {
				if(filterDomainCode.equals(aln.getDomainCode())) {
					found = true;
					break;
				}
			}
			if(!found) return false;
		}
		return true;
	}

	public static void main(String[] args) throws Exception {
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		String outPrefix = args[1];
		DeNovoTransposableElementsFinderWindowSearch instance = new DeNovoTransposableElementsFinderLTR();
		instance.setNumThreads(8);
		List<TransposableElementAnnotation> anns = instance.findTransposons(genome);
		try (PrintStream out=new PrintStream(outPrefix+"_regions.txt")) {
			for(TransposableElementAnnotation ann:anns) {
				out.print(""+ann.getSequenceName()+"\t"+ann.getFirst()+"\t"+ann.getLast()+"\t"+ann.isPositiveStrand()+"\t"+ann.getInferredFamily()+"\t"+ann.getLeftEndRepeat()+"\t"+ann.getRightStartRepeat()+"\t"+ann.getOrientation()+"\t"+ann.getTsd());
				if(ann.getDomainAlignments()!=null && ann.getDomainAlignments().size()>0) {
					out.print("\t");
					for(TransposonDomainAlignment daln:ann.getDomainAlignments()) {
						out.print(daln.getDomainCode()+":"+daln.getStart()+":"+daln.getEvalue()+";");
					}
				}
				out.println();
			}
		}
		FastaSequencesHandler handler = new FastaSequencesHandler();
		try (PrintStream out=new PrintStream(outPrefix+"_sequences.fa")) {
			for(TransposableElementAnnotation ann:anns) {
				CharSequence seq = genome.getReference(ann);
				if(ann.isNegativeStrand()) seq = DNAMaskedSequence.getReverseComplement(seq);
				QualifiedSequence qseq = new QualifiedSequence(""+ann.getSequenceName()+"_"+ann.getFirst()+"_"+ann.getLast()+"_"+ann.isPositiveStrand()+"#"+ann.getInferredFamily(),seq);
				handler.saveSequence(qseq, out, 100);
			}
		}
	}

}
