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

public class DeNovoTransposableElementsFinderLTR extends DeNovoTransposableElementsFinderWindowSearch {
	private TransposableElementFamily filterOrder = TransposableElementFamily.LTR_UNKNOWN;
	private HMMTransposonDomainsFinder baseFinder;
	
	public DeNovoTransposableElementsFinderLTR () {
		baseFinder = new HMMTransposonDomainsFinder();
		Set<String> domainsLTR = new HashSet<>();
		domainsLTR.add("GAG");
		domainsLTR.add("INT");
		domainsLTR.add("AP");
		domainsLTR.add("ENV");
		domainsLTR.add("RNASEH");
		domainsLTR.add("RT");
		
		baseFinder.loadHMMsFromClasspath(domainsLTR);
	}
	
	protected List<TransposableElementAnnotation> findTransposons(QualifiedSequence seq, int start, int end) {
		Logger log = getLog();
		int kmerLength = getKmerLength();
		int windowLength = getWindowLength();
		int step = getStep();
		log.info("Processing sequence. "+seq.getName()+" from "+start+" to "+end);
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		
		String segmentDNA = seq.getCharacters().subSequence(start, end).toString();
		int n = segmentDNA.length();
		Map<Integer,Long> kmersMapForward = KmersExtractor.extractDNAKmerCodesAsMap(segmentDNA, kmerLength, 0, n,true);
		Map<Long, List<Integer>> reverseMapF = getReverseMap(kmersMapForward);
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
		PairwiseAlignerSimpleGap pwa = new PairwiseAlignerSimpleGap(2000);
		pwa.setForceEnd1(false);
		pwa.setForceEnd2(false);
		HMMTransposonDomainsFinder domainsFinder = baseFinder.clone();
		Collections.sort(clusters,(c1,c2)->c2.getCountKmerHitsCluster()-c1.getCountKmerHitsCluster());
		for(UngappedSearchHitsCluster cluster:clusters) {
			TransposableElementAnnotation ltrAnn = inferTEFromEndAlignment(seq, start, segmentDNA, cluster, pwa);
			if(ltrAnn==null) continue;
			//Assign family
			assignFamily(ltrAnn, seq, domainsFinder);
			if(passFilters(ltrAnn)) answer.add(ltrAnn);
		}
		//TODO: Discard possible hits due to tandem repeats
		return answer;
	}
	

	private TransposableElementAnnotation inferTEFromEndAlignment(QualifiedSequence seq, int start, String segmentDNA, UngappedSearchHitsCluster cluster, PairwiseAlignerSimpleGap pwa) {
		System.err.println("Checking hit of sequence. Start: "+start);
		//Find start alignment before hit
		int evidenceStart1 = start + cluster.getQueryEvidenceStart();
		int evidenceEnd1 = start + cluster.getQueryEvidenceEnd();
		int evidenceStart2 = start + cluster.getSubjectEvidenceStart();
		int evidenceEnd2 = start + cluster.getSubjectEvidenceEnd();
		
		
		DNAMaskedSequence dna = (DNAMaskedSequence) seq.getCharacters();
		DNAMaskedSequence leftSegment = ((DNAMaskedSequence)dna.subSequence(Math.max(0, evidenceStart1-300), evidenceStart1)).getReverseComplement();
		DNAMaskedSequence rightSegment = ((DNAMaskedSequence)dna.subSequence(Math.max(0, evidenceStart2-300), evidenceStart2)).getReverseComplement();
		
		
		PairwiseAlignment alnBeforeHit = pwa.calculateAlignment(leftSegment, rightSegment);
		int start2 = Math.max(0, evidenceStart1 - alnBeforeHit.getEnd1());
		int internalRight = Math.max(0, evidenceStart2-alnBeforeHit.getEnd2());
		
		//Find end alignment after hit
		leftSegment = (DNAMaskedSequence)dna.subSequence(evidenceStart1, Math.min(seq.getLength(), evidenceStart1+2000));
		rightSegment = (DNAMaskedSequence)dna.subSequence(evidenceStart2, Math.min(seq.getLength(), evidenceStart2+2000));
		PairwiseAlignment alnAfterHit = pwa.calculateAlignment(leftSegment, rightSegment);
		int internalLeft = evidenceStart1+alnAfterHit.getEnd1();
		int end2 = Math.min(seq.getLength(), evidenceStart2+alnAfterHit.getEnd2());
		
		System.err.println("Checked borders. Start: "+start2+" internal left "+internalLeft+" internal right: "+internalRight+" end "+end2);
		if(internalLeft<evidenceEnd1) return null;
		if(internalRight<=internalLeft) return null;
		if(end2 < evidenceEnd2) return null;
		//Create annotation
		TransposableElementAnnotation ann = new TransposableElementAnnotation(seq.getName(), start2+1, end2);
		ann.setRepeatLimits(internalLeft+1, internalRight+1, TransposableElementFamily.REPEAT_ORIENTATION_FF);
		return ann;
	}
	private void assignFamily(TransposableElementAnnotation ann, QualifiedSequence seq, HMMTransposonDomainsFinder domainsFinder) {
		DNAMaskedSequence dnaTE = (DNAMaskedSequence) seq.getCharacters().subSequence(ann.getLeftEndRepeat()-1, ann.getRightStartRepeat());
		domainsFinder.assignFamily(ann,dnaTE);
		if(ann.getInferredFamily()==null) {
			domainsFinder.assignFamily(ann,dnaTE.getReverseComplement());
			if(ann.getInferredFamily()!=null) ann.setNegativeStrand(true);
		}
	}
	
	
	private boolean passFilters(TransposableElementAnnotation ann) {
		if(ann.getInferredFamily()==null) return false;
		if(filterOrder!=null) {
			TransposableElementFamily family = ann.getInferredFamily();
			if(family==null || !family.getOrder().equals(filterOrder.getOrder())) return false;
		}
		return true;
	}

	public static void main(String[] args) throws Exception {
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		DeNovoTransposableElementsFinderWindowSearch instance = new DeNovoTransposableElementsFinderLTR();
		instance.setNumThreads(8);
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
