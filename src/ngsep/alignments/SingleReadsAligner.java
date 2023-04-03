package ngsep.alignments;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;

public class SingleReadsAligner {
	private ReferenceGenome genome;
	private UngappedSearchHitsClustersFinder hitClustersFinder;
	private UngappedSearchHitsClusterAligner aligner;
	private double minProportionBestCount = 0.2;
	private double minWeightedCount = 1;
	private int maxAlnsPerRead=1;
	
	
	public SingleReadsAligner(ReferenceGenome genome, UngappedSearchHitsClustersFinder hitClustersFinder, UngappedSearchHitsClusterAligner aligner) {
		super();
		this.genome = genome;
		this.hitClustersFinder = hitClustersFinder;
		this.aligner = aligner;
	}
	public double getMinProportionBestCount() {
		return minProportionBestCount;
	}
	public void setMinProportionBestCount(double minProportionBestCount) {
		this.minProportionBestCount = minProportionBestCount;
	}
	public double getMinWeightedCount() {
		return minWeightedCount;
	}
	public void setMinWeightedCount(double minWeightedCount) {
		this.minWeightedCount = minWeightedCount;
	}
	
	public int getMaxAlnsPerRead() {
		return maxAlnsPerRead;
	}
	public void setMaxAlnsPerRead(int maxAlnsPerRead) {
		this.maxAlnsPerRead = maxAlnsPerRead;
	}
	public List<ReadAlignment> alignRead (QualifiedSequence read) {
		
		List<ReadAlignment> alignments = new ArrayList<>();
		String readSeq = read.getCharacters().toString();
		String qual = read.getQualityScores();
		String reverseQS = null;
		if(qual == null || qual.length()!=readSeq.length()) {
			qual = RawRead.generateFixedQSString('5', readSeq.length());
			reverseQS = qual;
		} else {
			reverseQS = new StringBuilder(qual).reverse().toString();
		}
		String reverseComplement =  DNAMaskedSequence.getReverseComplement(readSeq).toString();
		
		alignments.addAll(alignQuerySequence(readSeq));
		//System.out.println("Read: "+read.getName()+" Forward inexact alignments: "+alignments.size());
		if(reverseComplement!=null) {
			List<ReadAlignment> alnsR = alignQuerySequence(reverseComplement);
			//System.out.println("Read: "+read.getName()+" Reverse inexact alignments: "+alnsR.size());
			for (ReadAlignment aln:alnsR) aln.setNegativeStrand(true);
			alignments.addAll(alnsR);
		}
		
		//System.out.println("Read: "+read.getName()+" total alignments: "+alignments.size());
		for(ReadAlignment aln:alignments) {
			aln.setReadName(read.getName());
			if(!aln.isNegativeStrand()) aln.setQualityScores(qual);
			else aln.setQualityScores(reverseQS);
		}
		return filterAlignments(alignments);
	}
	public List<ReadAlignment> alignQuerySequence(CharSequence query) {
		List<UngappedSearchHitsCluster> clusters = hitClustersFinder.findHitClusters(query);
		return buildAlignments(query, clusters);
	}
	
	public List<ReadAlignment> buildAlignments(CharSequence query, List<UngappedSearchHitsCluster> clusters) {
		double maxCount = summarize(clusters);
		Collections.sort(clusters, (o1,o2)-> ((int)o2.getWeightedCount())-((int)o1.getWeightedCount()));
		double limitCount = Math.min(minWeightedCount, minProportionBestCount*maxCount);
		int limitClusters = Math.min(clusters.size(), Math.max(5, 3*maxAlnsPerRead));
		List<ReadAlignment> answer = new ArrayList<ReadAlignment>();
		//System.out.println("Filtering clusters. Max alns per read: "+maxAlnsPerRead);
		for (int i=0;i<limitClusters;i++) {
			UngappedSearchHitsCluster cluster = clusters.get(i);
			int sequenceIdx = cluster.getSubjectIdx();
			double wc = cluster.getWeightedCount();
			//System.out.println("Qlen: "+query.length()+" next cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" hits "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount());
			if(i>0 && wc<limitCount) break;
			QualifiedSequence refSeq = genome.getSequenceByIndex(sequenceIdx);
			ReadAlignment aln = aligner.buildAlignment(query, refSeq.getCharacters(), cluster);
			//System.out.println("Qlen: "+query.length()+" next cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" hits "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount()+" aln "+aln);
			if(aln!=null) {
				aln.setSequenceName(refSeq.getName());
				answer.add(aln);
			}
		}
		return filterAlignments(answer);
	}
	
	private double summarize(List<UngappedSearchHitsCluster> clusters) {
		double maxCount = 0;
		for (UngappedSearchHitsCluster cluster:clusters) {
			cluster.summarize();
			maxCount = Math.max(maxCount,cluster.getWeightedCount());
			//System.out.println("Summarizing clusters. Next cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+cluster.getNumDifferentKmers()+" count: "+cluster.getWeightedCount()+" maxCount: "+maxCount);
		}
		return maxCount;
	}
	
	private List<ReadAlignment> filterAlignments(List<ReadAlignment> alignments) {
		if (alignments.size()==0) return alignments;
		Collections.sort(alignments, (aln1,aln2) -> aln2.getAlignmentQuality() - aln1.getAlignmentQuality());
		short bestQual = alignments.get(0).getAlignmentQuality();
		//TODO. Investigate alignment score
		int threshold = (int) (0.8*bestQual);
		List<ReadAlignment> filteredAlignments = new ArrayList<>();
		for (int i=0;i<alignments.size();i++) {
			ReadAlignment aln = alignments.get(i);	
			//System.out.println("Aln: "+aln+" qual: "+aln.getAlignmentQuality()+" threshold "+threshold);
			if(aln.getAlignmentQuality()<=threshold) break;
			if(i>0) aln.setSecondary(true);
			filteredAlignments.add(aln);
		}
		int n = filteredAlignments.size();
		
		if(n>1) {
			for(ReadAlignment aln:filteredAlignments) {
				aln.setAlignmentQuality((byte) Math.round(0.3*aln.getAlignmentQuality()/n));
			}
		}
		int limit = Math.min(maxAlnsPerRead, filteredAlignments.size());
		List<ReadAlignment> finalAlignments = new ArrayList<>();
		for (int i=0;i<limit;i++) finalAlignments.add(filteredAlignments.get(i));
		
		//System.out.println("Initial alignments: "+alignments.size()+" final: "+filteredAlignments.size());
		return finalAlignments;
	}

}
