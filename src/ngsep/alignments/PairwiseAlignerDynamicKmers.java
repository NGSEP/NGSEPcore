package ngsep.alignments;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.sequences.KmersExtractor;
import ngsep.sequences.UngappedSearchHit;

public class PairwiseAlignerDynamicKmers implements PairwiseAligner {
	private static int debugLength = -1;
	private PairwiseAlignerStaticBanded alignerBanded = new PairwiseAlignerStaticBanded();
	@Override
	public String[] calculateAlignment(CharSequence sequence1, CharSequence sequence2) {
		int n1 = sequence1.length();
		int n2 = sequence2.length();
		if(n2 == debugLength) System.out.println("Aligning sequences\n"+sequence1+"\n"+sequence2);
		if(n1<n2) {
			String[] alnRev = calculateAlignment(sequence2, sequence1);
			if(alnRev == null) return null;
			String [] answer = {alnRev[1].toString(),alnRev[0].toString()};
			return answer;
		}
		if(n1<50 && n2==0) return (new PairwiseAlignerNaive(false)).calculateAlignment(sequence1, sequence2);
		if(n1<100 && n2<20) {
			PairwiseAlignerSimpleGap aligner = new PairwiseAlignerSimpleGap();
			aligner.setForceStart1(false);
			aligner.setForceEnd1(false);
			return aligner.calculateAlignment(sequence1, sequence2);
		}
		if(0.3*n1>n2) {
			//System.err.println("WARN: Unbalanced lengths for global alignment: "+n1+" "+n2);
			return null;
		}
		int kmerLength = Math.max(11, Math.min(n1, n2)/50);
		kmerLength = Math.min(31, kmerLength);
		Map<Integer,Long> codesSubject = KmersExtractor.extractDNAKmerCodesAsMap(sequence1, kmerLength, 0, sequence1.length());
		Map<Integer,Long> codesQuery = KmersExtractor.extractDNAKmerCodesAsMap(sequence2, kmerLength, 0, sequence2.length());
		UngappedSearchHitsCluster bestCluster = MinimizersUngappedSearchHitsClustersFinder.findBestKmersCluster(n1, codesSubject, n2, codesQuery, kmerLength);
		if(bestCluster==null) {
			//System.err.println("WARN: Null cluster for alignment of sequences with lengths: "+n1+" "+n2+" kmer length: "+kmerLength+" kmers: "+codesSubject.size()+" "+codesQuery.size()+" sequences\n"+sequence1+"\n"+sequence2);
			return null;
		}
		List<UngappedSearchHit> kmerHits = bestCluster.getHitsByQueryIdx();
		int subjectNext = 0;
		if(n2 == debugLength) System.out.println("S1 length: "+n1+". S2 length: "+n2+" kmer hits: "+kmerHits.size()+" subject next: "+subjectNext+" cluster predicted start: "+bestCluster.getSubjectPredictedStart()+" kmer length: "+kmerLength);
		int queryNext = 0;
		int alnStart = -1;
		StringBuilder aln1 = new StringBuilder();
		StringBuilder aln2 = new StringBuilder();
		int nextMatchLength = 0;
		for(UngappedSearchHit kmerHit:kmerHits) {
			if(n2 == debugLength) System.out.println("Processing Kmer hit at pos: "+kmerHit.getQueryStart()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getSubjectStart());
			if(alnStart==-1) {
				//Inconsistent kmer hit
				if (kmerHit.getSubjectStart()<bestCluster.getSubjectPredictedStart()) continue;
				alnStart = kmerHit.getSubjectStart();
				String seq1Fragment = sequence1.subSequence(0,alnStart).toString();
				String seq2Fragment = sequence2.subSequence(0,kmerHit.getQueryStart()).toString();
				if(seq1Fragment.length()>0 || seq2Fragment.length()>0) {
					if(n2 == debugLength) System.out.println("Aligning "+seq1Fragment+" with "+seq2Fragment+" lengths: "+seq1Fragment.length()+" "+seq2Fragment.length());
					String [] alignedFragments = calculateAlignment(seq1Fragment, seq2Fragment);
					if(n2 == debugLength && alignedFragments==null) System.out.println("Null start alignment between "+seq1Fragment+" and "+seq2Fragment);
					if(alignedFragments==null) return null;
					aln1.append(alignedFragments[0]);
					aln2.append(alignedFragments[1]);
				}
				nextMatchLength+=kmerLength;
				subjectNext = kmerHit.getSubjectStart()+kmerLength;
				queryNext = kmerHit.getQueryStart()+kmerLength;
			} else if(queryNext < kmerHit.getQueryStart() && subjectNext<kmerHit.getSubjectStart()) {
				//Kmer does not overlap with already aligned segments
				int subjectNextLength = kmerHit.getSubjectStart()-subjectNext;
				int queryNextLength = kmerHit.getQueryStart()-queryNext;
				if(subjectNextLength==queryNextLength) {
					nextMatchLength+=subjectNextLength;
				} else {
					int minLength = Math.min(subjectNextLength, queryNextLength);
					int maxLength = Math.max(subjectNextLength, queryNextLength);
					int diffLength = maxLength-minLength;
					if(diffLength>3 && 0.95*maxLength>minLength) {
						//Possible invalid kmer hit. Delay alignment
						if(n2 == debugLength) System.out.println("Possible invalid kmer hit. Kmer hit at pos: "+kmerHit.getQueryStart()+" subject hit start: "+kmerHit.getSubjectStart()+" Subject length "+subjectNextLength+" query length "+queryNextLength);
						continue;
					}
					if(nextMatchLength>0) {
						aln1.append(sequence1,subjectNext-nextMatchLength, subjectNext);
						aln2.append(sequence2,queryNext-nextMatchLength, queryNext);
						nextMatchLength = 0;
					}
					
					String seq1Fragment = sequence1.subSequence(subjectNext,kmerHit.getSubjectStart()).toString();
					String seq2Fragment = sequence2.subSequence(queryNext,kmerHit.getQueryStart()).toString();
					if(n2 == debugLength)  System.out.println("Aligning segment of length "+subjectNextLength+" of subject with total length: "+n1+" to segment with length "+queryNextLength+" of query with total length: "+n2);
					String [] alignedFragments = calculateAlignment(seq1Fragment,seq2Fragment);
					if(n2 == debugLength && alignedFragments==null) System.out.println("Null middle alignment between "+seq1Fragment+" and "+seq2Fragment);  
					if(alignedFragments==null && diffLength<=10) {
						//Try with the static band if the length difference is small
						alignedFragments = alignerBanded.calculateAlignment(seq1Fragment, seq2Fragment);
					}
					if(alignedFragments==null) return null;
					aln1.append(alignedFragments[0]);
					aln2.append(alignedFragments[1]);
				}
				nextMatchLength+=kmerLength;
				subjectNext = kmerHit.getSubjectStart()+kmerLength;
				queryNext = kmerHit.getQueryStart()+kmerLength;
			} else {
				int kmerSubjectNext = kmerHit.getSubjectStart()+kmerLength;
				int diffSubject = kmerSubjectNext - subjectNext;
				int kmerQueryNext = kmerHit.getQueryStart()+kmerLength;
				int diffQuery = kmerQueryNext - queryNext;
				if (diffSubject>0 && diffSubject==diffQuery) {
					//Kmer consistent with current alignment. Augment match with difference
					nextMatchLength+=diffSubject;
					subjectNext = kmerSubjectNext;
					queryNext = kmerQueryNext;
				}
			}
			
			//System.out.println("Processed Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext);
		}
		if(nextMatchLength>0) {
			aln1.append(sequence1,subjectNext-nextMatchLength, subjectNext);
			aln2.append(sequence2,queryNext-nextMatchLength, queryNext);
			nextMatchLength = 0;
		}
		if(subjectNext<n1 || queryNext<n2) {
			int subjectNextLength = n1-subjectNext;
			int queryNextLength = n2-queryNext;
			if(subjectNextLength>0 || queryNextLength>0) {
				String seq1Fragment = (subjectNextLength>0)?sequence1.subSequence(subjectNext,n1).toString():"";
				String seq2Fragment = (queryNextLength>0)?sequence2.subSequence(queryNext,n2).toString():"";
				if(n2 == debugLength) System.out.println("Aligning segment of length "+subjectNextLength+" of subject with total length: "+n1+" to segment with length "+queryNextLength+" of query with total length: "+n2);
				String [] alignedFragments = calculateAlignment(seq1Fragment,seq2Fragment);
				if(n2 == debugLength && alignedFragments==null) System.out.println("Null end alignment between "+seq1Fragment+" and "+seq2Fragment);
				if(alignedFragments==null) return null;
				aln1.append(alignedFragments[0]);
				aln2.append(alignedFragments[1]);
			}
		}
		String [] answer = {aln1.toString(),aln2.toString()};
		//System.out.println("Segment alignment\n"+aln1+"\n"+aln2);
		return answer;
	}

}
