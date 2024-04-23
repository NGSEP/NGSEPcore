/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.alignments;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.LimitedSequence;
import ngsep.sequences.UngappedSearchHit;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class LongReadsUngappedSearchHitsClusterAligner implements UngappedSearchHitsClusterAligner {
	private int maxLengthFullPairwiseAlignment = 10000;
	private int maxLengthEndsPairwiseAlignment = 10000;
	private PairwiseAligner alignerCenter;
	private PairwiseAligner alignerStart;
	private PairwiseAligner alignerEnd;
	private boolean trySimpleGap = false;
	
	private int subjectIdxDebug = -1;
	private int queryLengthDebug = -1;
	
	public LongReadsUngappedSearchHitsClusterAligner () {
		this(ALIGNMENT_ALGORITHM_AFFINE_GAP);
	}
	public LongReadsUngappedSearchHitsClusterAligner (int alignmentAlgorithm) {
		if(alignmentAlgorithm == ALIGNMENT_ALGORITHM_AFFINE_GAP) {
			maxLengthEndsPairwiseAlignment = 500;
			maxLengthFullPairwiseAlignment = 5000;
			alignerCenter = new PairwiseAlignerAffineGap(maxLengthFullPairwiseAlignment+1);
			alignerStart = new PairwiseAlignerAffineGap(maxLengthEndsPairwiseAlignment);
			alignerEnd = new PairwiseAlignerAffineGap(maxLengthEndsPairwiseAlignment);
			((PairwiseAlignerAffineGap)alignerStart).setForceStart2(false);
			((PairwiseAlignerAffineGap)alignerEnd).setForceEnd2(false);
		}
		if(alignmentAlgorithm == ALIGNMENT_ALGORITHM_DYNAMIC_KMERS) {
			alignerCenter = new PairwiseAlignerDynamicKmers();
			alignerStart = new PairwiseAlignerDynamicKmers();
			alignerEnd = new PairwiseAlignerDynamicKmers();
		}
		if(alignmentAlgorithm == ALIGNMENT_ALGORITHM_STATIC_BAND) {
			alignerCenter = new PairwiseAlignerStaticBanded();
			alignerStart = new PairwiseAlignerStaticBanded();
			alignerEnd = new PairwiseAlignerStaticBanded();
		}
		if(alignmentAlgorithm == ALIGNMENT_ALGORITHM_SIMPLE_GAP) {
			maxLengthEndsPairwiseAlignment = 500;
			maxLengthFullPairwiseAlignment = 5000;
			alignerCenter = new PairwiseAlignerSimpleGap(maxLengthFullPairwiseAlignment+1);
			alignerStart = new PairwiseAlignerSimpleGap(maxLengthEndsPairwiseAlignment);
			alignerEnd = new PairwiseAlignerSimpleGap(maxLengthEndsPairwiseAlignment);
		}
		if(alignmentAlgorithm == ALIGNMENT_ALGORITHM_NAIVE) {
			alignerCenter = new PairwiseAlignerNaive(false);
			alignerStart = new PairwiseAlignerNaive(true);
			alignerEnd = new PairwiseAlignerNaive(false);
		}
	}
	
	public ReadAlignment buildAlignment(CharSequence query, CharSequence subject, UngappedSearchHitsCluster kmerHitsCluster) { 
		int subjectIdx = kmerHitsCluster.getSubjectIdx();
		List<UngappedSearchHit> kmerHits = kmerHitsCluster.getHitsByQueryIdx();
		String queryS = query.toString();
		int queryLength= query.length();
		//System.out.println("Query length: "+queryLength);
		HammingSequenceDistanceMeasure hamming = new HammingSequenceDistanceMeasure();
		int subjectNext = -1;
		short numMismatches = 0;
		int coverageSharedKmers = 0;
		double weightedCoverageSharedKmers = 0;
		if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligning long read. Subject length: "+subject.length()+". Query length: "+query.length()+" kmer hits: "+kmerHits.size()+" cluster predicted limits: "+kmerHitsCluster.getSubjectPredictedStart()+ " "+kmerHitsCluster.getSubjectPredictedEnd());
		int queryNext = 0;
		int alnStart = -1;
		int queryStart = -1;
		List<Integer> alignmentEncoding = new ArrayList<Integer>();
		int nextMatchLength = 0;
		for(UngappedSearchHit kmerHit:kmerHits) {
			int hitLength = kmerHit.getHitLength();
			int hitPredictedStart = kmerHit.getSubjectStart()-kmerHit.getQueryStart();
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Processing Kmer hit at pos: "+kmerHit.getQueryStart()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getSubjectStart()+" Predicted start subject by kmer: "+hitPredictedStart);
			if(alnStart==-1) {
				
				int absDiffPredictedStarts = Math.abs(kmerHitsCluster.getSubjectPredictedStart()-hitPredictedStart);
				if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Candidate hit start. QueryPos: "+kmerHit.getQueryStart()+" subject start: "+kmerHit.getSubjectStart()+" Predicted start kmer: "+hitPredictedStart);
				//Inconsistent kmer hit
				if (kmerHit.getSubjectStart()<kmerHitsCluster.getSubjectPredictedStart() || absDiffPredictedStarts>30) continue;
				alnStart = kmerHit.getSubjectStart();
				queryStart = kmerHit.getQueryStart();
				boolean startAligned = queryStart<=0;
				if(!startAligned && queryStart<kmerHit.getSubjectStart()) {
					String queryStr = queryS.substring(0,queryStart);
					//TODO: Take into account cluster predicted subject start
					int possibleAlnStart = Math.max(0, kmerHit.getSubjectStart()-queryStart-10);
					String subjectStr = subject.subSequence(possibleAlnStart,kmerHit.getSubjectStart()).toString();
					if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Hit start. Query segment: "+queryStr+" subject segment: "+subjectStr);
					String [] alignedFragments = alignFragments(queryStr, subjectStr, subjectIdx, queryLength, alignerStart, 0);
					if(alignedFragments!=null) {
						
						String [] trimmedAln = trimStartDeletions(alignedFragments);
						int trimmedBP = alignedFragments[0].length()-trimmedAln[0].length();
						if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligned fragments: \n"+alignedFragments[0]+"\n"+alignedFragments[1]+"\ntrimmed bp: "+trimmedBP);
						alignedFragments = trimmedAln;
						alignmentEncoding.addAll(ReadAlignment.encodePairwiseAlignment(alignedFragments));
						double mismatchesSegment = hamming.calculateDistance(alignedFragments[0], alignedFragments[1]); 
						numMismatches+=mismatchesSegment;
						if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligned and trimmed fragments: \n"+alignedFragments[0]+"\n"+alignedFragments[1]+"\nSegment mismatches: "+mismatchesSegment+" total mismatches: "+numMismatches);
						startAligned = true;
						queryStart=0;
						alnStart = possibleAlnStart+trimmedBP;
					}
					
				} 
				if(!startAligned) alignmentEncoding.add(ReadAlignment.getAlnValue(queryStart, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
				nextMatchLength+=hitLength;
				coverageSharedKmers+=hitLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)hitLength*weight);
				subjectNext = kmerHit.getSubjectStart()+hitLength;
				queryNext = kmerHit.getQueryStart()+hitLength;
			} else if(kmerHit.getQueryStart() > queryNext && subjectNext<kmerHit.getSubjectStart()) {
				//Kmer does not overlap with already aligned segments
				int subjectNextLength = kmerHit.getSubjectStart()-subjectNext;
				int queryNextLength = kmerHit.getQueryStart()-queryNext;
				boolean goodMatch = subjectNextLength==queryNextLength;
				double hammingDistance = goodMatch?hamming.calculateDistance(subject.subSequence(subjectNext, kmerHit.getSubjectStart()), queryS.substring(queryNext, kmerHit.getQueryStart())):0;
				goodMatch = goodMatch && hammingDistance<5;
				if(goodMatch) {
					nextMatchLength+=subjectNextLength;
					numMismatches+=hammingDistance;
					if (subjectIdx == subjectIdxDebug  && queryLength==queryLengthDebug) System.out.println("Aligning equal length strings. Kmer hit at pos: "+kmerHit.getQueryStart()+" subject hit start: "+kmerHit.getSubjectStart()+" Subject length "+subjectNextLength+" query length "+queryNextLength+" total mismatches: "+numMismatches);
				}
				else {
					int minLength = Math.min(subjectNextLength, queryNextLength);
					int maxLength = Math.max(subjectNextLength, queryNextLength);
					if(maxLength>minLength+3 && 0.95*maxLength>minLength) {
						//Possible invalid kmer hit. Delay alignment
						if (subjectIdx == subjectIdxDebug  && queryLength==queryLengthDebug) System.out.println("Possible invalid kmer hit. Kmer hit at pos: "+kmerHit.getQueryStart()+" subject hit start: "+kmerHit.getSubjectStart()+" Subject length "+subjectNextLength+" query length "+queryNextLength);
						continue;
					}
					if(nextMatchLength>0) {
						if (subjectIdx == subjectIdxDebug  && queryLength==queryLengthDebug) System.out.println("Found internal segment for possible alignment. Kmer hit at pos: "+kmerHit.getQueryStart()+" subject hit start: "+kmerHit.getSubjectStart()+" Subject length "+subjectNextLength+" query length "+queryNextLength+" current match length: "+nextMatchLength);
						alignmentEncoding.add(ReadAlignment.getAlnValue(nextMatchLength, ReadAlignment.ALIGNMENT_MATCH));
						nextMatchLength = 0;
					}
					String subjectStr = subject.subSequence(subjectNext,kmerHit.getSubjectStart()).toString();
					String queryStr = queryS.substring(queryNext,kmerHit.getQueryStart()).toString();
					String [] alignedFragments = alignFragments(queryStr, subjectStr, subjectIdx, queryLength, alignerCenter, 1);
					if(alignedFragments==null) {
						if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Null alignment. Trying default alignment");
						if(maxLength>0.2*query.length()) return null;
						alignmentEncoding.add(ReadAlignment.getAlnValue(minLength, ReadAlignment.ALIGNMENT_MISMATCH));
						if(subjectNextLength>queryNextLength) alignmentEncoding.add(ReadAlignment.getAlnValue(subjectNextLength-queryNextLength, ReadAlignment.ALIGNMENT_DELETION));
						else if (subjectNextLength<queryNextLength) alignmentEncoding.add(ReadAlignment.getAlnValue(queryNextLength-subjectNextLength, ReadAlignment.ALIGNMENT_INSERTION));
						if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug)  System.out.println("Default alignmet for query coords "+queryNext+" "+kmerHit.getQueryStart()+" length: "+queryStr.length()+" subject coords: "+subjectNext+" " +kmerHit.getSubjectStart());
						numMismatches+=Math.max(subjectNextLength,queryNextLength);
					} else {
						alignmentEncoding.addAll(ReadAlignment.encodePairwiseAlignment(alignedFragments));
						double mismatchesSegment = hamming.calculateDistance(alignedFragments[0], alignedFragments[1]); 
						numMismatches+=mismatchesSegment;
						if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligned fragments: \n"+alignedFragments[0]+"\n"+alignedFragments[1]+"\nSegment mismatches: "+mismatchesSegment+" total mismatches: "+numMismatches);
					}
					
				}
				nextMatchLength+=hitLength;
				coverageSharedKmers+=hitLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)hitLength*weight);
				subjectNext = kmerHit.getSubjectStart()+hitLength;
				queryNext = kmerHit.getQueryStart()+hitLength;
			} else {
				int kmerSubjectNext = kmerHit.getSubjectStart()+hitLength;
				int diffSubject = kmerSubjectNext - subjectNext;
				int kmerQueryNext = kmerHit.getQueryStart()+hitLength;
				int diffQuery = kmerQueryNext - queryNext;
				if (diffSubject>0 && diffSubject==diffQuery) {
					//Kmer consistent with current alignment. Augment match with difference
					nextMatchLength+=diffSubject;
					subjectNext = kmerSubjectNext;
					queryNext = kmerQueryNext;
					coverageSharedKmers+=Math.min(diffQuery, hitLength);
					double weight = kmerHit.getWeight();
					weightedCoverageSharedKmers+=((double)Math.min(diffQuery, hitLength)*weight);
				}
			}
			
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Processed Kmer hit at pos: "+kmerHit.getQueryStart()+" query next: "+queryNext+" subject next: "+subjectNext);
		}
		if(subjectNext==-1) {
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("WARN. Alignment could not be started for aln to: "+kmerHitsCluster.getSubjectIdx()+" query length: "+queryLength+" alnEncoding: "+alignmentEncoding);
			return null;
		}
		if(nextMatchLength>0) {
			alignmentEncoding.add(ReadAlignment.getAlnValue(nextMatchLength, ReadAlignment.ALIGNMENT_MATCH));
			nextMatchLength = 0;
		}
		int alnEnd = subjectNext;
		
		int remainder = query.length()-queryNext;
		boolean endAligned = remainder<=0;
		if(!endAligned && remainder+5 < maxLengthEndsPairwiseAlignment) {
			int end = Math.min(subjectNext+remainder+10, subject.length());
			if(subject.length()-subjectNext>=remainder) {
				String queryStr = queryS.substring(queryNext,query.length()).toString();
				String subjectStr = subject.subSequence(subjectNext,end).toString();
				if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligning end "+subjectStr+" of subject subsequence with total length: "+subject.length()+" to end "+queryStr+" of query with total length: "+query.length());
				String [] alignedFragments = alignFragments(queryStr, subjectStr, subjectIdx, queryLength, alignerEnd, 2);
				if(alignedFragments!=null) {
					int numUnalignedEndSubject = calculateUnaligned(subjectStr,alignedFragments[1]);
					String [] trimmedAlignedFragments = trimEndDeletions(alignedFragments);
					int trimmedBP = alignedFragments[0].length()-trimmedAlignedFragments[0].length();
					if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligned end fragments: \n"+alignedFragments[0]+"\n"+alignedFragments[1]+"\ntrimmed bp: "+trimmedBP+" Last bp1: "+alignedFragments[0].substring(alignedFragments[0].length()-10)+" lastbp2: "+alignedFragments[1].substring(alignedFragments[1].length()-10));
					alignedFragments = trimmedAlignedFragments;
					alignmentEncoding.addAll(ReadAlignment.encodePairwiseAlignment(alignedFragments));
					double mismatchesSegment = hamming.calculateDistance(alignedFragments[0], alignedFragments[1]); 
					numMismatches+=mismatchesSegment;
					alnEnd = end-numUnalignedEndSubject-trimmedBP;
					if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligned and trimmed end fragments: \n"+alignedFragments[0]+"\n"+alignedFragments[1]+"\nSegment mismatches: "+mismatchesSegment+" total mismatches: "+numMismatches+" aln end: "+alnEnd);
					endAligned = true;
				}	
			}
		}
		if(!endAligned) {
			//Ignore last bp
			alignmentEncoding.add(ReadAlignment.getAlnValue(remainder, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
		}
		ReadAlignment finalAlignment = new ReadAlignment(subjectIdx, alnStart+1, alnEnd, query.length(), 0);
		finalAlignment.setReadCharacters(query);
		finalAlignment.setAlignment(alignmentEncoding);
		String initialCIGAR = finalAlignment.getCigarString(); 
		finalAlignment.setNumMismatches(numMismatches);
		finalAlignment.collapseComplementaryIndels();
		if (subjectIdx == subjectIdxDebug  && queryLength==queryLengthDebug) System.out.println("FinalAlignment. CIGAR before: "+initialCIGAR+" Alignment after: "+finalAlignment);
		finalAlignment.setCoverageSharedKmers(coverageSharedKmers);
		finalAlignment.setWeightedCoverageSharedKmers((int)Math.round(weightedCoverageSharedKmers));
		//TODO: Define better alignment quality
		double covScore = 1.0*(queryNext-queryStart)/query.length();
		double estMaxMismatches = queryLength/5;
		double mismatchScore = Math.max(0, 1.0*(estMaxMismatches-numMismatches)/estMaxMismatches);
		finalAlignment.setAlignmentQuality((byte) Math.round(100*covScore*mismatchScore));
		finalAlignment.clipBorders(5);
		return finalAlignment;
	}
	private int calculateUnaligned(String sequence, String alignedSequence) {
		int n1 = sequence.length();
		int n2 = alignedSequence.length();
		int i=0;
		int j=0;
		while(i<n1 && j<n2) {
			if(alignedSequence.charAt(j)==LimitedSequence.GAP_CHARACTER) {
				j++;
			} else if (sequence.charAt(i)==alignedSequence.charAt(j)) {
				i++;
				j++;
			} else {
				throw new RuntimeException("Inconsistent base pair between sequence and alignment. Positions: "+i+" "+j+" sequences\n"+sequence+"\n"+alignedSequence);
			}
		}
		return n1-i;
	}
	private String[] trimStartDeletions(String[] alignedFragments) {
		int n = alignedFragments[0].length();
		int i;
		for(i=0;i<n && alignedFragments[0].charAt(i)==LimitedSequence.GAP_CHARACTER;i++);
		String[] answer = new String[2];
		answer[0]=alignedFragments[0].substring(i);
		answer[1]=alignedFragments[1].substring(i);
		return answer;
	}
	private String[] trimEndDeletions(String[] alignedFragments) {
		int n = alignedFragments[0].length();
		int i;
		for(i=n-1;i>=0 && alignedFragments[0].charAt(i)==LimitedSequence.GAP_CHARACTER;i--);
		String[] answer = new String[2];
		//System.out.println("Last bp1: "+alignedFragments[0].charAt(n-1)+" lastbp2: "+alignedFragments[1].charAt(n-1)+" i: "+i);
		answer[0]=alignedFragments[0].substring(0,i+1);
		answer[1]=alignedFragments[1].substring(0,i+1);
		return answer;
	}
	private String [] alignFragments(String queryStr, String subjectStr, int subjectIdx, int queryLength, PairwiseAligner aligner, int type) {
		int ql = queryStr.length();
		int sl = subjectStr.length();
		int minLength = Math.min(ql, sl);
		int maxLength = Math.max(ql, sl);
		if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligning segment of length "+ql+" of query to segment with length "+sl+" of subject \n"+queryStr+"\n"+subjectStr);
		PairwiseAlignment alignment = null;
		if(minLength<6 && maxLength>=minLength+11 && type == 1) {
			//Try to find a (quasi) perfect match
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Trying substring alignment");
			alignment = (new PairwiseAlignerSubstring()).calculateAlignment(queryStr, subjectStr);
		}
		if(alignment==null && minLength<0.1*maxLength) {
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Very large difference between segments. Trying naive alignment");
			alignment = (new PairwiseAlignerNaive(type==0)).calculateAlignment(queryStr, subjectStr);
		}
		if(alignment==null && maxLength < minLength+11 ) {
			//Try first with banded alignment if lengths are close
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Trying banded alignment");
			PairwiseAlignerStaticBanded alignerB = new PairwiseAlignerStaticBanded();
			alignerB.setForceStart2(type!=0);
			alignerB.setForceEnd2(type!=2);
			alignment = (alignerB.calculateAlignment(queryStr, subjectStr));
			if(alignment!=null) { 
				double mismatchesAln = alignment.getMismatches();
				if (type==0) {
					String[] seqs = {alignment.getAlignedSequence1(),alignment.getAlignedSequence2()}; 
					String [] trimmed = trimStartDeletions(seqs);
					int nt = alignment.getAlignedSequence1().length()-trimmed[0].length();
					mismatchesAln-=nt;
				} else if (type ==2) {
					String[] seqs = {alignment.getAlignedSequence1(),alignment.getAlignedSequence2()}; 
					String [] trimmed = trimEndDeletions(seqs);
					int nt = alignment.getAlignedSequence1().length()-trimmed[0].length();
					mismatchesAln-=nt;
				}
				if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Banded alignment done. Alignment\n"+alignment.getAlignedSequence1()+"\n"+alignment.getAlignedSequence2()+"\nNum mismatches: "+mismatchesAln);
				if(mismatchesAln>0.1*minLength) alignment = null;
			}
		}
		if(alignment==null) {
			alignment = aligner.calculateAlignment(queryStr,subjectStr);
		}
		if(alignment==null && trySimpleGap && maxLength<maxLengthFullPairwiseAlignment ) {
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Null alignment. Trying single gap alignment");
			alignment = (new PairwiseAlignerSimpleGap().calculateAlignment(queryStr, subjectStr));
			double mismatchesAln = alignment.getMismatches();
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Single gap alignment done. Num mismatches: "+mismatchesAln);
			if(mismatchesAln>0.2*minLength) alignment = null;
		}
		if(alignment==null) return null;
		String [] alignedFragments= {alignment.getAlignedSequence1(),alignment.getAlignedSequence2()};
		return alignedFragments;
	}
	public static int[] simulateAlignment(int subjectSeqIdx, int subjectLength, int querySeqIdx, int queryLength, UngappedSearchHitsCluster kmerHitsCluster) {
		int debugIdxS = -1;
		int debugIdxQ = -1;
		List<UngappedSearchHit> kmerHits = kmerHitsCluster.getHitsByQueryIdx();
		if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subject id "+subjectSeqIdx+" Subject length: "+subjectLength+". Query length: "+queryLength+" kmer hits: "+kmerHits.size()+ " cluster last "+kmerHitsCluster.getSubjectPredictedEnd());
		int coverageSharedKmers = 0;
		double weightedCoverageSharedKmers = 0;
		List<Integer> indelStarts = new ArrayList<Integer>();
		List<Integer> indelCalls = new ArrayList<Integer>();
		int initialNumIndels = 0;
		int subjectNext = -1;
		int queryNext = 0;
		for(UngappedSearchHit kmerHit:kmerHits) {
			int hitLength = kmerHit.getHitLength();
			if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subject id "+subjectSeqIdx+" Processing Kmer hit at pos: "+kmerHit.getQueryStart()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getSubjectStart());
			if(subjectNext==-1) {
				//Inconsistent kmer hit
				if (kmerHit.getSubjectStart()<kmerHitsCluster.getSubjectPredictedStart()) continue;
				coverageSharedKmers+=hitLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)hitLength*weight);
				if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" hit length: "+hitLength+" cov shared: "+coverageSharedKmers+" weight: "+weight+" wcov: "+weightedCoverageSharedKmers+" partial indels estimation: "+initialNumIndels);
				subjectNext = kmerHit.getSubjectStart()+hitLength;
				queryNext = kmerHit.getQueryStart()+hitLength;
			} else if(kmerHit.getQueryStart() > queryNext && subjectNext<kmerHit.getSubjectStart()) {
				//Kmer does not overlap with already aligned segments
				int subjectNextLength = kmerHit.getSubjectStart()-subjectNext;
				int queryNextLength = kmerHit.getQueryStart()-queryNext;
				int minLength = Math.min(subjectNextLength, queryNextLength);
				int maxLength = Math.max(subjectNextLength, queryNextLength);
				if(maxLength>minLength+3 && 0.9*maxLength>minLength) {
					//Possible invalid kmer hit. Delay alignment
					if (subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("Possible invalid kmer hit. Kmer hit at pos: "+kmerHit.getQueryStart()+" subject hit start: "+kmerHit.getSubjectStart()+" Subject length "+subjectNextLength+" query length "+queryNextLength);
					continue;
				}
				
				//Penalize up to 3 bp for each inconsistency
				//if(subjectNextLength!=queryNextLength) numIndels+=Math.abs(queryNextLength-subjectNextLength);
				int diff = Math.abs(queryNextLength-subjectNextLength);
				if(diff>2) {
					indelStarts.add(kmerHit.getQueryStart());
					indelCalls.add(queryNextLength-subjectNextLength);
					initialNumIndels+=diff;
				}
				//TODO: use 2*window length
				int nextCSK = Math.min(100, subjectNextLength) + hitLength;
				if(diff > 1) nextCSK = Math.max(hitLength, nextCSK-diff);
				coverageSharedKmers+=nextCSK;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)nextCSK*weight);
				if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subid "+subjectSeqIdx+" subnext: "+subjectNext+" subhit: "+kmerHit.getSubjectStart()+" subdist: "+subjectNextLength+" qnext: "+queryNext+" qhit: "+kmerHit.getQueryStart()+" qdist: "+queryNextLength+" cov shared: "+coverageSharedKmers+" weight: "+weight+" wcov: "+weightedCoverageSharedKmers+" partial indels estimation: "+initialNumIndels);
				subjectNext = kmerHit.getSubjectStart()+hitLength;
				queryNext = kmerHit.getQueryStart()+hitLength;
			} else {
				int kmerSubjectNext = kmerHit.getSubjectStart()+hitLength;
				int diffSubject = kmerSubjectNext - subjectNext;
				int kmerQueryNext = kmerHit.getQueryStart()+hitLength;
				int diffQuery = kmerQueryNext - queryNext;
				if (diffSubject>0 && diffSubject==diffQuery) {
					//Kmer consistent with current alignment. Augment match with difference
					subjectNext = kmerSubjectNext;
					queryNext = kmerQueryNext;
					coverageSharedKmers+=Math.min(diffQuery, hitLength);
					double weight = kmerHit.getWeight();
					weightedCoverageSharedKmers+=((double)Math.min(diffQuery, hitLength)*weight);
					if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subid "+subjectSeqIdx+" subnext: "+subjectNext+" subhit: "+kmerHit.getSubjectStart()+" qnext: "+queryNext+" qhit: "+kmerHit.getQueryStart()+" diff query: "+diffQuery+" hit length: "+hitLength+" cov shared: "+coverageSharedKmers+" weight: "+weight+" wcov: "+weightedCoverageSharedKmers+" partial indels estimation: "+initialNumIndels);
				} else {
					if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" inconsistent kmer alignment. diff query: "+diffQuery+" diffsubject "+diffSubject+" hit length: "+hitLength+" cov shared: "+coverageSharedKmers+" wcov: "+weightedCoverageSharedKmers+" partial indels estimation: "+initialNumIndels);
				}
				
			}
			
			if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ)  System.out.println("subject id "+subjectSeqIdx+" Processed Kmer hit at pos: "+kmerHit.getQueryStart()+" query next: "+queryNext+" subject next: "+subjectNext);
		}
		int numIndels = initialNumIndels;
		//int numIndels = calculateTotalIndels(subjectSeqIdx, querySeqIdx,indelStarts, indelCalls);
		//if(subjectSeqIdx==0) System.out.println("query length: "+queryLength+" cov: "+coverageSharedKmers+" wcov: "+weightedCoverageSharedKmers);
		int [] answer = {coverageSharedKmers, (int)Math.round(weightedCoverageSharedKmers),numIndels};
		return answer;
	}
	public static int calculateTotalIndels(int subjectSeqIdx, int querySeqIdx, List<Integer> indelStarts, List<Integer> indelCalls) {
		int calls = 0;
		int debugIdxS = -1;
		int debugIdxQ = -1;
		int n = indelCalls.size();
		if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("Calculating indels. calls: "+indelCalls);
		boolean [] processed = new boolean[n];
		Arrays.fill(processed, false);
		while (true) {
			//Find maximum
			int idxMax = 0;
			int maxAbs = 0;
			for(int i=0;i<n;i++) {
				if(processed[i]) continue;
				int next = Math.abs(indelCalls.get(i));
				if(next<10) processed[i] = true;
				if(next>maxAbs) {
					idxMax = i;
					maxAbs = next;
				}
			}
			if(maxAbs<6) break;
			//Try to find balancing indel
			int value = indelCalls.get(idxMax);
			int idxMinPair = idxMax;
			int minPair = maxAbs;
			for(int i=Math.max(0, idxMax-1);i<n && i<idxMax+2;i++) {
				int pairValue = Math.abs(value+indelCalls.get(i));
				if(pairValue<minPair) {
					idxMinPair = i;
					minPair = pairValue;
				}
			}
			processed[idxMax] = true;
			if(minPair>0.2*maxAbs) {
				if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("Keeping high value of site "+idxMax+" "+indelStarts.get(idxMax)+" value "+value+" min difference: "+minPair);
				continue;
			}
			processed[idxMinPair] = true;
			if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("Replacing values of sites "+idxMax+" "+indelStarts.get(idxMax)+" and "+idxMinPair+" "+indelStarts.get(idxMinPair)+" original "+value+" "+indelCalls.get(idxMinPair)+" difference: "+minPair);
			indelCalls.set(idxMax,minPair);
			indelCalls.set(idxMinPair,0);
		}
		for(int i=0;i<n;i++) {
			int nextPos = indelStarts.get(i);
			int nextCall = indelCalls.get(i);
			int absNextCall = Math.abs(nextCall);
			if(absNextCall<=1) continue;
			if(absNextCall<6 ) {
				calls+= Math.max(0, absNextCall-2);
				if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("Adding indels. Idx: "+i+" query pos: "+nextPos+" next call: "+nextCall+" indels: "+calls);
				continue;
			}
			//Try to find a balancing indel
			int j;
			int sumRange = 0;
			for(j=i;j<n && j<i+3;j++) {
				int nextCall2 = indelCalls.get(j);
				sumRange+=nextCall2;
				if(Math.abs(sumRange)< 0.2*absNextCall) {
					break;
				}
			}
			if(i>1 && j<n-2) {
				calls+=Math.max(0, Math.abs(sumRange)-2);
				if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("Adding indels group from "+i+" to "+j+". query pos: "+nextPos+" sum range. "+sumRange+" indels: "+calls);
			} else {
				if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("Ignoring indels group from "+i+" to "+j+". query pos: "+nextPos+" sum range. "+sumRange+" indels: "+calls);
			}
			
			i=j;
		}
		return calls;
	}
}