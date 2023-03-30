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
import java.util.List;

import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.UngappedSearchHit;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class LongReadsUngappedSearchHitsClusterAligner implements UngappedSearchHitsClusterAligner {
	public static final int ALIGNMENT_ALGORITHM_AFFINE_GAP = 1;
	public static final int ALIGNMENT_ALGORITHM_DYNAMIC_KMERS = 2;
	public static final int ALIGNMENT_ALGORITHM_SIMPLE_GAP = 3;
	public static final int ALIGNMENT_ALGORITHM_NAIVE = 4;
	private HammingSequenceDistanceMeasure hamming = new HammingSequenceDistanceMeasure();
	private int maxLengthFullPairwiseAlignment = 4000;
	private int maxLengthEndsPairwiseAlignment = 500;
	private PairwiseAligner alignerCenter;
	private PairwiseAligner alignerStart;
	private PairwiseAligner alignerEnd;
	
	public LongReadsUngappedSearchHitsClusterAligner () {
		this(ALIGNMENT_ALGORITHM_AFFINE_GAP);
	}
	public LongReadsUngappedSearchHitsClusterAligner (int alignmentAlgorithm) {
		if(alignmentAlgorithm == ALIGNMENT_ALGORITHM_AFFINE_GAP) {
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
		if(alignmentAlgorithm == ALIGNMENT_ALGORITHM_SIMPLE_GAP) {
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
		int subjectIdxDebug = -1;
		int queryLengthDebug = -1;
		int subjectIdx = kmerHitsCluster.getSubjectIdx();
		List<UngappedSearchHit> kmerHits = kmerHitsCluster.getHitsByQueryIdx();
		String queryS = query.toString();
		int queryLength= query.length();
		
		int subjectNext = -1;
		short numMismatches = 0;
		int coverageSharedKmers = 0;
		double weightedCoverageSharedKmers = 0;
		if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Subject length: "+subject.length()+". Query length: "+query.length()+" kmer hits: "+kmerHits.size()+" subject next: "+subjectNext+ " cluster last "+kmerHitsCluster.getSubjectPredictedEnd());
		int queryNext = 0;
		int alnStart = -1;
		int queryStart = -1;
		List<Integer> alignmentEncoding = new ArrayList<Integer>();
		int nextMatchLength = 0;
		for(UngappedSearchHit kmerHit:kmerHits) {
			int hitLength = kmerHit.getHitLength();
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Processing Kmer hit at pos: "+kmerHit.getQueryStart()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getSubjectStart());
			if(alnStart==-1) {
				//Inconsistent kmer hit
				if (kmerHit.getSubjectStart()<kmerHitsCluster.getSubjectPredictedStart()) continue;
				alnStart = kmerHit.getSubjectStart();
				queryStart = kmerHit.getQueryStart();
				boolean startAligned = queryStart<=0;
				if(!startAligned && queryStart<kmerHit.getSubjectStart()) {
					String queryStr = queryS.substring(0,queryStart);
					int possibleAlnStart = Math.max(0, kmerHit.getSubjectStart()-queryStart-5);
					String subjectStr = subject.subSequence(possibleAlnStart,kmerHit.getSubjectStart()).toString();
					if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Hit start. Query segment: "+queryStr+" subject segment: "+subjectStr);
					String [] alignedFragments = null;
					if(queryStr.length()<=5 || subjectStr.length()<=5) {
						alignedFragments = (new PairwiseAlignerNaive(true)).calculateAlignment(queryStr, subjectStr);
					} else if(queryStr.length()<maxLengthEndsPairwiseAlignment && subjectStr.length()<maxLengthEndsPairwiseAlignment){
						alignedFragments = alignerStart.calculateAlignment(queryStr, subjectStr);
					}
					if(alignedFragments!=null) {
						alignmentEncoding.addAll(ReadAlignment.encodePairwiseAlignment(alignedFragments));
						numMismatches+=hamming.calculateDistance(alignedFragments[0], alignedFragments[1]);
						startAligned = true;
						queryStart=0;
						alnStart = possibleAlnStart;
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
				boolean goodMatch = subjectNextLength==queryNextLength && subjectNextLength<50;
				double hammingDistance = goodMatch?hamming.calculateDistance(subject.subSequence(subjectNext, kmerHit.getSubjectStart()), queryS.substring(queryNext, kmerHit.getQueryStart())):0;
				goodMatch = goodMatch && hammingDistance<0.03*queryNextLength;
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
					if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligning segment of length "+subjectNextLength+" of subject with total length: "+subject.length()+" to segment with length "+queryNextLength+" of query with total length: "+query.length()+"\n"+subjectStr+"\n"+queryStr);
					String [] alignedFragments = alignerCenter.calculateAlignment(queryStr,subjectStr);
					if(alignedFragments==null && minLength<0.1*maxLength ) {
						//Possible large indel event
						if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Null alignment. Trying naive alignment");
						alignedFragments = (new PairwiseAlignerNaive(true)).calculateAlignment(queryStr, subjectStr);
					}
					if(alignedFragments==null ) {
						if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Null alignment. Trying single gap alignment");
						alignedFragments = (new PairwiseAlignerSimpleGap().calculateAlignment(queryStr, subjectStr));
						numMismatches+=hamming.calculateDistance(alignedFragments[0], alignedFragments[1]);
						if(numMismatches>0.2*minLength) alignedFragments = null; 
					}
					if(alignedFragments==null) {
						if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Null alignment. Trying default alignment");
						if(maxLength>0.2*query.length()) return null;
						alignmentEncoding.add(ReadAlignment.getAlnValue(minLength, ReadAlignment.ALIGNMENT_MISMATCH));
						if(subjectNextLength>queryNextLength) alignmentEncoding.add(ReadAlignment.getAlnValue(subjectNextLength-queryNextLength, ReadAlignment.ALIGNMENT_DELETION));
						else if (subjectNextLength<queryNextLength) alignmentEncoding.add(ReadAlignment.getAlnValue(queryNextLength-subjectNextLength, ReadAlignment.ALIGNMENT_INSERTION));
						numMismatches+=Math.max(subjectNextLength,queryNextLength);
					} else {
						alignmentEncoding.addAll(ReadAlignment.encodePairwiseAlignment(alignedFragments));
						numMismatches+=hamming.calculateDistance(alignedFragments[0], alignedFragments[1]);
					}
					if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) {
						if(alignedFragments==null) System.out.println("Default alignmet for query coords "+queryNext+" "+kmerHit.getQueryStart()+" length: "+queryStr.length()+" subject coords: "+subjectNext+" " +kmerHit.getSubjectStart());
						else System.out.println("Aligned fragments: \n"+alignedFragments[0]+"\n"+alignedFragments[1]+"\ntotal mismatches: "+numMismatches);
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
			
			//System.out.println("Processed Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext);
		}
		if(nextMatchLength>0) {
			alignmentEncoding.add(ReadAlignment.getAlnValue(nextMatchLength, ReadAlignment.ALIGNMENT_MATCH));
			nextMatchLength = 0;
		}
		int alnEnd = subjectNext;
		
		int remainder = query.length()-queryNext;
		boolean endAligned = remainder<=0;
		if(!endAligned && remainder+5 < maxLengthEndsPairwiseAlignment) {
			int end = Math.min(subjectNext+remainder+5, subject.length());
			if(subject.length()-subjectNext>=remainder) {
				String queryStr = queryS.substring(queryNext,query.length()).toString();
				String subjectStr = subject.subSequence(subjectNext,end).toString();
				if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligning end "+subjectStr+" of subject subsequence with total length: "+subject.length()+" to end "+queryStr+" of query with total length: "+query.length());
				String [] alignedFragments = alignerEnd.calculateAlignment(queryStr, subjectStr);
				if(alignedFragments!=null) {
					alignmentEncoding.addAll(ReadAlignment.encodePairwiseAlignment(alignedFragments));
					numMismatches+=hamming.calculateDistance(alignedFragments[0], alignedFragments[1]);
					alnEnd = end;
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
		double cov = 1.0*(queryNext-queryStart)/query.length();
		finalAlignment.setAlignmentQuality((byte) Math.round(100*cov));
		finalAlignment.clipBorders(5);
		return finalAlignment;
	}
}
