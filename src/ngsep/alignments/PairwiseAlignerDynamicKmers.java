package ngsep.alignments;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.sequences.AbstractLimitedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.SimpleEditDistanceMeasure;
import ngsep.sequences.UngappedSearchHit;

public class PairwiseAlignerDynamicKmers implements PairwiseAligner {
	@Override
	public String[] calculateAlignment(CharSequence sequence1, CharSequence sequence2) {
		int n1 = sequence1.length();
		int n2 = sequence2.length();
		if(n1<n2) {
			String[] alnRev = calculateAlignment(sequence2, sequence1);
			if(alnRev == null) return null;
			String [] answer = {alnRev[1].toString(),alnRev[0].toString()};
			return answer;
		}
		if(n2==0) return (new PairwiseAlignerNaive(false)).calculateAlignment(sequence1, sequence2);
		if(0.3*n1>n2 || n2<20) {
			//System.err.println("WARN: Unbalanced lengths for global alignment: "+n1+" "+n2);
			PairwiseAlignerSimpleGap aligner = new PairwiseAlignerSimpleGap();
			aligner.setForceStart1(false);
			aligner.setForceEnd1(false);
			return aligner.calculateAlignment(sequence1, sequence2);
		}
		int kmerLength = Math.max(5, Math.min(n1, n2)/50);
		kmerLength = Math.min(31, kmerLength);
		Map<Integer,Long> codesSubject = KmersExtractor.extractDNAKmerCodes(sequence1, kmerLength, 0, sequence1.length());
		Map<Integer,Long> codesQuery = KmersExtractor.extractDNAKmerCodes(sequence2, kmerLength, 0, sequence2.length());
		UngappedSearchHitsCluster bestCluster = findBestKmersCluster(n1, codesSubject, n2, codesQuery, kmerLength);
		if(bestCluster==null) {
			System.err.println("WARN: Null cluster for alignment of sequences with lengths: "+n1+" "+n2+" kmer length: "+kmerLength+" kmers: "+codesSubject.size()+" "+codesQuery.size()+" sequences\n"+sequence1+"\n"+sequence2);
			return null;
		}
		List<UngappedSearchHit> kmerHits = bestCluster.getHitsByQueryIdx();
		int subjectNext = 0;
		//System.out.println("S1 length: "+n1+". S2 length: "+n2+" kmer hits: "+kmerHits.size()+" subject next: "+subjectNext+" cluster predicted start: "+bestCluster.getSubjectPredictedStart());
		int queryNext = 0;
		int alnStart = -1;
		StringBuilder aln1 = new StringBuilder();
		StringBuilder aln2 = new StringBuilder();
		int nextMatchLength = 0;
		for(UngappedSearchHit kmerHit:kmerHits) {
			//System.out.println("Processing Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getStart());
			if(alnStart==-1) {
				//Inconsistent kmer hit
				if (kmerHit.getStart()<bestCluster.getSubjectPredictedStart()) continue;
				alnStart = kmerHit.getStart();
				String seq1Fragment = sequence1.subSequence(0,alnStart).toString();
				String seq2Fragment = sequence2.subSequence(0,kmerHit.getQueryIdx()).toString();
				if(seq1Fragment.length()>0 || seq2Fragment.length()>0) {
					String [] alignedFragments = calculateAlignment(seq1Fragment, seq2Fragment);
					if(alignedFragments==null) return null;
					aln1.append(alignedFragments[0]);
					aln2.append(alignedFragments[1]);
				}
				nextMatchLength+=kmerLength;
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else if(queryNext <= kmerHit.getQueryIdx() && subjectNext<=kmerHit.getStart()) {
				//Kmer does not overlap with already aligned segments
				int subjectNextLength = kmerHit.getStart()-subjectNext;
				int queryNextLength = kmerHit.getQueryIdx()-queryNext;
				if(subjectNextLength==queryNextLength) {
					nextMatchLength+=subjectNextLength;
				} else {
					if(nextMatchLength>0) {
						aln1.append(sequence1,subjectNext-nextMatchLength, subjectNext);
						aln2.append(sequence2,queryNext-nextMatchLength, queryNext);
						nextMatchLength = 0;
					}
					if(subjectNextLength>0 || queryNextLength>0) {
						String seq1Fragment = sequence1.subSequence(subjectNext,kmerHit.getStart()).toString();
						String seq2Fragment = sequence2.subSequence(queryNext,kmerHit.getQueryIdx()).toString();
						if (subjectNextLength>10 || queryNextLength>10) System.out.println("Aligning segment of length "+subjectNextLength+" of subject with total length: "+n1+" to segment with length "+queryNextLength+" of query with total length: "+n2);
						String [] alignedFragments = calculateAlignment(seq1Fragment,seq2Fragment);
						if(alignedFragments==null) return null;
						aln1.append(alignedFragments[0]);
						aln2.append(alignedFragments[1]);
					}
				}
				nextMatchLength+=kmerLength;
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else {
				int kmerSubjectNext = kmerHit.getStart()+kmerLength;
				int diffSubject = kmerSubjectNext - subjectNext;
				int kmerQueryNext = kmerHit.getQueryIdx()+kmerLength;
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
				String seq1Fragment = sequence1.subSequence(subjectNext,n1).toString();
				String seq2Fragment = sequence2.subSequence(queryNext,n2).toString();
				//if (subjectNextLength>10 || queryNextLength>10) System.out.println("Aligning segment of length "+subjectNextLength+" of subject with total length: "+subject.length()+" to segment with length "+queryNextLength+" of query with total length: "+query.length());
				String [] alignedFragments = calculateAlignment(seq1Fragment,seq2Fragment);
				if(alignedFragments==null) return null;
				aln1.append(alignedFragments[0]);
				aln2.append(alignedFragments[1]);
			}
		}
		String [] answer = {aln1.toString(),aln2.toString()};
		return answer;
	}
	
	public static UngappedSearchHitsCluster findBestKmersCluster (int subjectLength, Map<Integer, Long> codesSubject, int queryLength, Map<Integer, Long> codesQuery, int kmerLength) {
		List<UngappedSearchHit> initialKmerHits = alignKmerCodes(codesSubject, codesQuery, kmerLength);
		//System.out.println("Number of kmer hits: "+initialKmerHits.size());
		if(initialKmerHits.size()==0) return null;
		
		List<UngappedSearchHitsCluster> clusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(queryLength, subjectLength, initialKmerHits, 0);
		//System.out.println("Number of clusters: "+clusters.size());
		//printClusters(clusters);
		if(clusters.size()>1) Collections.sort(clusters, (o1,o2)->o2.getNumDifferentKmers()-o1.getNumDifferentKmers());	
		else if (clusters.size()==0) return null;
		return clusters.get(0);
	}
	private static List<UngappedSearchHit> alignKmerCodes(Map<Integer, Long> codesSubject, Map<Integer, Long> codesQuery, int kmerLength) {
		Map<Long,List<Integer>> reverseSubjectMap = new HashMap<Long, List<Integer>>();
		for(Map.Entry<Integer, Long> entry:codesSubject.entrySet()) {
			List<Integer> starts = reverseSubjectMap.computeIfAbsent(entry.getValue(), v->new ArrayList<Integer>());
			starts.add(entry.getKey());
		}
		List<UngappedSearchHit> initialKmerHits = new ArrayList<UngappedSearchHit>();
		for(int i:codesQuery.keySet()) {
			Long codeRead = codesQuery.get(i);
			List<Integer> subjectPosList = reverseSubjectMap.get(codeRead);
			if(subjectPosList==null) continue;
			for(int subjectPos:subjectPosList) {
				CharSequence kmerRead = new String(AbstractLimitedSequence.getSequence(codeRead, kmerLength, DNASequence.EMPTY_DNA_SEQUENCE));
				UngappedSearchHit hit = new UngappedSearchHit(kmerRead, 0 , subjectPos);
				hit.setQueryIdx(i);
				initialKmerHits.add(hit);
			}
		}
		return initialKmerHits;
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
			if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subject id "+subjectSeqIdx+" Processing Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getStart());
			int kmerLength = kmerHit.getQuery().length();
			if(subjectNext==-1) {
				//Inconsistent kmer hit
				if (kmerHit.getStart()<kmerHitsCluster.getSubjectPredictedStart()) continue;
				coverageSharedKmers+=kmerLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)kmerLength*weight);
				if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" kmerLength: "+kmerLength+" cov shared: "+coverageSharedKmers+" weight: "+weight+" wcov: "+weightedCoverageSharedKmers+" partial indels estimation: "+initialNumIndels);
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else if(kmerHit.getQueryIdx() >= queryNext && subjectNext<=kmerHit.getStart()) {
				//Kmer does not overlap with already aligned segments
				int subjectNextLength = kmerHit.getStart()-subjectNext;
				int queryNextLength = kmerHit.getQueryIdx()-queryNext;
				//Penalize up to 3 bp for each inconsistency
				//if(subjectNextLength!=queryNextLength) numIndels+=Math.abs(queryNextLength-subjectNextLength);
				int diff = Math.abs(queryNextLength-subjectNextLength);
				if(diff>2) {
					indelStarts.add(kmerHit.getQueryIdx());
					indelCalls.add(queryNextLength-subjectNextLength);
					initialNumIndels+=diff;
				}
				coverageSharedKmers+=kmerLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)kmerLength*weight);
				if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" kmerLength: "+kmerLength+" cov shared: "+coverageSharedKmers+" weight: "+weight+" wcov: "+weightedCoverageSharedKmers+" partial indels estimation: "+initialNumIndels);
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else {
				int kmerSubjectNext = kmerHit.getStart()+kmerLength;
				int diffSubject = kmerSubjectNext - subjectNext;
				int kmerQueryNext = kmerHit.getQueryIdx()+kmerLength;
				int diffQuery = kmerQueryNext - queryNext;
				if (diffSubject>0 && diffSubject==diffQuery) {
					//Kmer consistent with current alignment. Augment match with difference
					subjectNext = kmerSubjectNext;
					queryNext = kmerQueryNext;
					coverageSharedKmers+=Math.min(diffQuery, kmerLength);
					double weight = kmerHit.getWeight();
					weightedCoverageSharedKmers+=((double)Math.min(diffQuery, kmerLength)*weight);
					if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" diff query: "+diffQuery+" kmerLength: "+kmerLength+" cov shared: "+coverageSharedKmers+" weight: "+weight+" wcov: "+weightedCoverageSharedKmers+" partial indels estimation: "+initialNumIndels);
				} else {
					if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" inconsistent kmer alignment. diff query: "+diffQuery+" diffsubject "+diffSubject+" kmerLength: "+kmerLength+" cov shared: "+coverageSharedKmers+" wcov: "+weightedCoverageSharedKmers+" partial indels estimation: "+initialNumIndels);
				}
				
			}
			
			if(subjectSeqIdx==debugIdxS && querySeqIdx==debugIdxQ)  System.out.println("subject id "+subjectSeqIdx+" Processed Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext);
		}
		int numIndels = calculateTotalIndels(subjectSeqIdx, querySeqIdx,indelStarts, indelCalls);
		//if(subjectSeqIdx==0) System.out.println("query length: "+queryLength+" cov: "+coverageSharedKmers+" wcov: "+weightedCoverageSharedKmers);
		int [] answer = {coverageSharedKmers, (int)Math.round(weightedCoverageSharedKmers),numIndels};
		return answer;
	}
	private static int calculateTotalIndels(int subjectSeqIdx, int querySeqIdx, List<Integer> indelStarts, List<Integer> indelCalls) {
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
