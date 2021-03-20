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
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMapAnalyzer;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.PairwiseAlignmentAffineGap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;

/**
 * @author Jorge Duitama
 */
public class MinimizersTableReadAlignmentAlgorithm implements ReadAlignmentAlgorithm {

	private Logger log = Logger.getLogger(MinimizersTableReadAlignmentAlgorithm.class.getName());
	private int maxLengthFullPairwiseAlignment = 4000;
	private int maxLengthEndsPairwiseAlignment = 500;
	private PairwiseAlignmentAffineGap alignerCenter = new PairwiseAlignmentAffineGap(maxLengthFullPairwiseAlignment+1);
	private PairwiseAlignmentAffineGap alignerStart = new PairwiseAlignmentAffineGap(maxLengthEndsPairwiseAlignment);
	private PairwiseAlignmentAffineGap alignerEnd = new PairwiseAlignmentAffineGap(maxLengthEndsPairwiseAlignment);
	private int maxAlnsPerRead = 3;
	private ReferenceGenome genome;
	private MinimizersTable minimizersTable;
	private boolean onlyPositiveStrand = false;
	private HammingSequenceDistanceMeasure hamming = new HammingSequenceDistanceMeasure();
	
	public MinimizersTableReadAlignmentAlgorithm() {
		alignerStart.setForceStart2(false);
		alignerEnd.setForceEnd2(false);
	}
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	
	public int getMaxAlnsPerRead() {
		return maxAlnsPerRead;
	}

	public void setMaxAlnsPerRead(int maxAlnsPerRead) {
		this.maxAlnsPerRead = maxAlnsPerRead;
	}

	public void loadGenome(ReferenceGenome genome, int kmerLength, int windowLength, int numThreads) throws InterruptedException {
		this.genome = genome;
		int n = genome.getNumSequences();
		//log.info("Calculating kmers distribution");
		KmersExtractor extractor = new KmersExtractor();
		extractor.setNumThreads(numThreads);
		log.info("Extracting kmers from reference sequence");
		extractor.processQualifiedSequences(genome.getSequencesList());
		//KmersMapAnalyzer analyzer = new KmersMapAnalyzer(extractor.getKmersMap(), true);
		log.info("Creating minimizers table for genome with "+n+" sequences loaded from file: "+genome.getFilename());
		//minimizersTable = new MinimizersTable(analyzer, kmerLength, windowLength);
		minimizersTable = new MinimizersTable(kmerLength, windowLength);
		minimizersTable.setKmersMap(extractor.getKmersMap());
		minimizersTable.setMaxCountPerMinimizer(20);
		minimizersTable.setKeepSingletons(true);
		minimizersTable.setLog(log);
		ThreadPoolManager poolMinimizers = new ThreadPoolManager(numThreads, n);
		poolMinimizers.setSecondsPerTask(60);
		for (int i=0;i<n;i++) {
			final int seqId = i;
			poolMinimizers.queueTask(()->minimizersTable.addSequence(seqId, genome.getSequenceCharacters(seqId)));
		}
		poolMinimizers.terminatePool();
		minimizersTable.calculateDistributionHits().printDistribution(System.err);
		log.info("Calculated minimizers. Total: "+minimizersTable.size());
	}
	
	
	public MinimizersTable getMinimizersTable() {
		return minimizersTable;
	}
	public void setMinimizersTable(ReferenceGenome genome, MinimizersTable minimizersTable) {
		this.genome = genome;
		this.minimizersTable = minimizersTable;
	}
	@Override
	public List<ReadAlignment> alignRead (RawRead read) {
		List<ReadAlignment> alignments = new ArrayList<>();
		String readSeq = read.getSequenceString();
		String qual = read.getQualityScores();
		String reverseQS = null;
		if(qual == null || qual.length()!=readSeq.length()) {
			qual = RawRead.generateFixedQSString('5', readSeq.length());
			reverseQS = qual;
		} else if (!onlyPositiveStrand) {
			reverseQS = new StringBuilder(qual).reverse().toString();
		}
		String reverseComplement = null;
		if(!onlyPositiveStrand) {
			reverseComplement = DNAMaskedSequence.getReverseComplement(readSeq).toString();	
		}
		
		alignments.addAll(alignQueryToReference(readSeq));
		//System.out.println("Read: "+read.getName()+" Forward inexact alignments: "+alignments.size());
		if(reverseComplement!=null) {
			List<ReadAlignment> alnsR = alignQueryToReference(reverseComplement);
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
		return alignments;
	}
	public List<ReadAlignment> alignQueryToReference(CharSequence query) {
		int queryLength = query.length();
		List<ReadAlignment> answer = new ArrayList<ReadAlignment>();
		Map<Integer,List<UngappedSearchHit>> hitsByReference = minimizersTable.match(-1,query);
		List<KmerHitsCluster> clusters = new ArrayList<KmerHitsCluster>();
		for (int sequenceIdx:hitsByReference.keySet()) {
			int sequenceLength = genome.getSequenceByIndex(sequenceIdx).getLength();
			List<UngappedSearchHit> totalHitsSubject = hitsByReference.get(sequenceIdx);
			//System.out.println("Reference id: "+sequenceIdx+" total raw hits: "+totalHitsSubject.size());
			//for(UngappedSearchHit hit: totalHitsSubject) if(hit.getQueryIdx()>21000 && hit.getQueryIdx()<27000) System.out.println("Next hit. "+hit.getQueryIdx()+" "+hit.getQuery()+" "+hit.getStart()+" "+hit.getWeight()+" "+(hit.getStart()-hit.getQueryIdx()));
			Collections.sort(totalHitsSubject, (h1,h2)->h1.getStart()-h2.getStart());
			List<UngappedSearchHit> rawClusterKmers = new ArrayList<UngappedSearchHit>();
			KmerHitsCluster cluster = null;
			for(UngappedSearchHit hit:totalHitsSubject) {
				if(hit.getWeight()<0.01) {
					//System.out.println("Hit with low weight. Pos: "+hit.getQueryIdx()+" weight: "+hit.getWeight());
					continue;
				}
				if(cluster==null) {
					cluster = new KmerHitsCluster(queryLength, sequenceLength, hit);
				} else if (!cluster.addKmerHit(hit, 0)) {
					if (rawClusterKmers.size()>=0.01*query.length()) {
						List<KmerHitsCluster> regionClusters = KmerHitsCluster.clusterRegionKmerAlns(queryLength, sequenceLength, rawClusterKmers, 0);
						//System.out.println("Qlen: "+query.length()+" next raw cluster inside "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+rawClusterKmers.size()+" subclusters "+regionClusters.size());
						clusters.addAll(regionClusters);
					}
					cluster = new KmerHitsCluster(queryLength, sequenceLength, hit);
					rawClusterKmers.clear();
				}
				rawClusterKmers.add(hit);
			}
			if(cluster!=null && rawClusterKmers.size()>=0.01*query.length()) {
				List<KmerHitsCluster> regionClusters = KmerHitsCluster.clusterRegionKmerAlns(queryLength, sequenceLength, rawClusterKmers, 0);
				//System.out.println("Qlen: "+query.length()+" next raw cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+cluster.getNumDifferentKmers()+" subclusters "+regionClusters.size());
				clusters.addAll(regionClusters);
			}
		}
		
		double maxCount = 0;
		for (KmerHitsCluster cluster:clusters) {
			cluster.summarize();
			//System.out.println("Qlen: "+query.length()+" next cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+cluster.getNumDifferentKmers());
			maxCount = Math.max(maxCount,cluster.getWeightedCount());
		}
		Collections.sort(clusters, (o1,o2)-> (int)(o2.getWeightedCount()-o1.getWeightedCount()));
		for (int i=0;i<clusters.size() && i<maxAlnsPerRead;i++) {
			KmerHitsCluster cluster = clusters.get(i);
			int sequenceIdx = cluster.getSubjectIdx();
			if(cluster.getWeightedCount()<0.2*maxCount) break;
			QualifiedSequence refSeq = genome.getSequenceByIndex(sequenceIdx);
			ReadAlignment aln = buildCompleteAlignment(sequenceIdx, refSeq.getCharacters(), query, cluster);
			//ReadAlignment aln = alignRead(sequenceIdx, refSeq.getCharacters(),query,subjectStart,subjectEnd,0.3);
			//System.out.println("Qlen: "+query.length()+" next cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" hits "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount()+" aln "+aln);
			if(aln!=null) {
				aln.setSequenceName(refSeq.getName());
				answer.add(aln);
			}
		}
		//System.out.println("Found "+answer.size()+" alignments");
		return answer;
	}

	

	public void printClusters(List<KmerHitsCluster> clusters) {
		System.out.println("Clusters: "+clusters.size());
		for(KmerHitsCluster cluster:clusters) {
			System.out.println("kmers: "+cluster.getNumDifferentKmers()+" predicted limits: "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd()+" query limits "+cluster.getQueryEvidenceStart()+"-"+cluster.getQueryEvidenceEnd());
		}
		
	}

	public ReadAlignment buildCompleteAlignment(int subjectIdx, CharSequence subject, CharSequence query, KmerHitsCluster kmerHitsCluster) { 
		List<UngappedSearchHit> kmerHits = kmerHitsCluster.getHitsByQueryIdx();
		int subjectNext = -1;
		short numMismatches = 0;
		int coverageSharedKmers = 0;
		double weightedCoverageSharedKmers = 0;
		//System.out.println("Subject length: "+subject.length()+". Query length: "+query.length()+" kmer hits: "+kmerHits.size()+" subject next: "+subjectNext+ " cluster last "+kmerHitsCluster.getSubjectPredictedEnd());
		int queryNext = 0;
		int alnStart = -1;
		int queryStart = -1;
		LinkedList<Integer> alignmentEncoding = new LinkedList<Integer>();
		int nextMatchLength = 0;
		for(UngappedSearchHit kmerHit:kmerHits) {
			//System.out.println("Processing Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getStart()+" cigar length: "+cigar.length());
			int kmerLength = kmerHit.getQuery().length();
			if(alnStart==-1) {
				//Inconsistent kmer hit
				if (kmerHit.getStart()<kmerHitsCluster.getSubjectPredictedStart()) continue;
				alnStart = kmerHit.getStart();
				queryStart = kmerHit.getQueryIdx();
				if(queryStart>0 && queryStart<maxLengthEndsPairwiseAlignment && queryStart<kmerHit.getStart()) {
					String queryStr = query.subSequence(0,queryStart).toString();
					alnStart = Math.max(0, kmerHit.getStart()-queryStart-5);
					queryStart=0;
					String subjectStr = subject.subSequence(alnStart,kmerHit.getStart()).toString();
					String [] alignedFragments = alignerStart.getAlignment(queryStr, subjectStr);
					alignmentEncoding.addAll(ReadAlignment.encodePairwiseAlignment(alignedFragments));
					numMismatches+=hamming.calculateDistance(alignedFragments[0], alignedFragments[1]);
				} else if (queryStart>0) alignmentEncoding.add(ReadAlignment.getAlnValue(queryStart, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
				nextMatchLength+=kmerLength;
				coverageSharedKmers+=kmerLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)kmerLength*weight);
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else if(kmerHit.getQueryIdx() >= queryNext && subjectNext<=kmerHit.getStart()) {
				//Kmer does not overlap with already aligned segments
				int subjectNextLength = kmerHit.getStart()-subjectNext;
				int queryNextLength = kmerHit.getQueryIdx()-queryNext;
				if(subjectNextLength==queryNextLength && (subjectNextLength<10 || subjectNextLength>maxLengthFullPairwiseAlignment)) {
					nextMatchLength+=subjectNextLength;
					if(subjectNextLength>0) numMismatches+=hamming.calculateDistance(subject.subSequence(subjectNext, kmerHit.getStart()), query.subSequence(queryNext, kmerHit.getQueryIdx()));
				} else {
					if(nextMatchLength>0 && (subjectNextLength>0 || queryNextLength>0)) {
						//System.out.println("Found internal segment for possible alignment. Subject length "+subjectStr.length()+" query length "+queryStr.length()+" current match length: "+nextMatchLength);
						alignmentEncoding.add(ReadAlignment.getAlnValue(nextMatchLength, ReadAlignment.ALIGNMENT_MATCH));
						nextMatchLength = 0;
					}
					if(subjectNextLength>0 && queryNextLength>0) {		
						if (subjectNextLength>maxLengthFullPairwiseAlignment || queryNextLength>maxLengthFullPairwiseAlignment) return null;
						String subjectStr = subject.subSequence(subjectNext,kmerHit.getStart()).toString();
						String queryStr = query.subSequence(queryNext,kmerHit.getQueryIdx()).toString();
						//if (subjectNextLength>10 || queryNextLength>10) System.out.println("Aligning segment of length "+subjectNextLength+" of subject with total length: "+subject.length()+" to segment with length "+queryNextLength+" of query with total length: "+query.length());
						String [] alignedFragments = alignerCenter.getAlignment(queryStr,subjectStr);
						alignmentEncoding.addAll(ReadAlignment.encodePairwiseAlignment(alignedFragments));
						numMismatches+=hamming.calculateDistance(alignedFragments[0], alignedFragments[1]);
					} else if (subjectNextLength>0) {
						alignmentEncoding.add(ReadAlignment.getAlnValue(subjectNextLength, ReadAlignment.ALIGNMENT_DELETION));
						numMismatches+=subjectNextLength;
					} else if (queryNextLength>0) {
						alignmentEncoding.add(ReadAlignment.getAlnValue(queryNextLength, ReadAlignment.ALIGNMENT_INSERTION));
						numMismatches+=queryNextLength;
					}
				}
				nextMatchLength+=kmerLength;
				coverageSharedKmers+=kmerLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)kmerLength*weight);
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
					coverageSharedKmers+=Math.min(diffQuery, kmerLength);
					double weight = kmerHit.getWeight();
					weightedCoverageSharedKmers+=((double)Math.min(diffQuery, kmerLength)*weight);
				}
			}
			
			//System.out.println("Processed Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext);
		}
		if(nextMatchLength>0) {
			alignmentEncoding.add(ReadAlignment.getAlnValue(nextMatchLength, ReadAlignment.ALIGNMENT_MATCH));
			nextMatchLength = 0;
		}
		int alnEnd = subjectNext;
		if(queryNext<query.length()) {
			int remainder = query.length()-queryNext;
			int end = Math.min(subjectNext+remainder+5, subject.length());
			if(subject.length()-subjectNext>=remainder && remainder<=maxLengthEndsPairwiseAlignment && (end-subjectNext)<=maxLengthEndsPairwiseAlignment) {
				String queryStr = query.subSequence(queryNext,query.length()).toString();
				String subjectStr = subject.subSequence(subjectNext,end).toString();
				//System.out.println("Aligning end of length "+subjectStr.length()+" of subject subsequence with total length: "+subject.length()+" to end with length "+queryStr.length()+" of query with total length: "+query.length());
				String [] alignedFragments = alignerEnd.getAlignment(queryStr, subjectStr);
				alignmentEncoding.addAll(ReadAlignment.encodePairwiseAlignment(alignedFragments));
				numMismatches+=hamming.calculateDistance(alignedFragments[0], alignedFragments[1]);
				alnEnd = end;
			} else {
				//Ignore last bp
				alignmentEncoding.add(ReadAlignment.getAlnValue(remainder, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
			}
		}
		ReadAlignment finalAlignment = new ReadAlignment(subjectIdx, alnStart+1, alnEnd, query.length(), 0);
		finalAlignment.setReadCharacters(query);
		finalAlignment.setAlignment(alignmentEncoding);
		finalAlignment.setNumMismatches(numMismatches);
		finalAlignment.setCoverageSharedKmers(coverageSharedKmers);
		finalAlignment.setWeightedCoverageSharedKmers((int)Math.round(weightedCoverageSharedKmers));
		//TODO: Define better alignment quality
		double cov = 1.0*(queryNext-queryStart)/query.length();
		finalAlignment.setAlignmentQuality((byte) Math.round(100*cov));
		finalAlignment.clipBorders(10);
		return finalAlignment;
	}
	public static int[] simulateAlignment(int subjectSeqIdx, int subjectLength, int querySeqIdx, int queryLength, KmerHitsCluster kmerHitsCluster) {
		int debugIdx = -1;
		List<UngappedSearchHit> kmerHits = kmerHitsCluster.getHitsByQueryIdx();
		if(querySeqIdx==debugIdx) System.out.println("subject id "+subjectSeqIdx+" Subject length: "+subjectLength+". Query length: "+queryLength+" kmer hits: "+kmerHits.size()+ " cluster last "+kmerHitsCluster.getSubjectPredictedEnd());
		int coverageSharedKmers = 0;
		double weightedCoverageSharedKmers = 0;
		int numIndels = 0;
		int subjectNext = -1;
		int queryNext = 0;
		for(UngappedSearchHit kmerHit:kmerHits) {
			if(querySeqIdx==debugIdx) System.out.println("subject id "+subjectSeqIdx+" Processing Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getStart());
			int kmerLength = kmerHit.getQuery().length();
			if(subjectNext==-1) {
				//Inconsistent kmer hit
				if (kmerHit.getStart()<kmerHitsCluster.getSubjectPredictedStart()) continue;
				coverageSharedKmers+=kmerLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)kmerLength*weight);
				if(querySeqIdx==debugIdx) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" kmerLength: "+kmerLength+" cov shared: "+coverageSharedKmers+" weight: "+weight+" wcov: "+weightedCoverageSharedKmers+" indels: "+numIndels);
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else if(kmerHit.getQueryIdx() >= queryNext && subjectNext<=kmerHit.getStart()) {
				//Kmer does not overlap with already aligned segments
				int subjectNextLength = kmerHit.getStart()-subjectNext;
				int queryNextLength = kmerHit.getQueryIdx()-queryNext;
				//Penalize up to 3 bp for each inconsistency
				//if(subjectNextLength!=queryNextLength) numIndels+=Math.abs(queryNextLength-subjectNextLength);
				if(subjectNextLength!=queryNextLength) numIndels+=Math.min(Math.abs(queryNextLength-subjectNextLength),3);
				coverageSharedKmers+=kmerLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)kmerLength*weight);
				if(querySeqIdx==debugIdx) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" kmerLength: "+kmerLength+" cov shared: "+coverageSharedKmers+" weight: "+weight+" wcov: "+weightedCoverageSharedKmers+" indels: "+numIndels);
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
					if(querySeqIdx==debugIdx) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" diff query: "+diffQuery+" kmerLength: "+kmerLength+" cov shared: "+coverageSharedKmers+" weight: "+weight+" wcov: "+weightedCoverageSharedKmers+" indels: "+numIndels);
				} else {
					//numIndels++;
					if(querySeqIdx==debugIdx) System.out.println("subject id "+subjectSeqIdx+" subject next: "+subjectNext+" inconsistent kmer alignment. diff query: "+diffQuery+" diffsubject "+diffSubject+" kmerLength: "+kmerLength+" cov shared: "+coverageSharedKmers+" wcov: "+weightedCoverageSharedKmers+" indels: "+numIndels);
				}
				
			}
			
			if(querySeqIdx==debugIdx)  System.out.println("subject id "+subjectSeqIdx+" Processed Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext);
		}
		//if(subjectSeqIdx==0) System.out.println("query length: "+queryLength+" cov: "+coverageSharedKmers+" wcov: "+weightedCoverageSharedKmers);
		int [] answer = {coverageSharedKmers, (int)Math.round(weightedCoverageSharedKmers),numIndels};
		return answer;
	}
}
