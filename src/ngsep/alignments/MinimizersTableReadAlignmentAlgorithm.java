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
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.ThreadPoolManager;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;
import ngsep.sequences.ShortKmerCodesTable;

/**
 * @author Jorge Duitama
 */
public class MinimizersTableReadAlignmentAlgorithm implements ReadAlignmentAlgorithm {

	public static final int ALIGNMENT_ALGORITHM_AFFINE_GAP = 1;
	public static final int ALIGNMENT_ALGORITHM_DYNAMIC_KMERS = 2;
	public static final int ALIGNMENT_ALGORITHM_SIMPLE_GAP = 3;
	public static final int ALIGNMENT_ALGORITHM_NAIVE = 4;
	private Logger log = Logger.getLogger(MinimizersTableReadAlignmentAlgorithm.class.getName());
	private HammingSequenceDistanceMeasure hamming = new HammingSequenceDistanceMeasure();
	private int maxLengthFullPairwiseAlignment = 4000;
	private int maxLengthEndsPairwiseAlignment = 500;
	private PairwiseAligner alignerCenter;
	private PairwiseAligner alignerStart;
	private PairwiseAligner alignerEnd;
	private int maxAlnsPerRead = 3;
	private double maxProportionBestCount = 0.2;
	private ReferenceGenome genome;
	private ShortKmerCodesTable kmerCodesTable;
	private boolean onlyPositiveStrand = false;
	
	public MinimizersTableReadAlignmentAlgorithm() {
		alignerCenter = new PairwiseAlignerAffineGap(maxLengthFullPairwiseAlignment+1);
		alignerStart = new PairwiseAlignerAffineGap(maxLengthEndsPairwiseAlignment);
		alignerEnd = new PairwiseAlignerAffineGap(maxLengthEndsPairwiseAlignment);
		((PairwiseAlignerAffineGap)alignerStart).setForceStart2(false);
		((PairwiseAlignerAffineGap)alignerEnd).setForceEnd2(false);
	}
	public MinimizersTableReadAlignmentAlgorithm(int alignmentAlgorithm) {
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
	
	public double getMaxProportionBestCount() {
		return maxProportionBestCount;
	}
	public void setMaxProportionBestCount(double maxProportionBestCount) {
		this.maxProportionBestCount = maxProportionBestCount;
	}
	
	public void loadGenome(ReferenceGenome genome, int kmerLength, int windowLength, int numThreads) {
		loadGenome(genome, kmerLength, windowLength, numThreads,false);
	}
	public void loadGenome(ReferenceGenome genome, int kmerLength, int windowLength, int numThreads, boolean buildKmersTable) {
		this.genome = genome;
		int n = genome.getNumSequences();
		//KmersMapAnalyzer analyzer = new KmersMapAnalyzer(extractor.getKmersMap(), true);
		log.info("Creating kmer codes table for genome with "+n+" sequences loaded from file: "+genome.getFilename());
		kmerCodesTable = new ShortKmerCodesTable(kmerLength, windowLength);
		
		//log.info("Calculating kmers distribution");
		if(buildKmersTable) {
			KmersExtractor extractor = new KmersExtractor();
			extractor.setLog(log);
			extractor.setNumThreads(numThreads);
			log.info("Extracting kmers from reference sequence");
			extractor.processQualifiedSequences(genome.getSequencesList());
			kmerCodesTable.setKmersMap(extractor.getKmersMap());
		}
		kmerCodesTable.setLog(log);
		kmerCodesTable.setMaxHitsKmerCode(1000);
		ThreadPoolManager poolTable = new ThreadPoolManager(numThreads, n);
		poolTable.setSecondsPerTask(60);
		for (int i=0;i<n;i++) {
			final int seqId = i;
			try {
				poolTable.queueTask(()->kmerCodesTable.addSequence(seqId, genome.getSequenceCharacters(seqId)));
			} catch (InterruptedException e) {
				e.printStackTrace();
				throw new RuntimeException("Concurrence error creating minimizers table",e);
			}
		}
		try {
			poolTable.terminatePool();
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException("Concurrence error creating codes table",e);
		}
		//minimizersTable.calculateDistributionHits().printDistribution(System.err);
		log.info("Calculated kmer codes. Total: "+kmerCodesTable.size());
	}
	
	public ShortKmerCodesTable getKmerCodesTable() {
		return kmerCodesTable;
	}
	public void setKmerCodesTable(ReferenceGenome genome, ShortKmerCodesTable kmerCodesTable) {
		this.genome = genome;
		this.kmerCodesTable = kmerCodesTable;
	}
	@Override
	public List<ReadAlignment> alignRead (RawRead read) {
		return alignRead((QualifiedSequence)read);
	}
	public List<ReadAlignment> alignRead (QualifiedSequence read) {
		List<ReadAlignment> alignments = new ArrayList<>();
		String readSeq = read.getCharacters().toString();
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
		Map<Integer,List<UngappedSearchHit>> hitsByReference = kmerCodesTable.match(-1,query);
		List<UngappedSearchHitsCluster> clusters = new ArrayList<UngappedSearchHitsCluster>();
		for (int sequenceIdx:hitsByReference.keySet()) {
			int sequenceLength = genome.getSequenceByIndex(sequenceIdx).getLength();
			List<UngappedSearchHit> totalHitsSubject = hitsByReference.get(sequenceIdx);
			//System.out.println("Reference id: "+sequenceIdx+" total raw hits: "+totalHitsSubject.size());
			//for(UngappedSearchHit hit: totalHitsSubject) if(hit.getQueryIdx()>21000 && hit.getQueryIdx()<27000) System.out.println("Next hit. "+hit.getQueryIdx()+" "+hit.getQuery()+" "+hit.getStart()+" "+hit.getWeight()+" "+(hit.getStart()-hit.getQueryIdx()));
			Collections.sort(totalHitsSubject, (h1,h2)->h1.getStart()-h2.getStart());
			List<UngappedSearchHit> rawClusterKmers = new ArrayList<UngappedSearchHit>();
			UngappedSearchHitsCluster cluster = null;
			for(UngappedSearchHit hit:totalHitsSubject) {
				if(hit.getWeight()<0.01) {
					//System.out.println("Hit with low weight. Pos: "+hit.getQueryIdx()+" weight: "+hit.getWeight());
					continue;
				}
				if(cluster==null) {
					cluster = new UngappedSearchHitsCluster(queryLength, sequenceLength, hit);
				} else if (!cluster.addKmerHit(hit, 0)) {
					if (rawClusterKmers.size()>=0.01*query.length()) {
						List<UngappedSearchHitsCluster> regionClusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(queryLength, sequenceLength, rawClusterKmers, 0);
						//System.out.println("Qlen: "+query.length()+" next raw cluster inside "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+rawClusterKmers.size()+" subclusters "+regionClusters.size());
						clusters.addAll(regionClusters);
					}
					cluster = new UngappedSearchHitsCluster(queryLength, sequenceLength, hit);
					rawClusterKmers.clear();
				}
				rawClusterKmers.add(hit);
			}
			if(cluster!=null && rawClusterKmers.size()>=0.01*query.length()) {
				List<UngappedSearchHitsCluster> regionClusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(queryLength, sequenceLength, rawClusterKmers, 0);
				//System.out.println("Qlen: "+query.length()+" next raw cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+cluster.getNumDifferentKmers()+" subclusters "+regionClusters.size());
				clusters.addAll(regionClusters);
			}
		}
		
		double maxCount = 0;
		for (UngappedSearchHitsCluster cluster:clusters) {
			cluster.summarize();
			maxCount = Math.max(maxCount,cluster.getWeightedCount());
			//System.out.println("Qlen: "+query.length()+" next cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+cluster.getNumDifferentKmers()+" count: "+cluster.getWeightedCount()+" maxCount: "+maxCount);
		}
		Collections.sort(clusters, (o1,o2)-> ((int)o2.getWeightedCount())-((int)o1.getWeightedCount()));
		for (int i=0;i<clusters.size() && i<maxAlnsPerRead;i++) {
			UngappedSearchHitsCluster cluster = clusters.get(i);
			int sequenceIdx = cluster.getSubjectIdx();
			if(cluster.getWeightedCount()<maxProportionBestCount*maxCount) break;
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
	
	public ReadAlignment buildCompleteAlignment(int subjectIdx, CharSequence subject, CharSequence query, UngappedSearchHitsCluster kmerHitsCluster) { 
		int subjectIdxDebug = -1;
		int queryLengthDebug = -1;
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
			if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Processing Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getStart());
			int kmerLength = kmerHit.getQuery().length();
			if(alnStart==-1) {
				//Inconsistent kmer hit
				if (kmerHit.getStart()<kmerHitsCluster.getSubjectPredictedStart()) continue;
				alnStart = kmerHit.getStart();
				queryStart = kmerHit.getQueryIdx();
				boolean startAligned = queryStart<=0;
				if(!startAligned && queryStart<kmerHit.getStart()) {
					String queryStr = queryS.substring(0,queryStart);
					int possibleAlnStart = Math.max(0, kmerHit.getStart()-queryStart-5);
					String subjectStr = subject.subSequence(possibleAlnStart,kmerHit.getStart()).toString();
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
				nextMatchLength+=kmerLength;
				coverageSharedKmers+=kmerLength;
				double weight = kmerHit.getWeight();
				weightedCoverageSharedKmers+=((double)kmerLength*weight);
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else if(kmerHit.getQueryIdx() > queryNext && subjectNext<kmerHit.getStart()) {
				//Kmer does not overlap with already aligned segments
				int subjectNextLength = kmerHit.getStart()-subjectNext;
				int queryNextLength = kmerHit.getQueryIdx()-queryNext;
				boolean goodMatch = subjectNextLength==queryNextLength && subjectNextLength<50;
				double hammingDistance = goodMatch?hamming.calculateDistance(subject.subSequence(subjectNext, kmerHit.getStart()), queryS.substring(queryNext, kmerHit.getQueryIdx())):0;
				goodMatch = goodMatch && hammingDistance<0.03*queryNextLength;
				if(goodMatch) {
					nextMatchLength+=subjectNextLength;
					numMismatches+=hammingDistance;
					if (subjectIdx == subjectIdxDebug  && queryLength==queryLengthDebug) System.out.println("Aligning equal length strings. Kmer hit at pos: "+kmerHit.getQueryIdx()+" subject hit start: "+kmerHit.getStart()+" Subject length "+subjectNextLength+" query length "+queryNextLength+" total mismatches: "+numMismatches);
				}
				else {
					int minLength = Math.min(subjectNextLength, queryNextLength);
					int maxLength = Math.max(subjectNextLength, queryNextLength);
					if(maxLength>minLength+3 && 0.95*maxLength>minLength) {
						//Possible invalid kmer hit. Delay alignment
						if (subjectIdx == subjectIdxDebug  && queryLength==queryLengthDebug) System.out.println("Possible invalid kmer hit. Kmer hit at pos: "+kmerHit.getQueryIdx()+" subject hit start: "+kmerHit.getStart()+" Subject length "+subjectNextLength+" query length "+queryNextLength);
						continue;
					}
					if(nextMatchLength>0) {
						if (subjectIdx == subjectIdxDebug  && queryLength==queryLengthDebug) System.out.println("Found internal segment for possible alignment. Kmer hit at pos: "+kmerHit.getQueryIdx()+" subject hit start: "+kmerHit.getStart()+" Subject length "+subjectNextLength+" query length "+queryNextLength+" current match length: "+nextMatchLength);
						alignmentEncoding.add(ReadAlignment.getAlnValue(nextMatchLength, ReadAlignment.ALIGNMENT_MATCH));
						nextMatchLength = 0;
					}
					String subjectStr = subject.subSequence(subjectNext,kmerHit.getStart()).toString();
					String queryStr = queryS.substring(queryNext,kmerHit.getQueryIdx()).toString();
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
					if (subjectIdx == subjectIdxDebug && queryLength==queryLengthDebug) System.out.println("Aligned fragments: \n"+alignedFragments[0]+"\n"+alignedFragments[1]+"\ntotal mismatches: "+numMismatches);
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

	public void printClusters(List<UngappedSearchHitsCluster> clusters) {
		System.out.println("Clusters: "+clusters.size());
		for(UngappedSearchHitsCluster cluster:clusters) {
			System.out.println("kmers: "+cluster.getNumDifferentKmers()+" predicted limits: "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd()+" query limits "+cluster.getQueryEvidenceStart()+"-"+cluster.getQueryEvidenceEnd());
		}
		
	}
}
