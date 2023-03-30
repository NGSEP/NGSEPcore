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
import ngsep.sequences.KmerSearchResultsCompressedTable;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;
import ngsep.sequences.ShortKmerCodesTable;

/**
 * @author Jorge Duitama
 */
public class MinimizersTableReadAlignmentAlgorithm implements ReadAlignmentAlgorithm {

	public static final int ALIGNMENT_ALGORITHM_AFFINE_GAP = LongReadsUngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_AFFINE_GAP;
	public static final int ALIGNMENT_ALGORITHM_DYNAMIC_KMERS = LongReadsUngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS;
	public static final int ALIGNMENT_ALGORITHM_SIMPLE_GAP = LongReadsUngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_SIMPLE_GAP;
	public static final int ALIGNMENT_ALGORITHM_NAIVE = LongReadsUngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_NAIVE;
	public static final int ALIGNMENT_ALGORITHM_SHORTREADS = 100;

	private Logger log = Logger.getLogger(MinimizersTableReadAlignmentAlgorithm.class.getName());
	
	private UngappedSearchHitsClusterAligner aligner;
	private int maxAlnsPerRead = 3;
	private double minWeightedCount = 5;
	private double minProportionBestCount = 0.2;
	private double minProportionReadLength = 0.01;
	
	private ReferenceGenome genome;
	private ShortKmerCodesTable kmerCodesTable;
	private boolean onlyPositiveStrand = false;
	
	public MinimizersTableReadAlignmentAlgorithm() {
		this (ALIGNMENT_ALGORITHM_AFFINE_GAP);
	}
	public MinimizersTableReadAlignmentAlgorithm(int alignmentAlgorithm) {
		if(alignmentAlgorithm == ALIGNMENT_ALGORITHM_SHORTREADS) {
			aligner = new ShortReadsUngappedSearchHitsClusterAligner();
		} else {
			aligner = new LongReadsUngappedSearchHitsClusterAligner(alignmentAlgorithm);
		}
	}
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public UngappedSearchHitsClusterAligner getAligner() {
		return aligner;
	}
	public void setAligner(UngappedSearchHitsClusterAligner aligner) {
		this.aligner = aligner;
	}
	public int getMaxAlnsPerRead() {
		return maxAlnsPerRead;
	}
	public void setMaxAlnsPerRead(int maxAlnsPerRead) {
		this.maxAlnsPerRead = maxAlnsPerRead;
	}
	
	public double getMinWeightedCount() {
		return minWeightedCount;
	}
	public void setMinWeightedCount(double minWeightedCount) {
		this.minWeightedCount = minWeightedCount;
	}
	public double getMinProportionBestCount() {
		return minProportionBestCount;
	}
	public void setMinProportionBestCount(double minProportionBestCount) {
		this.minProportionBestCount = minProportionBestCount;
	}
	public double getMinProportionReadLength() {
		return minProportionReadLength;
	}
	public void setMinProportionReadLength(double minProportionReadLength) {
		this.minProportionReadLength = minProportionReadLength;
	}
	public void loadGenome(ReferenceGenome genome, int kmerLength, int windowLength, int numThreads) {
		loadGenome(genome, kmerLength, windowLength, numThreads,false);
	}
	public void loadGenome(ReferenceGenome genome, int kmerLength, int windowLength, int numThreads, boolean buildKmersTable) {
		this.genome = genome;
		int n = genome.getNumSequences();
		int longestSequenceLengthMbp = 1+genome.getLongestSequenceLength()/1000000;
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
		//TODO: Fix introduced bias
		kmerCodesTable.setMaxHitsKmerCode(1000);
		kmerCodesTable.setLimitHitsPerSequence(1000);
		ThreadPoolManager poolTable = new ThreadPoolManager(numThreads, n);
		poolTable.setSecondsPerTask(longestSequenceLengthMbp);
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
		List<UngappedSearchHitsCluster> clusters = buildHitClusters(query,false);
		List<ReadAlignment> answer = buildAlignments(query, clusters);
		//System.out.println("Found "+answer.size()+" alignments");
		return answer;
	}
	public List<UngappedSearchHitsCluster> buildHitClusters (CharSequence query, boolean extensiveKmersSearch) {
		int queryLength = query.length();
		Map<Integer,List<UngappedSearchHit>> hitsByReference;
		if(extensiveKmersSearch) {
			Map<Integer, Long> codes = KmersExtractor.extractDNAKmerCodes(query.toString(), kmerCodesTable.getKmerLength() , 0, query.length());
			KmerSearchResultsCompressedTable results = kmerCodesTable.matchCompressed(-1, query.length(), codes, -1);
			hitsByReference = results.getAllHits();	
		}
		else hitsByReference = kmerCodesTable.match(-1,query);
		List<UngappedSearchHitsCluster> clusters = new ArrayList<UngappedSearchHitsCluster>();
		double minRawHitsSize = Math.max(2*minWeightedCount, minProportionReadLength*query.length());
		for (int sequenceIdx:hitsByReference.keySet()) {
			int sequenceLength = genome.getSequenceByIndex(sequenceIdx).getLength();
			List<UngappedSearchHit> totalHitsSubject = hitsByReference.get(sequenceIdx);
			//System.out.println("Reference id: "+sequenceIdx+" total raw hits: "+totalHitsSubject.size());
			//for(UngappedSearchHit hit: totalHitsSubject) if(hit.getQueryIdx()>21000 && hit.getQueryIdx()<27000) System.out.println("Next hit. "+hit.getQueryIdx()+" "+hit.getQuery()+" "+hit.getStart()+" "+hit.getWeight()+" "+(hit.getStart()-hit.getQueryIdx()));
			//for(UngappedSearchHit hit: totalHitsSubject) if(query.length()==11805 && hit.getStart()>0 && hit.getStart()<1000000) System.out.println("Next hit. "+hit.getQueryIdx()+" "+hit.getQuery()+" "+hit.getStart()+" "+hit.getWeight()+" "+(hit.getStart()-hit.getQueryIdx()));
			Collections.sort(totalHitsSubject, (h1,h2)->h1.getSubjectStart()-h2.getSubjectStart());
			List<UngappedSearchHit> rawClusterKmers = new ArrayList<UngappedSearchHit>();
			UngappedSearchHitsCluster cluster = null;
			for(UngappedSearchHit hit:totalHitsSubject) {
				if(hit.getWeight()<0.01) {
					//System.out.println("Hit with low weight. Pos: "+hit.getQueryStart()+" weight: "+hit.getWeight());
					continue;
				}
				if(cluster==null) {
					cluster = new UngappedSearchHitsCluster(queryLength, sequenceIdx, sequenceLength, hit);
				} else if (!cluster.addKmerHit(hit, 0)) {
					if (rawClusterKmers.size()>=minRawHitsSize) {
						List<UngappedSearchHitsCluster> regionClusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(queryLength, sequenceIdx, sequenceLength, rawClusterKmers, 0);
						//System.out.println("Qlen: "+query.length()+" next raw cluster inside "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+rawClusterKmers.size()+" subclusters "+regionClusters.size());
						clusters.addAll(regionClusters);
					} 
					//else System.out.println("Qlen: "+query.length()+" next raw small cluster discarded "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+rawClusterKmers.size()+" limit: "+minRawHitsSize);
					cluster = new UngappedSearchHitsCluster(queryLength, sequenceIdx, sequenceLength, hit);
					rawClusterKmers.clear();
				}
				rawClusterKmers.add(hit);
			}
			if(cluster!=null && rawClusterKmers.size()>=minRawHitsSize) {
				List<UngappedSearchHitsCluster> regionClusters = (new UngappedSearchHitsClusterBuilder()).clusterRegionKmerAlns(queryLength, sequenceIdx, sequenceLength, rawClusterKmers, 0);
				//System.out.println("Qlen: "+query.length()+" next raw cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+cluster.getNumDifferentKmers()+" subclusters "+regionClusters.size());
				clusters.addAll(regionClusters);
			}
			//System.out.println("Reference id: "+sequenceIdx+" total clusters: "+clusters.size());
		}
		//System.out.println("Final numnber of clusters: "+clusters.size());
		return clusters;
	}
	
	

	public void printClusters(List<UngappedSearchHitsCluster> clusters) {
		System.out.println("Clusters: "+clusters.size());
		for(UngappedSearchHitsCluster cluster:clusters) {
			System.out.println("kmers: "+cluster.getNumDifferentKmers()+" predicted limits: "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd()+" query limits "+cluster.getQueryEvidenceStart()+"-"+cluster.getQueryEvidenceEnd());
		}
		
	}
	public List<ReadAlignment> buildAlignments(CharSequence query, List<UngappedSearchHitsCluster> clusters) {
		Collections.sort(clusters, (o1,o2)-> ((int)o2.getWeightedCount())-((int)o1.getWeightedCount()));
		double maxCount = summarize(clusters);
		List<ReadAlignment> answer = new ArrayList<ReadAlignment>();
		//System.out.println("Filtering clusters. Max alns per read: "+maxAlnsPerRead);
		for (int i=0;i<clusters.size() && i<maxAlnsPerRead;i++) {
			UngappedSearchHitsCluster cluster = clusters.get(i);
			int sequenceIdx = cluster.getSubjectIdx();
			double wc = cluster.getWeightedCount();
			//System.out.println("Qlen: "+query.length()+" next cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" hits "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount());
			if(i>0 && (wc<minWeightedCount || wc<minProportionBestCount*maxCount)) break;
			QualifiedSequence refSeq = genome.getSequenceByIndex(sequenceIdx);
			ReadAlignment aln = aligner.buildAlignment(query, refSeq.getCharacters(), cluster);
			//System.out.println("Qlen: "+query.length()+" next cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" hits "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount()+" aln "+aln);
			if(aln!=null) {
				aln.setSequenceName(refSeq.getName());
				answer.add(aln);
			}
		}
		return answer;
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
}

