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
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.assembly.AlignmentConstantGap;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.PairwiseAlignmentAffineGap;
import ngsep.sequences.QualifiedSequence;

/**
 * @author Jorge Duitama
 */
public class LongReadsAligner {

	private Logger log = Logger.getLogger(LongReadsAligner.class.getName());
	private int maxLengthFullPairwiseAlignment = 4000;
	private PairwiseAlignmentAffineGap aligner = null;
	private ReferenceGenome genome;
	private MinimizersTable minimizersTable;
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	public void loadGenome(ReferenceGenome genome, int kmerLength, int windowLength) {
		this.genome = genome;
		int n = genome.getNumSequences();
		log.info("Creating minimizers table for genome with "+n+" sequences loaded from file: "+genome.getFilename());
		minimizersTable = new MinimizersTable(kmerLength, windowLength);
		minimizersTable.setMaxAbundanceMinimizer(20);
		minimizersTable.setSaveRepeatedMinimizersWithinSequence(true);
		minimizersTable.setLog(log);
		for (int i=0;i<n;i++) {
			minimizersTable.addSequence(i, genome.getSequenceCharacters(i));
		}
		minimizersTable.calculateDistributionHits().printDistribution(System.out);
		minimizersTable.clearOverrepresentedMinimizers();
		log.info("Calculated minimizers. Total: "+minimizersTable.getTotalMinimizers());
	}
	public List<ReadAlignment> alignQueryToReference(CharSequence query) {
		List<ReadAlignment> answer = new ArrayList<ReadAlignment>();
		Map<Integer,List<UngappedSearchHit>> hitsByReference = minimizersTable.match(query);
		List<KmerHitsCluster> clusters = new ArrayList<KmerHitsCluster>();
		for (int sequenceIdx:hitsByReference.keySet()) {
			List<UngappedSearchHit> totalHitsSubject = hitsByReference.get(sequenceIdx);
			Collections.sort(totalHitsSubject, (h1,h2)->h1.getStart()-h2.getStart());
			KmerHitsCluster cluster = null;
			for(UngappedSearchHit hit:totalHitsSubject) {
				if(hit.getTotalHitsQuery()>5) continue;
				if(cluster==null) {
					cluster = new KmerHitsCluster(query, hit);
				} else if (!cluster.addKmerHit(hit, 0)) {
					if (cluster.getNumDifferentKmers()>=0.01*query.length()) {
						List<KmerHitsCluster> regionClusters = KmerHitsCluster.clusterRegionKmerAlns(query, cluster.getHitsByQueryIdx(), 0.3);
						//System.out.println("Qlen: "+query.length()+" next raw cluster "+cluster.getSequenceIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" hits: "+cluster.getNumDifferentKmers()+" subclusters "+regionClusters.size());
						clusters.addAll(regionClusters);
					}
					cluster = new KmerHitsCluster(query, hit);
				}
			}
			if(cluster!=null && cluster.getNumDifferentKmers()>=0.01*query.length()) {
				List<KmerHitsCluster> regionClusters = KmerHitsCluster.clusterRegionKmerAlns(query, cluster.getHitsByQueryIdx(), 0.3);
				//System.out.println("Qlen: "+query.length()+" next raw cluster "+cluster.getSequenceIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" hits: "+cluster.getNumDifferentKmers()+" subclusters "+regionClusters.size());
				clusters.addAll(regionClusters);
			}
		}
		
		double maxCount = 0;
		for (KmerHitsCluster cluster:clusters) {
			cluster.summarize(1,false);
			maxCount = Math.max(maxCount,cluster.getWeightedCount());
		}
		Collections.sort(clusters, (o1,o2)-> (int)(o2.getWeightedCount()-o1.getWeightedCount()));
		//TODO: Parameter
		for (int i=0;i<clusters.size() && i<5;i++) {
			KmerHitsCluster cluster = clusters.get(i);
			int sequenceIdx = cluster.getSequenceIdx();
			if(cluster.getWeightedCount()<0.2*maxCount) break;
			QualifiedSequence refSeq = genome.getSequenceByIndex(sequenceIdx);
			ReadAlignment aln = buildCompleteAlignment(sequenceIdx, refSeq.getCharacters(), query, cluster);
			//ReadAlignment aln = alignRead(sequenceIdx, refSeq.getCharacters(),query,subjectStart,subjectEnd,0.3);
			//System.out.println("Qlen: "+query.length()+" next cluster "+cluster.getSequenceIdx()+": "+subjectStart+" "+subjectEnd+" hits "+cluster.getNumDifferentKmers()+" weighted count: "+cluster.getWeightedCount()+" aln "+aln);
			if(aln!=null) {
				aln.setSequenceName(refSeq.getName());
				answer.add(aln);
			}
			
		}
		
		//System.out.println("Found "+answer.size()+" alignments");
		return answer;
	}

	public Map<CharSequence, Integer> extractUniqueKmers(CharSequence sequence, int start, int end) {
		Map<Integer, CharSequence> rawKmers = KmersExtractor.extractKmersAsMap(sequence, 15, 1, start, end, true, true, true);
		Map<CharSequence, Integer> answer = new LinkedHashMap<CharSequence, Integer>();
		Map<CharSequence, Integer> reverseMap = new HashMap<CharSequence,Integer>();
		Set<Integer> multiple = new HashSet<>();
		for(int i=start;i<end;i++) {
			CharSequence kmer = rawKmers.get(i);
			if(kmer == null) {
				multiple.add(i);
				continue;
			}
			Integer previousStart = reverseMap.get(kmer);
			if(previousStart!=null) {
				multiple.add(i);
				continue;
			}
			reverseMap.put(kmer,i);
		}
		for(int i=start;i<end;i++) {
			CharSequence kmer = rawKmers.get(i);
			if(kmer!=null && !multiple.contains(i)) {
				answer.put(kmer.toString(),i);
			}
			
		}
		return answer;
	}

	public ReadAlignment alignRead(int subjectIdx, CharSequence subject, CharSequence read, int start, int end, double minQueryCoverage) {
		Map<CharSequence, Integer> uniqueKmersSubject = extractUniqueKmers(subject,start,end);
		//System.out.println("Number of unique k-mers subject: "+uniqueKmersSubject.size());
		return alignRead(subjectIdx, subject, read, uniqueKmersSubject, minQueryCoverage);
	}
	public ReadAlignment alignRead(int subjectIdx, CharSequence subject, CharSequence read, Map<CharSequence, Integer> uniqueKmersSubject, double minQueryCoverage) {
		Map<CharSequence, Integer> uniqueKmersRead = extractUniqueKmers(read,0,read.length());
		//System.out.println("Number of unique k-mers read: "+uniqueKmersRead.size());
		List<UngappedSearchHit> initialKmerHits = alignUniqueKmers(-1,subject.length(),uniqueKmersSubject, uniqueKmersRead);
		if(initialKmerHits.size()==0) return null;
		List<KmerHitsCluster> clusters = KmerHitsCluster.clusterRegionKmerAlns(read, initialKmerHits, minQueryCoverage);
		//printClusters(clusters);
		if(clusters.size()>1) {
			Collections.sort(clusters, (o1,o2)->o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
			KmerHitsCluster c1 = clusters.get(0);
			KmerHitsCluster c2 = clusters.get(1);
			int overlap = GenomicRegionSpanComparator.getInstance().getSpanLength(c1.getSubjectPredictedStart(), c1.getSubjectPredictedEnd(), c2.getSubjectPredictedStart(), c2.getSubjectPredictedEnd());
			int c1Length = c1.getSubjectPredictedEnd()-c1.getSubjectPredictedStart();
			int c2Length = c2.getSubjectPredictedEnd()-c2.getSubjectPredictedStart();
			if((overlap <0.9*c1Length || overlap < 0.9*c2Length) && c1.getNumDifferentKmers()<0.9*initialKmerHits.size()) {
				return null;
			}	
		} else if (clusters.size()==0) return null;
		KmerHitsCluster bestCluster = clusters.get(0);
		//System.out.println("Number of clusters: "+clusters.size()+" best cluster kmers: "+bestCluster.getNumDifferentKmers()+" first "+bestCluster.getFirst()+" last "+bestCluster.getLast());
		return buildCompleteAlignment(subjectIdx, subject, read, bestCluster);
	}
	private List<UngappedSearchHit> alignUniqueKmers(int subjectIdx, int subjectLength, Map<CharSequence, Integer> uniqueKmersSubject, Map<CharSequence, Integer> uniqueKmersQuery) {
		List<UngappedSearchHit> initialKmerHits = new ArrayList<UngappedSearchHit>();
		for(CharSequence kmerRead:uniqueKmersQuery.keySet()) {
			Integer subjectPos = uniqueKmersSubject.get(kmerRead);
			if(subjectPos==null) continue;
			UngappedSearchHit hit = new UngappedSearchHit(kmerRead, subjectIdx , subjectPos);
			hit.setSequenceLength(subjectLength);
			hit.setQueryIdx(uniqueKmersQuery.get(kmerRead));
			initialKmerHits.add(hit);
		}
		return initialKmerHits;
	}

	public void printClusters(List<KmerHitsCluster> clusters) {
		System.out.println("Clusters: "+clusters.size());
		for(KmerHitsCluster cluster:clusters) {
			System.out.println("kmers: "+cluster.getNumDifferentKmers()+" predicted limits: "+cluster.getSubjectPredictedStart()+" - "+cluster.getSubjectPredictedEnd()+" query limits "+cluster.getQueryEvidenceStart()+"-"+cluster.getQueryEvidenceEnd());
		}
		
	}

	private ReadAlignment buildCompleteAlignment(int subjectIdx, CharSequence subject, CharSequence query, KmerHitsCluster kmerHitsCluster) {
		if (aligner == null) aligner = new PairwiseAlignmentAffineGap(1, 2, 1, 1, 4100);
		List<UngappedSearchHit> kmerHits = kmerHitsCluster.getHitsByQueryIdx();
		int subjectNext = Math.max(0, kmerHitsCluster.getSubjectPredictedStart());
		//System.out.println("Subject length: "+subject.length()+". Query length: "+query.length()+" kmer hits: "+kmerHits.size()+" subject next: "+subjectNext+ " cluster last "+kmerHitsCluster.getSubjectPredictedEnd());
		int queryNext = 0;
		int alnStart = -1;
		int queryStart = -1;
		char matchOp = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		char insertionOp = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION);
		char deletionOp = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_DELETION);
		char softClipOp = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_SKIPFROMREAD);
		StringBuilder cigar = new StringBuilder();
		int nextMatchLength = 0;
		for(UngappedSearchHit kmerHit:kmerHits) {
			//System.out.println("Processing Kmer hit at pos: "+kmerHit.getQueryIdx()+" query next: "+queryNext+" subject next: "+subjectNext+" subject hit start: "+kmerHit.getStart()+" cigar length: "+cigar.length());
			int kmerLength = kmerHit.getQuery().length();
			if(alnStart==-1) {
				alnStart = kmerHit.getStart();
				queryStart = kmerHit.getQueryIdx();
				if(kmerHit.getQueryIdx()>0) cigar.append(""+kmerHit.getQueryIdx()+""+softClipOp);
				nextMatchLength+=kmerLength;
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else if(kmerHit.getQueryIdx() >= queryNext && subjectNext<=kmerHit.getStart()) {
				//Kmer does not overlap with already aligned segments
				int subjectNextLength = kmerHit.getStart()-subjectNext;
				int queryNextLength = kmerHit.getQueryIdx()-queryNext;
				if(subjectNextLength==queryNextLength && (subjectNextLength<10 || subjectNextLength>maxLengthFullPairwiseAlignment)) {
					nextMatchLength+=subjectNextLength;
				} else {
					if(nextMatchLength>0 && (subjectNextLength>0 || queryNextLength>0)) {
						//System.out.println("Found internal segment for possible alignment. Subject length "+subjectStr.length()+" query length "+queryStr.length()+" current match length: "+nextMatchLength);
						cigar.append(""+nextMatchLength+""+matchOp);
						nextMatchLength = 0;
					}
					if(subjectNextLength>0 && queryNextLength>0) {		
						if (subjectNextLength>maxLengthFullPairwiseAlignment || queryNextLength>maxLengthFullPairwiseAlignment) return null;
						String subjectStr = subject.subSequence(subjectNext,kmerHit.getStart()).toString();
						String queryStr = query.subSequence(queryNext,kmerHit.getQueryIdx()).toString();
						//if (subjectNextLength>10 || queryNextLength>10) System.out.println("Aligning segment of length "+subjectNextLength+" of subject with total length: "+subject.length()+" to segment with length "+queryNextLength+" of query with total length: "+query.length());
						String [] alignedFragments = aligner.getAlignment(subjectStr, queryStr);
						String cigarSegment = buildCigar(alignedFragments[0],alignedFragments[1]);
						cigar.append(cigarSegment);
						
					} else if (subjectNextLength>0) {
						cigar.append(""+subjectNextLength+""+deletionOp);
					} else if (queryNextLength>0) {
						cigar.append(""+queryNextLength+""+insertionOp);
					}
				}
				nextMatchLength+=kmerLength;
				subjectNext = kmerHit.getStart()+kmerLength;
				queryNext = kmerHit.getQueryIdx()+kmerLength;
			} else {
				int kmerSubjectNext = kmerHit.getStart()+kmerLength;
				int diffSubject = kmerSubjectNext - subjectNext;
				int kmerQueryNext = kmerHit.getQueryIdx()-queryNext;
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
			cigar.append(""+nextMatchLength+""+matchOp);
			nextMatchLength = 0;
		}
		//int alnFirst = alnStart+1;
		int alnLast = subjectNext;
		//System.out.println("Aligned query. first: "+alnStart+" last: "+alnLast+" CIGAR: "+cigar+" query next: "+queryNext+" query length: "+query.length());
		if(queryNext<query.length()) {
			//TODO: check if it is worth to align the sequence end
			/*
			if(subjectNext<subject.length()) {
				String queryStr = query.substring(queryNext);
				String subjectStr = subject.substring(subjectNext,Math.min(clusterLast, subject.length()));
				//TODO: Softclip if it starts with indel
				System.out.println("Aligning end of length "+subjectStr.length()+" of subject subsequence with total length: "+subject.length()+" to end with length "+queryStr.length()+" of query with total length: "+query.length());
				String [] alignedFragments = aligner.getAlignment(subjectStr, queryStr);
				String cigarSegment = buildCigar(alignedFragments[0],alignedFragments[1]);
				
				cigar.append(cigarSegment);	
			} else {*/
				//Ignore last bp
			int remainder = query.length()-queryNext;
			cigar.append(""+remainder+""+softClipOp);
			//}
		}
		ReadAlignment finalAlignment = new ReadAlignment(subjectIdx, alnStart+1, alnLast, query.length(), 0);
		finalAlignment.setReadCharacters(query);
		finalAlignment.setCigarString(cigar.toString());
		//TODO: Define better alignment quality
		double cov = 1.0*(queryNext-queryStart)/query.length();
		finalAlignment.setAlignmentQuality((byte) Math.round(100*cov));
		return finalAlignment;
	}

	private String buildCigar(String subjectAln, String queryAln) {
		StringBuilder cigar = new StringBuilder();
		char nextOperator = 0;
		int nextLength = 0;
		for(int i=0;i<subjectAln.length();i++) {
			char subjectChar = subjectAln.charAt(i);
			char queryChar = queryAln.charAt(i);
			char op = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
			if(subjectChar == AlignmentConstantGap.GAP_CHARACTER) {
				op = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION);
			} else if(queryChar == AlignmentConstantGap.GAP_CHARACTER) {
				op = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_DELETION);
			}
			if(op != nextOperator) {
				if(nextLength>0) {
					cigar.append(""+nextLength+""+nextOperator);
				}
				nextOperator = op;
				nextLength = 0;
			}
			nextLength++;
		}
		if (nextLength>0) cigar.append(""+nextLength+""+nextOperator);
		return cigar.toString();
	}
}
