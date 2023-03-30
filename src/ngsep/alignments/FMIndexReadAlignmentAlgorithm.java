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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.genome.GenomicRegion;
import ngsep.genome.ReferenceGenomeFMIndex;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.RawRead;
import ngsep.sequences.UngappedSearchHit;

/**
 *
 * @author German Andrade
 * @author Jorge Duitama
 *
 */
public class FMIndexReadAlignmentAlgorithm implements ReadAlignmentAlgorithm {
	
	private int kmerLength;
	private ReferenceGenomeFMIndex fMIndex;
	private ShortReadsUngappedSearchHitsClusterAligner aligner = new ShortReadsUngappedSearchHitsClusterAligner();
	
	private Set<CharSequence> repetitiveKmers = new HashSet<CharSequence>();
	
	private boolean onlyPositiveStrand = false;
	

	public FMIndexReadAlignmentAlgorithm(ReferenceGenomeFMIndex fMIndex, int kmerLength) {
		this.fMIndex = fMIndex;
		this.kmerLength = kmerLength;
	}
	public ReferenceGenomeFMIndex getFMIndex() {
		return fMIndex;
	}
	public void setFMIndex(ReferenceGenomeFMIndex fMIndex) {
		this.fMIndex = fMIndex;
	}

	public ShortReadsUngappedSearchHitsClusterAligner getAligner() {
		return aligner;
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
		if(readSeq.length()<500) {
			int maxMismatches = 2;
			alignments.addAll(fewMismatchesSingleStrandSearch(readSeq,maxMismatches));
			//System.out.println("Read: "+read.getName()+" Forward exact alignments: "+alignments.size());
			if(reverseComplement!=null) {
				List<ReadAlignment> alnsR = fewMismatchesSingleStrandSearch(reverseComplement,maxMismatches);
				//System.out.println("Read: "+read.getName()+" Reverse exact alignments: "+alnsR.size());
				for (ReadAlignment aln:alnsR) aln.setNegativeStrand(true);
				alignments.addAll(alnsR);
			}
		}
		if(alignments.size()==0) {
			alignments.addAll(alignQueryToReference(readSeq));
			//System.out.println("Read: "+read.getName()+" Forward inexact alignments: "+alignments);
			if(reverseComplement!=null) {
				List<ReadAlignment> alnsR = alignQueryToReference(reverseComplement);
				//System.out.println("Read: "+read.getName()+" Reverse inexact alignments: "+alnsR);
				for (ReadAlignment aln:alnsR) aln.setNegativeStrand(true);
				alignments.addAll(alnsR);
			}
		}
		
		//System.out.println("Read: "+read.getName()+" total alignments: "+alignments.size());
		for(ReadAlignment aln:alignments) {
			aln.setReadName(read.getName());
			if(!aln.isNegativeStrand()) aln.setQualityScores(qual);
			else aln.setQualityScores(reverseQS);
		}
		return alignments;
	}
	public List<ReadAlignment> fewMismatchesSingleStrandSearch(String query, int maxMismatches) {
		List<ReadAlignment> alns = new ArrayList<ReadAlignment>();
		List<Integer> alignment = new ArrayList<Integer>(1);
		alignment.add(ReadAlignment.getAlnValue(query.length(), ReadAlignment.ALIGNMENT_MATCH));
		// Whole read exact search
		List<UngappedSearchHit> readHits=fMIndex.exactSearch(query);
		for(int i=0;i<readHits.size();i++) {
			UngappedSearchHit hit = readHits.get(i);
			String seqName = fMIndex.getReferenceName(hit.getSubjectIdx());
			CharSequence subject = fMIndex.getSequence(seqName);
			ReadAlignment aln = aligner.buildAln (query, hit.getSubjectIdx(), subject, hit.getSubjectStart()+1, hit.getSubjectStart()+query.length(), alignment);
			if(aln==null) continue;
			aln.setSequenceName(seqName);
			aln.setAlignmentQuality((byte) 100);
			alns.add(aln);
		}
		if (alns.size()>0) return alns;
		//One mismatch search
		int middle = query.length()/2;
		if(middle < 50) return alns;
		String firstPart = query.substring(0,middle);
		readHits=fMIndex.exactSearch(firstPart);
		for(UngappedSearchHit hit: readHits) {
			String seqName = fMIndex.getReferenceName(hit.getSubjectIdx());
			CharSequence subject = fMIndex.getSequence(seqName);
			ReadAlignment aln = aligner.buildAln (query, hit.getSubjectIdx(), subject,  hit.getSubjectStart()+1, hit.getSubjectStart()+query.length(), alignment);
			if (aln==null) continue;
			int[] mismatches = aligner.countMismatches(query, subject, aln);
			if(mismatches==null) continue;
			//System.out.println("first half. Next aln: "+aln.getSequenceName()+":"+aln.getFirst()+" mismatches: "+mismatches[0]+" CIGAR: "+aln.getCigarString()+" clip: "+mismatches[1]+" "+mismatches[2]);
			if(mismatches[0]<=maxMismatches) {
				if (mismatches[1]+mismatches[2]>0) {
					aln = aligner.buildAln(query, hit.getSubjectIdx(), subject, hit.getSubjectStart()+1+mismatches[1], hit.getSubjectStart()+query.length()-mismatches[2], aligner.encodeAlignment(query.length(),mismatches));
				}
				if(aln==null) continue;
				aln.setAlignmentQuality((byte) (100-5*mismatches[0]));
				aln.setNumMismatches((short) mismatches[0]);
				
				alns.add(aln);
			}
		}
		if (alns.size()>0) return alns;
		String secondPart = query.substring(middle);
		readHits=fMIndex.exactSearch(secondPart);
		for(UngappedSearchHit hit: readHits) {
			String seqName = fMIndex.getReferenceName(hit.getSubjectIdx());
			CharSequence subject = fMIndex.getSequence(seqName);
			int start = hit.getSubjectStart()-middle;
			if(start<0) continue;
			ReadAlignment aln = aligner.buildAln (query, hit.getSubjectIdx(), "", start+1, start+query.length(), alignment);
			if (aln==null) continue;
			int[] mismatches = aligner.countMismatches(query, subject, aln);
			if(mismatches==null) continue;
			//System.out.println("Second half. Next aln: "+aln.getSequenceName()+":"+aln.getFirst()+" mismatches: "+mismatches[0]+" CIGAR: "+aln.getCigarString()+" clip: "+mismatches[1]+" "+mismatches[2]);
			if(mismatches[0]<=maxMismatches) {
				if (mismatches[1]+mismatches[2]>0) {
					aln = aligner.buildAln(query, hit.getSubjectIdx(), subject, start+1+mismatches[1], start+query.length()-mismatches[2], aligner.encodeAlignment(query.length(),mismatches));
				}
				if(aln==null) continue;
				aln.setAlignmentQuality((byte) (100-5*mismatches[0]));
				aln.setNumMismatches((short) mismatches[0]);
				alns.add(aln);
			}
		}
		return alns;
	}
	
	
	public List<ReadAlignment> alignQueryToReference (String query) {
		
		return kmerBasedSingleStrandInexactSearchAlgorithm(query);
	}
	
	/**
	 * Inexact search of kmers to an FM-index
	 * It iterates the Sequences of the genome and if there is at least MIN_ACCURACY percentage of the kmers
	 * it allow the alignment with the first an the last position of the kmers ocurrence
	 * Only tries to align the given quey in the positive strand
	 * @return List<ReadAlignment>
	 */
	private List<ReadAlignment> kmerBasedSingleStrandInexactSearchAlgorithm (String query) 
	{
		Map<Integer,String> kmersMap = KmersExtractor.extractKmersAsMap(query, kmerLength, 15, true, false, false);
		List<ReadAlignment> finalAlignments =  new ArrayList<>();
		//System.out.println("Read: "+query+" length "+query.length()+" kmers: "+kmersMap.size());
		int kmersCount=kmersMap.size();
		if(kmersCount==0) return finalAlignments;
		List<UngappedSearchHit> initialKmerHits = searchKmers (kmersMap);
		List<UngappedSearchHitsCluster> clusteredKmerHits = clusterKmerHits(query, initialKmerHits);
		if(clusteredKmerHits.size()==0) return finalAlignments;
		//System.out.println("Initial kmer hits: "+initialKmerHits.size()+" Clusters: "+clusteredKmerHits.size());
		Collections.sort(clusteredKmerHits, (o1, o2) -> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
		//KmerHitsCluster cluster = clusteredKmerHits.get(0);
		//ReadAlignment readAln = createNewAlignmentFromConsistentKmers(cluster, query);
		//if(readAln!=null) finalAlignments.add(readAln);
		int kmersMaxCluster = 0;
		for (int i=0;i<clusteredKmerHits.size();i++) {
			UngappedSearchHitsCluster cluster = clusteredKmerHits.get(i);
			int numKmers = cluster.getNumDifferentKmers();
			//System.out.println("Processing cluster "+i+" spanning "+cluster.getSubjectIdx()+":"+cluster.getSubjectPredictedStart()+"-"+cluster.getSubjectPredictedEnd()+" Num kmers: "+cluster.getNumDifferentKmers()+" consistent: "+cluster.isAllConsistent()+" max kmers: "+kmersMaxCluster);
			if(i==0) kmersMaxCluster = numKmers;
			else if (finalAlignments.size()>0 && (numKmers<2 || numKmers< 0.5*kmersMaxCluster)) break;
			ReadAlignment readAln = aligner.buildAlignment(query, fMIndex.getSequence(fMIndex.getReferenceName(cluster.getSubjectIdx())), cluster);
			if(readAln!=null) finalAlignments.add(readAln);
		}
		//System.out.println("Found "+finalAlignments.size()+" alignments for query: "+query);
		return finalAlignments;
	}


	/**
	 * Searches the given kmers in the fmIndex 
	 * @param kmers to search
	 * @return List of alignments of each kmer. The read number of each alignment contains the kmer number.
	 */
	private List<UngappedSearchHit> searchKmers(Map<Integer,String> kmersMap) {
		List<UngappedSearchHit> answer = new ArrayList<>();
		for (int start:kmersMap.keySet()) {
			String kmer = kmersMap.get(start);
			CharSequence kmerP = KmersExtractor.pack(kmer);
			if(repetitiveKmers.contains(kmerP)) continue;
			List<UngappedSearchHit> kmerHits=fMIndex.exactSearch(kmer);
			//System.out.println("Kmer: "+kmer+" hits: "+kmerHits.size());
			if(kmerHits.size()>=ReferenceGenomeFMIndex.MAX_HITS_QUERY) {
				repetitiveKmers.add(kmerP);
				continue;
			}
			for(UngappedSearchHit hit:kmerHits) {
				//System.out.println("Next hit: "+hit.getSequenceName()+" "+hit.getStart());
				hit.setQueryStart(start);
				answer.add(hit);
			}
		}
		return answer;
	}

	private List<UngappedSearchHitsCluster> clusterKmerHits(String query, List<UngappedSearchHit> initialKmerHits) {
		List<UngappedSearchHitsCluster> clusters = new ArrayList<>();
		Map<Integer,List<UngappedSearchHit>> hitsBySubjectIdx = new LinkedHashMap<Integer, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:initialKmerHits) {
			List<UngappedSearchHit> hitsSeq = hitsBySubjectIdx.computeIfAbsent(hit.getSubjectIdx(), k -> new ArrayList<>());
			hitsSeq.add(hit);
		}
		for(int subjectIdx:hitsBySubjectIdx.keySet()) {
			int subjectLength = fMIndex.getReferenceLength(subjectIdx);
			List<UngappedSearchHit> hitsSeq = hitsBySubjectIdx.get(subjectIdx);
			Collections.sort(hitsSeq, (hit0,hit1)-> hit0.getSubjectStart()-hit1.getSubjectStart());
			clusters.addAll(clusterSequenceKmerAlns(query, subjectIdx, subjectLength, hitsSeq));
		}
		return clusters;
	}
	
	private List<UngappedSearchHitsCluster> clusterSequenceKmerAlns(String query, int subjectIdx, int sequenceLength, List<UngappedSearchHit> sequenceHits) {
		List<UngappedSearchHitsCluster> answer = new ArrayList<>();
		//System.out.println("Alns to cluster: "+sequenceAlns.size());
		UngappedSearchHitsCluster cluster=null;
		for(UngappedSearchHit kmerHit:sequenceHits) {
			if(cluster==null || !cluster.addKmerHit(kmerHit, 0)) {
				cluster = new UngappedSearchHitsCluster(query.length(), subjectIdx, sequenceLength, kmerHit);
				answer.add(cluster);
			}
		}
		return answer;
	}
	public void loadSTRsFile(String knownSTRsFile) throws IOException {
		aligner.loadSTRsFile(knownSTRsFile);
	}
	public Map<String, List<GenomicRegion>> getKnownSTRs() {
		return aligner.getKnownSTRs();
	}
	public void setKnownSTRs(Map<String, List<GenomicRegion>> knownSTRs) {
		aligner.setKnownSTRs(knownSTRs);
	}
	
	

	
}
