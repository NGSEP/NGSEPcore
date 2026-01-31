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
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.genome.ReferenceGenomeFMIndex;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.UngappedSearchHit;

/**
 *
 * @author German Andrade
 * @author Jorge Duitama
 *
 */
public class FMIndexUngappedSearchHitsClustersFinder implements UngappedSearchHitsClustersFinder {
	
	private int kmerLength;
	private ReferenceGenomeFMIndex fMIndex;
	
	private Set<CharSequence> repetitiveKmers = new HashSet<CharSequence>();
	

	public FMIndexUngappedSearchHitsClustersFinder(ReferenceGenomeFMIndex fMIndex, int kmerLength) {
		this.fMIndex = fMIndex;
		this.kmerLength = kmerLength;
	}
	public ReferenceGenomeFMIndex getFMIndex() {
		return fMIndex;
	}
	public void setFMIndex(ReferenceGenomeFMIndex fMIndex) {
		this.fMIndex = fMIndex;
	}
	
	@Override
	public List<UngappedSearchHitsCluster> findHitClusters(CharSequence query) {
		Map<Integer,String> kmersMap = KmersExtractor.extractKmersAsMap(query.toString(), kmerLength, 15, true, false, false);
		int kmersCount=kmersMap.size();
		if(kmersCount==0) return new ArrayList<>();
		List<UngappedSearchHit> initialKmerHits = searchKmers (kmersMap);
		return clusterKmerHits(query, initialKmerHits);
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
				//System.out.println("Next hit: "+hit.getSubjectIdx()+" "+hit.getSubjectStart());
				hit.setQueryStart(start);
				answer.add(hit);
			}
		}
		return answer;
	}

	private List<UngappedSearchHitsCluster> clusterKmerHits(CharSequence query, List<UngappedSearchHit> initialKmerHits) {
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
	
	private List<UngappedSearchHitsCluster> clusterSequenceKmerAlns(CharSequence query, int subjectIdx, int sequenceLength, List<UngappedSearchHit> sequenceHits) {
		List<UngappedSearchHitsCluster> answer = new ArrayList<>();
		//System.out.println("Alns to cluster: "+sequenceAlns.size());
		UngappedSearchHitsCluster cluster=null;
		for(UngappedSearchHit kmerHit:sequenceHits) {
			if(cluster==null || !cluster.addKmerHit(kmerHit, query.length()/2)) {
				cluster = new UngappedSearchHitsCluster(query.length(), subjectIdx, sequenceLength, kmerHit);
				answer.add(cluster);
			}
		}
		return answer;
	}
	public List<ReadAlignment> fewMismatchesSingleStrandSearch(String query, int maxMismatches) {
		ShortReadsUngappedSearchHitsClusterAligner aligner = new ShortReadsUngappedSearchHitsClusterAligner();
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
}
