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
import java.util.List;
import java.util.Map;

import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.RawRead;

public class PairedReadsAligner {
	private ReferenceGenome genome;
	private UngappedSearchHitsClustersFinder clustersFinder;
	private UngappedSearchHitsClusterAligner aligner;
	private int minInsertLength=0;
	private int maxInsertLength=1000;
	private int maxAlnsPerRead=1;
	private boolean createUnmappedReadRecords = false;
	//Alignment statistics
	private int numReadsAligned = 0;
	private boolean properPair = false;
	private boolean pair = false;
	private int numUniqueAlignments = 0;
	
	public PairedReadsAligner(ReferenceGenome genome, UngappedSearchHitsClustersFinder clustersFinder, UngappedSearchHitsClusterAligner aligner) {
		super();
		this.genome = genome;
		this.clustersFinder = clustersFinder;
		this.aligner = aligner;
	}
	
	

	public int getMaxAlnsPerRead() {
		return maxAlnsPerRead;
	}
	public void setMaxAlnsPerRead(int maxAlnsPerRead) {
		this.maxAlnsPerRead = maxAlnsPerRead;
	}

	public boolean isCreateUnmappedReadRecords() {
		return createUnmappedReadRecords;
	}
	public void setCreateUnmappedReadRecords(boolean createUnmappedReadRecords) {
		this.createUnmappedReadRecords = createUnmappedReadRecords;
	}
	
	public int getMinInsertLength() {
		return minInsertLength;
	}
	public void setMinInsertLength(int minInsertLength) {
		this.minInsertLength = minInsertLength;
	}

	public int getMaxInsertLength() {
		return maxInsertLength;
	}
	public void setMaxInsertLength(int maxInsertLength) {
		this.maxInsertLength = maxInsertLength;
	}



	public int getNumReadsAligned() {
		return numReadsAligned;
	}

	public boolean isProperPair() {
		return properPair;
	}

	public boolean isPair() {
		return pair;
	}

	public int getNumUniqueAlignments() {
		return numUniqueAlignments;
	}
	
	
	public List<ReadAlignment> alignReads(RawRead read1, RawRead read2) {
		CharSequence query1F = read1.getCharacters();
		CharSequence query1R = DNAMaskedSequence.getReverseComplement(query1F);
		CharSequence query2F = read2.getCharacters();
		CharSequence query2R = DNAMaskedSequence.getReverseComplement(query2F);
		List<UngappedSearchHitsCluster> clusters1F = clustersFinder.findHitClusters(query1F);
		List<UngappedSearchHitsCluster> clusters1R = clustersFinder.findHitClusters(query1R);
		List<UngappedSearchHitsCluster> clusters2F = clustersFinder.findHitClusters(query2F);
		//System.out.println("Searching reverse of: "+read2.getName());
		List<UngappedSearchHitsCluster> clusters2R = clustersFinder.findHitClusters(query2R);
		ReadAlignment unmapped1 = ReadAlignment.createMockAlignmentUnmappedRead(read1, true, true);
		ReadAlignment unmapped2 = ReadAlignment.createMockAlignmentUnmappedRead(read2, true, false); 
		//System.out.println("Clusters found: "+clusters1F.size()+" "+clusters1R.size()+" "+clusters2F.size()+" "+clusters2R.size());
		List<ReadAlignment> alns = new ArrayList<ReadAlignment>();
		if(clusters1F.size()+clusters1R.size()==0) {
			if (createUnmappedReadRecords) alns.add(unmapped1);
			List<ReadAlignment> alns2 = alignAsSingle(query2F, clusters2F, query2R, clusters2R);
			if(alns2.size()==0) {
				unmapped1.setMateUnmapped(true);
				unmapped2.setMateUnmapped(true);
				if (createUnmappedReadRecords) alns.add(unmapped2);
			} else {
				int n =Math.min(alns2.size(), maxAlnsPerRead);
				for(int i=0;i<n;i++) {
					ReadAlignment aln2 = alns2.get(i);
					aln2.setPaired(true);
					aln2.setSecondOfPair(true);
					aln2.setMateUnmapped(true);
					setMateInfo(aln2, unmapped1);
					if(i>0) aln2.setSecondary(true);
					else setMateInfo(unmapped1, aln2);
					if (aln2.getAlignmentQuality()>=10) aln2.setAlignmentQuality((byte) Math.round(0.5*aln2.getAlignmentQuality()));
					alns.add(aln2);
				}
				if(alns2.size()==1 && alns2.get(0).getAlignmentQuality()>=20) numUniqueAlignments++;
				numReadsAligned = 1;
			}
		} else if(clusters2F.size()+clusters2R.size()==0) {
			if (createUnmappedReadRecords) alns.add(unmapped2);
			List<ReadAlignment> alns1 = alignAsSingle(query1F, clusters1F, query1R, clusters1R);
			if(alns1.size()==0) {
				unmapped1.setMateUnmapped(true);
				unmapped2.setMateUnmapped(true);
				if (createUnmappedReadRecords) alns.add(unmapped1);
			} else {
				int n =Math.min(alns1.size(), maxAlnsPerRead);
				for(int i=0;i<n;i++) {
					ReadAlignment aln1 = alns1.get(i);
					aln1.setPaired(true);
					aln1.setFirstOfPair(true);
					aln1.setMateUnmapped(true);
					if(i>0) aln1.setSecondary(true);
					else setMateInfo(unmapped2, aln1);
					if (aln1.getAlignmentQuality()>=10) aln1.setAlignmentQuality((byte) Math.round(0.5*aln1.getAlignmentQuality()));
					alns.add(aln1);
				}
				if(alns1.size()==1 && alns1.get(0).getAlignmentQuality()>=20) numUniqueAlignments++;
				
				
				if(alns1.size()>0) numReadsAligned = 1;
			}
		} else {
			List<ReadAlignmentPair> pairedAlns = new ArrayList<>();
			pairedAlns.addAll( findPairs(query1F, clusters1F, false, query2R, clusters2R, true, true));
			pairedAlns.addAll( findPairs(query1R, clusters1R, true, query2F, clusters2F, false, true));
			//System.out.println("Pairs proper: "+pairedAlns.size());
			properPair=pairedAlns.size()>0;
			if(pairedAlns.isEmpty() ) {
				pairedAlns.addAll(findPairs(query1F, clusters1F, false, query2R, clusters2R, true, false));
				pairedAlns.addAll(findPairs(query1R, clusters1R, true, query2F, clusters2F, false, false));
			}
			//System.out.println("Non proper 1: "+pairedAlns.size());
			if(pairedAlns.isEmpty()) {
				pairedAlns.addAll(findPairs(query1F, clusters1F, false, query2F, clusters2F, false, false));
				pairedAlns.addAll(findPairs(query1R, clusters1R, true, query2R, clusters2R, true, false));	
			}
			//System.out.println("Non proper 2: "+pairedAlns.size());
			if(pairedAlns.isEmpty()) {
				List<ReadAlignment> alns1 = alignAsSingle(query1F, clusters1F, query1R, clusters1R);
				List<ReadAlignment> alns2 = alignAsSingle(query2F, clusters2F, query2R, clusters2R);
				//System.out.println("Alignments 1: "+alns1);
				//System.out.println("Alignments 2: "+alns2);
				if(alns1.size()>0) numReadsAligned++;
				if(alns2.size()>0) numReadsAligned++;
				int n = Math.min(alns1.size(), maxAlnsPerRead);
				for(int i=0;i<n;i++) {
					ReadAlignment aln1 = alns1.get(i);
					aln1.setPaired(true);
					aln1.setFirstOfPair(true);
					if(alns2.size()>0) {
						aln1.setMateDifferentSequence(true);
						setMateInfo(aln1, alns2.get(0));
					} else {
						aln1.setMateUnmapped(true);
						setMateInfo(aln1, unmapped2);
						if (createUnmappedReadRecords && i==0) {
							setMateInfo(unmapped2, aln1);
							alns.add(unmapped2);
						}
					}
					if(i>0) aln1.setSecondary(true);
					if (aln1.getAlignmentQuality()>=10) aln1.setAlignmentQuality((byte) Math.round(0.5*aln1.getAlignmentQuality()));
					alns.add(aln1);
				}
				if(alns1.size()==1 && alns1.get(0).getAlignmentQuality()>=20) numUniqueAlignments++;
				n =Math.min(alns2.size(), maxAlnsPerRead);
				for(int i=0;i<n;i++) {
					ReadAlignment aln2 = alns2.get(i);
					aln2.setPaired(true);
					aln2.setSecondOfPair(true);
					
					if(alns1.size()>0) {
						aln2.setMateDifferentSequence(true);
						setMateInfo(aln2, alns1.get(0));
					} else {
						aln2.setMateUnmapped(true);
						setMateInfo(aln2, unmapped1);
						if (createUnmappedReadRecords && i==0) {
							setMateInfo(unmapped1, aln2);
							alns.add(unmapped1);
						}
					}
					if(i>0) aln2.setSecondary(true);
					if (aln2.getAlignmentQuality()>=10) aln2.setAlignmentQuality((byte) Math.round(0.5*aln2.getAlignmentQuality()));
					alns.add(aln2);
				}
				if(alns2.size()==1 && alns2.get(0).getAlignmentQuality()>=20) numUniqueAlignments++;
			} else {
				numReadsAligned = 2;
				addPairedAlignments(alns, pairedAlns, clusters1F.size()+clusters1R.size(),clusters2F.size()+ clusters2R.size());
				pair=true;
				if(pairedAlns.size()==1) numUniqueAlignments=2;
			}
		}
		for(ReadAlignment aln:alns) {
			String qs;
			int readLength;
			if(aln.isFirstOfPair()) {
				aln.setReadName(read1.getName());
				qs = read1.getQualityScores();
				readLength = read1.getLength();
			}
			else {
				aln.setReadName(read2.getName());
				qs = read2.getQualityScores();
				readLength = read2.getLength();
			}
			if(qs!=null && qs.length()!=readLength) {
				System.err.println("WARNING. Inconsistent length of quality scores for read. "+aln.getReadName()+". Read length: "+readLength+" scores length: "+qs.length()+". Ignoring scores");
				qs = null;
			}
			if(qs!=null && aln.isNegativeStrand()) qs = new StringBuilder(qs).reverse().toString();
			aln.setQualityScores(qs);
		}
		return alns;
	}

	private List<ReadAlignmentPair> findPairs(CharSequence query1, List<UngappedSearchHitsCluster> clusters1, boolean r1, CharSequence query2, List<UngappedSearchHitsCluster> clusters2, boolean r2, boolean onlyProper ){
		List<UngappedSearchHitClusterPair> pairedClusters = new ArrayList<>();
		Map<Integer,List<UngappedSearchHitsCluster>> clustersBySubjectIdx1 = organizeBySubjectIdx(clusters1);
		Map<Integer,List<UngappedSearchHitsCluster>> clustersBySubjectIdx2 = organizeBySubjectIdx(clusters2);
		for(Map.Entry<Integer, List<UngappedSearchHitsCluster>> entry:clustersBySubjectIdx1.entrySet()) {
			int subjectIdx = entry.getKey();
			List<UngappedSearchHitsCluster> clustersSubject1 = entry.getValue();
			List<UngappedSearchHitsCluster> clustersSubject2 = clustersBySubjectIdx2.get(subjectIdx);
			if(clustersSubject2==null) continue;
			pairedClusters.addAll(findPairsSubject(clustersSubject1, r1, clustersSubject2, r2, onlyProper ));
		}
		//System.out.println("Paired clusters: "+pairedClusters.size()+" proper: "+onlyProper);
		return buildPairedAlignments(query1, query2, pairedClusters, aligner);
	}

	private Map<Integer, List<UngappedSearchHitsCluster>> organizeBySubjectIdx(List<UngappedSearchHitsCluster> clusters) {
		Map<Integer, List<UngappedSearchHitsCluster>> answer = new HashMap<>();
		for(UngappedSearchHitsCluster cluster:clusters) {
			List<UngappedSearchHitsCluster> clustersSeq = answer.computeIfAbsent(cluster.getSubjectIdx(), v->new ArrayList<>());
			clustersSeq.add(cluster);
		}
		return answer;
	}

	private List<UngappedSearchHitClusterPair> findPairsSubject(List<UngappedSearchHitsCluster> clusters1, boolean r1, List<UngappedSearchHitsCluster> clusters2, boolean r2, boolean onlyProper) {
		List<UngappedSearchHitClusterPair> pairs = new ArrayList<>();
		int n1 = clusters1.size();
		int n2 = clusters2.size();
		boolean [] clustersPaired = new boolean[n2];
		int initialIndex2 = 0;
		for (int i = 0; i < n1; i++) {
			UngappedSearchHitsCluster c1 = clusters1.get(i);
			if(onlyProper) {
				//Update initial index
				while(initialIndex2<n2) {
					UngappedSearchHitsCluster c2 = clusters2.get(initialIndex2);
					int distance = c1.getSubjectPredictedEnd()-c2.getSubjectPredictedStart();
					if(distance <=maxInsertLength) break;
					initialIndex2++;
				}
				if(initialIndex2==n2) break;
			}
			UngappedSearchHitClusterPair pair = findPair(c1,r1,clusters2,initialIndex2,r2,clustersPaired, onlyProper);
			if(pair!=null) pairs.add(pair);
		}
		return pairs;
	}

	private UngappedSearchHitClusterPair findPair(UngappedSearchHitsCluster c1, boolean r1, List<UngappedSearchHitsCluster> clusters2,int initialIndex2, boolean r2, boolean [] clustersPaired, boolean onlyProper) {
		List<Integer> candidatePositions = new ArrayList<>();
		int n = clusters2.size();
		for (int i = initialIndex2; i < n; i++) {
			UngappedSearchHitsCluster current =clusters2.get(i);
			//System.out.println("Next candidate "+current+" is paired: "+clustersPaired[i]);
			if(!clustersPaired[i]) {
				if(!onlyProper || isValidPair(c1,r1,current,r2)) {
					candidatePositions.add(i);
				}
			}
			if(onlyProper) {
				int distance = current.getSubjectPredictedEnd()-c1.getSubjectPredictedStart();
				if(distance>maxInsertLength) break;
			}
		}
		//System.out.println("Candidate clusters for pairing: "+candidatePositions);
		if(candidatePositions.size()==0) return null;
		int pos =  pickBestPairPosition(c1, clusters2, candidatePositions);
		clustersPaired[pos]=true;
		//System.out.println("Adding candidate. Aln 1: "+aln1+" Aln 2: "+current);
		return new UngappedSearchHitClusterPair(c1, r1, clusters2.get(pos),r2, onlyProper);
	}
	private boolean isValidPair(UngappedSearchHitsCluster c1, boolean r1, UngappedSearchHitsCluster c2,boolean r2) {
		if(c1.getSubjectIdx()!=c2.getSubjectIdx()) return false;
		int start1 = c1.getSubjectPredictedStart();
		int start2 = c2.getSubjectPredictedStart();
		int end1 = c1.getSubjectPredictedEnd();
		int end2 = c2.getSubjectPredictedEnd();
		boolean properDirection;
		int insertLength;
		if(start1 < end2) {
			insertLength = end2-start1+1;
			properDirection = !r1 && r2;
		} else {
			insertLength = end1-start2+1;
			properDirection = !r2 && r1;
		}
		//System.out.println("Min insert length: "+minInsertLength+" max insert length: "+maxInsertLength+" next: "+insertLength);
		return properDirection && insertLength>=minInsertLength && insertLength<=maxInsertLength;
	}
	private int pickBestPairPosition(UngappedSearchHitsCluster c1, List<UngappedSearchHitsCluster> clustersSubject2, List<Integer> candidatePositions) {
		int start1 = c1.getSubjectPredictedStart();
		int end1 = c1.getSubjectPredictedEnd();
		int answer = -1;
		int minDistance = -1;
		int middle = (maxInsertLength+minInsertLength)/2;
		for (int i:candidatePositions) {
			UngappedSearchHitsCluster c2 = clustersSubject2.get(i);
			int start2 = c2.getSubjectPredictedStart();
			int end2 = c2.getSubjectPredictedEnd();
			int insertLength;
			if(start1 < end2) {
				insertLength = end2-start1+1;
			} else {
				insertLength = end1-start2+1;
			}
			int d = Math.abs(middle - insertLength);
			if(minDistance==-1 || minDistance>d) {
				answer = i;
				minDistance = d;
			}
		}
		return answer;
	}
	
	private List<ReadAlignmentPair> buildPairedAlignments(CharSequence query1, CharSequence query2, List<UngappedSearchHitClusterPair> pairedClusters, UngappedSearchHitsClusterAligner aligner) {
		List<ReadAlignmentPair> answer = new ArrayList<>();
		Collections.sort(pairedClusters,(c1,c2)->c2.getScore()-c1.getScore());
		for(UngappedSearchHitClusterPair pair:pairedClusters) {
			CharSequence subject = genome.getSequenceCharacters(pair.getCluster1().getSubjectIdx());
			ReadAlignment aln1 = aligner.buildAlignment(query1, subject, pair.getCluster1());
			//System.out.println("Aln1: "+aln1);
			if(aln1==null) continue;
			String referenceName = genome.getSequenceByIndex(aln1.getSequenceIndex()).getName();
			aln1.setSequenceName(referenceName);
			aln1.setNegativeStrand(pair.isReverse1());
			ReadAlignment aln2 = aligner.buildAlignment(query2, subject, pair.getCluster2());
			//System.out.println("Aln2: "+aln2+" cluster predicted limits: "+pair.getCluster2().getSubjectPredictedStart()+" " +pair.getCluster2().getSubjectPredictedEnd());
			if(aln2==null) continue;
			aln2.setSequenceName(referenceName);
			aln2.setNegativeStrand(pair.isReverse2());
			answer.add(buildPair(aln1, aln2, pair.isProper()));
			if(answer.size()>3*maxAlnsPerRead) break;
		}
		return answer;
	}

	

	private ReadAlignmentPair buildPair(ReadAlignment aln1, ReadAlignment aln2, boolean proper) {
		aln1.setPaired(true);
		aln1.setFirstOfPair(true);
		aln1.setProperPair(proper);
		aln2.setPaired(true);
		aln2.setSecondOfPair(true);
		aln2.setProperPair(proper);
		int insertLength1 = aln1.getLast()-aln2.getFirst()+1;
		int insertLength2 = aln2.getLast()-aln1.getFirst()+1;
		if(insertLength1>insertLength2) {
			aln1.setInferredInsertSize(-insertLength1);
			aln2.setInferredInsertSize(insertLength1);
		} else {
			aln1.setInferredInsertSize(insertLength2);
			aln2.setInferredInsertSize(-insertLength2);
		}
		setMateInfo(aln1,aln2);
		setMateInfo(aln2,aln1);
		return new ReadAlignmentPair(aln1,aln2);
	}
	
	private void addPairedAlignments(List<ReadAlignment> alns, List<ReadAlignmentPair> pairAlns, int numAlnsUnpaired1, int numAlnsUnpaired2) {
		Collections.sort(pairAlns,(p1,p2)->p2.getQualitySum()-p1.getQualitySum());
		int n = Math.min(pairAlns.size(),maxAlnsPerRead);
		for (int i = 0; i < n; i++) {
			ReadAlignmentPair current = pairAlns.get(i);
			ReadAlignment aln1 = current.getAln1();
			ReadAlignment aln2 = current.getAln2();
			if(i>0) {
				aln1.setSecondary(true);
				aln2.setSecondary(true);
			}
			if (pairAlns.size()>1) {
				aln1.setAlignmentQuality((byte) Math.round(0.2*aln1.getAlignmentQuality()/(double)n));
				aln2.setAlignmentQuality((byte) Math.round(0.2*aln2.getAlignmentQuality()/(double)n));
			} else if (!aln1.isProperPair() || !aln2.isProperPair() || (numAlnsUnpaired1>1 && numAlnsUnpaired2>1)) {
				double div = Math.max(numAlnsUnpaired1+numAlnsUnpaired2-1,1);
				aln1.setAlignmentQuality((byte) Math.round(0.5*aln1.getAlignmentQuality()/div));
				aln2.setAlignmentQuality((byte) Math.round(0.5*aln2.getAlignmentQuality()/div));
			}
			alns.add(aln1);
			alns.add(aln2);
		}
	}
	private ReadAlignment setMateInfo(ReadAlignment alignment,ReadAlignment mate) {
		alignment.setMateFirst(mate.getFirst());
		alignment.setMateSequenceName(mate.getSequenceName());
		alignment.setMateNegativeStrand(mate.isNegativeStrand());
		//System.out.println("FIlled mate info for alignments: "+alignment+" "+mate);
		return alignment;
	}
	
	private List<ReadAlignment> alignAsSingle(CharSequence queryF, List<UngappedSearchHitsCluster> clustersF, CharSequence queryR, List<UngappedSearchHitsCluster> clustersR ) {
		List<ReadAlignment> answer = new ArrayList<ReadAlignment>();
		answer.addAll(alignAsSingle(queryF, false, clustersF));
		answer.addAll(alignAsSingle(queryR, true, clustersR));
		Collections.sort(answer, (aln1,aln2) -> aln2.getAlignmentQuality() - aln1.getAlignmentQuality());
		return answer;
	}
	private List<ReadAlignment> alignAsSingle(CharSequence query, boolean negativeStrand, List<UngappedSearchHitsCluster> clusters) {
		SingleReadsAligner singleReadsAligner = new SingleReadsAligner(genome, clustersFinder, aligner);
		singleReadsAligner.setMaxAlnsPerRead(maxAlnsPerRead);
		List<ReadAlignment> rawAlns = singleReadsAligner.buildAlignments(query, clusters, 1);
		List<ReadAlignment> answer = singleReadsAligner.filterAlignments(rawAlns);
		for(ReadAlignment aln:answer) {
			aln.setNegativeStrand(negativeStrand);
			aln.setPaired(true);
		}
		return answer;
	}
}
class UngappedSearchHitClusterPair {
	private UngappedSearchHitsCluster cluster1;
	private boolean reverse1;
	private UngappedSearchHitsCluster cluster2;
	private boolean reverse2;
	private boolean proper= false;
	
	public UngappedSearchHitClusterPair(UngappedSearchHitsCluster cluster1,	boolean r1, UngappedSearchHitsCluster cluster2, boolean r2, boolean proper) {
		super();
		this.cluster1 = cluster1;
		reverse1 = r1;
		this.cluster2 = cluster2;
		reverse2 = r2;
		this.proper = proper;
	}
	

	public UngappedSearchHitsCluster getCluster1() {
		return cluster1;
	}

	public UngappedSearchHitsCluster getCluster2() {
		return cluster2;
	}
	

	public boolean isReverse1() {
		return reverse1;
	}


	public boolean isReverse2() {
		return reverse2;
	}
	
	public boolean isProper() {
		return proper;
	}


	public int getScore() {
		return (int) (10*cluster1.getCountKmerHitsCluster()+10*cluster2.getCountKmerHitsCluster());
	}
}
