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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import ngsep.graphs.DisjointSets;
import ngsep.math.Distribution;
import ngsep.sequences.UngappedSearchHit;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class UngappedSearchHitsClusterBuilder {
	private int idxSubjectDebug = -1;
	private int queryLengthDebug = -1;
	private boolean debug = false;
	public List<UngappedSearchHitsCluster> clusterRegionKmerAlns(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> sequenceHits, double minQueryCoverage) {
		Map<Integer,Integer> countsByQueryIdx = new HashMap<Integer, Integer>();
		double minHits = Math.min(20,0.01*queryLength);
		debug = subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug;
		if(debug) System.out.println("Clustering hits: "+sequenceHits.size()+" minHits: "+minHits);
		if(sequenceHits.size()<minHits) return new ArrayList<UngappedSearchHitsCluster>();
		for(UngappedSearchHit hit:sequenceHits) {
			//if (debug && hit.getSubjectStart()<10000) System.out.println("Next qpos "+hit.getQueryStart()+" hit: "+hit.getSubjectStart());
			countsByQueryIdx.compute(hit.getQueryStart(), (k,v)->v==null?1:v+1);
		}
		//double estimatedClusters = 0.5*sequenceHits.size()/queryLength;
		double avg = 0;
		for(int count:countsByQueryIdx.values()) avg+=count;
		if(avg>0) avg/=countsByQueryIdx.size();
		double estimatedClusters = avg;
		 
		if(debug) System.out.println("Clustering hits: "+sequenceHits.size()+" query starts: "+countsByQueryIdx.size()+" estimatedClusters: "+estimatedClusters);
		
		//if(minQueryCoverage==0) 
		//List<List<UngappedSearchHit>> hitsClusters = clusterRegionKmerAlnsMultiple(queryLength, subjectIdx, subjectLength, sequenceHits, estimatedClusters);
		List<List<UngappedSearchHit>> hitsClusters = clusterRegionKmerAlnsKruskal(queryLength, subjectIdx, subjectLength, sequenceHits, estimatedClusters);
		return buildClusterObjectsFromClusteredHits(queryLength, subjectIdx, subjectLength, hitsClusters, minHits, sequenceHits.size());
		//if(estimatedClusters>1.5) return clusterRegionKmerAlnsMultiple(queryLength, subjectLength, sequenceHits, estimatedClusters);
		/*List<UngappedSearchHitsCluster> answer = new ArrayList<>();
		List<UngappedSearchHit> selectedHits = collapseAndSelectSortedHits(queryLength, subjectLength, sequenceHits);
		UngappedSearchHitsCluster uniqueCluster = new UngappedSearchHitsCluster(queryLength, subjectLength, selectedHits);
		uniqueCluster.setRawKmerHits(sequenceHits.size());
		uniqueCluster.setRawKmerHitsSubjectStartSD(1);
		//if(querySequenceId==idxDebug) System.out.println("Hits to cluster: "+sequenceKmerHits.size()+" target: "+uniqueCluster.getSequenceIdx()+" first: "+uniqueCluster.getFirst()+" last: "+uniqueCluster.getLast()+" kmers: "+uniqueCluster.getNumDifferentKmers());
		if (uniqueCluster.getNumDifferentKmers()<minHits || uniqueCluster.getQueryEvidenceEnd()-uniqueCluster.getQueryEvidenceStart()<minQueryCoverage*queryLength) return answer;
		answer.add(uniqueCluster);
		if(uniqueCluster.getNumDifferentKmers()>0.8*sequenceHits.size()) return answer;
		//Cluster remaining hits
		List<UngappedSearchHit> remainingHits = new ArrayList<UngappedSearchHit>();
		for(UngappedSearchHit hit:sequenceHits) {
			if (hit!=uniqueCluster.getKmerHit(hit.getQueryIdx())) remainingHits.add(hit);
		}
		selectedHits = collapseAndSelectSortedHits(queryLength, subjectLength, remainingHits);
		UngappedSearchHitsCluster cluster2 = new UngappedSearchHitsCluster(queryLength, subjectLength, selectedHits);
		uniqueCluster.setRawKmerHits(remainingHits.size());
		uniqueCluster.setRawKmerHitsSubjectStartSD(1);
		if(cluster2.getNumDifferentKmers()>=minHits && cluster2.getQueryEvidenceEnd()-cluster2.getQueryEvidenceStart()>=minQueryCoverage*queryLength) answer.add(cluster2);
		return answer;*/
	}
	private List<List<UngappedSearchHit>> clusterRegionKmerAlnsMultiple(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> sequenceHits, double estimatedClusters) {
		
		
		//UngappedSearchHit firstHit = sequenceHits.get(0);
		//int subjectIdx = firstHit.getSequenceIdx();
		//Initial clustering
		Map<Integer,List<UngappedSearchHit>> hitsByBin = new TreeMap<Integer, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:sequenceHits) {
			int estStart = estimateSubjectStart(hit);
			int bin = estStart/1000;
			if(estStart<0) bin--;
			//if(hit.getQueryIdx()==10851) System.out.println("Initial hit clustering. Next hit start: "+hit.getStart()+" estimated: "+estStart+" cluster bin: "+bin);
			List<UngappedSearchHit> hitsBin = hitsByBin.computeIfAbsent(bin, (v)-> new ArrayList<UngappedSearchHit>());
			hitsBin.add(hit);
		}
		if(debug) System.out.println("Clustering kmer hits. Initial bin starts: "+hitsByBin.keySet()+" estimated number of clusters: "+estimatedClusters);
		//Second cluster centered in best averages
		int maxClusters = (int) Math.ceil(2*estimatedClusters)+1;
		Set<Integer> estimatedStarts = calculateEstimatedStarts(hitsByBin.values(), maxClusters, false);
		
		if (estimatedStarts.size()==0) return new ArrayList<>();
		estimatedStarts = clumpByDistance(estimatedStarts);
		hitsByBin = clusterByEstimatedStarts(sequenceHits, estimatedStarts);
		if(debug) System.out.println("Clustering kmer hits. Clusters after first relocation");
		if(debug) printClusters(hitsByBin);
		estimatedStarts = calculateEstimatedStarts (hitsByBin.values(), maxClusters, true);
		hitsByBin = clusterByEstimatedStarts(sequenceHits, estimatedStarts);
		if(debug) System.out.println("Clustering kmer hits. Clusters after second relocation");
		if(debug) printClusters(hitsByBin);
		List<List<UngappedSearchHit>> hitsClusters=new ArrayList<List<UngappedSearchHit>>(hitsByBin.values());
		return hitsClusters;
	}
	
	private Set<Integer> calculateEstimatedStarts(Collection<List<UngappedSearchHit>> clusters, int maxClusters, boolean useModeDiff) {
		List<List<UngappedSearchHit>> sortedClusters = new ArrayList<List<UngappedSearchHit>>(clusters);
		Collections.sort(sortedClusters,(c1,c2)->c2.size()-c1.size());
		Set<Integer> answer = new HashSet<>();
		for(int i=0;i<sortedClusters.size() && i<maxClusters;i++) {
			List<UngappedSearchHit> cluster = sortedClusters.get(i);
			double [] stats = calculateStatsEstimatedSubjectStart(cluster);
			int avg = (int) Math.round(stats[0]);
			int modeDiff = (int) stats[4];
			if(useModeDiff && Math.abs(modeDiff)>100) {
				answer.add(avg+modeDiff);
			} else {
				answer.add(avg);
			}
		}
		return answer;
	}
	private Set<Integer> clumpByDistance(Set<Integer> clusterAverages) {
		Set<Integer> answer = new HashSet<>();
		List<Integer> avgsList = new ArrayList<>(clusterAverages);
		Collections.sort(avgsList, (c1,c2)->c1-c2);
		int next = avgsList.get(0);
		for(int average:avgsList) {
			if(average-next<500) next = (next+average)/2;
			else {
				answer.add(next);
				next = average;
			}
		}
		answer.add(next);
		return answer;
	}
	private Map<Integer, List<UngappedSearchHit>> clusterByEstimatedStarts(List<UngappedSearchHit> sequenceHits, Set<Integer> estimatedStarts) {
		Map<Integer, List<UngappedSearchHit>> hitsByBin = new TreeMap<>();
		for(int start:estimatedStarts) hitsByBin.put(start, new ArrayList<>());
		for(UngappedSearchHit hit:sequenceHits) {
			int estStart = estimateSubjectStart(hit);
			int minS = 0;
			int minD = 500;
			for(int start:hitsByBin.keySet()) {
				int d = Math.abs(estStart-start);
				if(d<minD) {
					minS = start;
					minD = d;
				}
			}
			//if(hit.getQueryIdx()==10851) System.out.println("Clustering kmer hits. Next hit: "+hit.getQueryIdx()+" estimated start "+estStart+" min s: "+minS+" minD: "+minD);
			//if(subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) 
			if(minD<500) hitsByBin.get(minS).add(hit);
		}
		return hitsByBin;
	}
	private List<List<UngappedSearchHit>> clusterRegionKmerAlnsKruskal(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> sequenceHits, double estimatedClusters) {
		int n = sequenceHits.size();
		DisjointSets sets = new DisjointSets(n);
		//int [] estimates = new int [n];
		//for(int i=0;i<n;i++) estimates[i] = estimateSubjectStart(sequenceHits.get(i));
		List<UngappedSearchHit> sortedHits = new ArrayList<>(sequenceHits);
		Collections.sort(sortedHits,(h1,h2)->estimateSubjectStart(h1)-estimateSubjectStart(h2));
		int [] distanceNext = new int [n];
		int estP = estimateSubjectStart(sortedHits.get(0));
		List<Integer> sortedPos = new ArrayList<>(n-1);
		for(int i=1;i<n;i++) {
			sortedPos.add(i-1);
			int estI = estimateSubjectStart(sortedHits.get(i));
			distanceNext[i-1] = estI-estP;
			estP = estI;
		}
 		Collections.sort(sortedPos, (p1,p2)->distanceNext[p1]-distanceNext[p2]);
 		for (int i:sortedPos) {
 			int d = distanceNext[i];
 			if(d>50) break;
 			if(!sets.sameSubsets(i, i+1)) sets.union(i, i+1);
 			//if(d>0 && sets.getNumSubsets()<=2*estimatedClusters) break;
 		}
 		List<List<UngappedSearchHit>> initialClusters = buildClustersFromSubsets(sortedHits,sets);
 		if(debug) System.out.println("Clustering kmer hits. Clusters kruskal");
		if(debug) printClusters(initialClusters);
		return initialClusters;
	}
	private List<List<UngappedSearchHit>> buildClustersFromSubsets(List<UngappedSearchHit> sortedHits,DisjointSets sets) {
		Map<Integer,Set<Integer>> subsets = sets.getSubsets();
		List<List<UngappedSearchHit>> answer = new ArrayList<>(subsets.size());
		for(Set<Integer> subset:subsets.values()) {
			List<UngappedSearchHit> cluster = new ArrayList<>(subset.size());
			for(int i:subset) {
				cluster.add(sortedHits.get(i));
			}
			answer.add(cluster);
		}
		return answer;
	}
	private List<UngappedSearchHitsCluster> buildClusterObjectsFromClusteredHits(int queryLength, int subjectIdx,int subjectLength, List<List<UngappedSearchHit>> hitsClusters, double minHits, int rawHits) {
		Collections.sort(hitsClusters,(l1,l2)->l2.size()-l1.size());
		List<UngappedSearchHitsCluster> answer = new ArrayList<>();
		for(List<UngappedSearchHit> hits:hitsClusters) {
			if(debug) System.out.println("Next candidate cluster size: "+hits.size()+" current limit: "+minHits);
			if(hits.size()<minHits) break;
			List<List<UngappedSearchHit>> subclusters = breakByQueryStart(hits);
			if(subclusters.size()>1) System.err.println("WARN. Cluster broken in subclusters by query starts. SubjectIdx: "+subjectIdx+" query length: "+queryLength);
			for(List<UngappedSearchHit> subcluster:subclusters) {
				if(subcluster.size()<minHits) continue;
				List<UngappedSearchHit> selectedHits = collapseAndSelectSortedHits(queryLength, subjectIdx, subjectLength, subcluster);
				if(selectedHits.size()<minHits) continue;
				UngappedSearchHitsCluster cluster = new UngappedSearchHitsCluster(queryLength, subjectIdx, subjectLength, selectedHits);
				cluster.setRawKmerHits(rawHits);
				cluster.setRawKmerHitsSubjectStartSD(1);
				if(debug) System.out.println("Next cluster subject predicted coords: "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" subject evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" query evidence: "+cluster.getQueryEvidenceStart()+" "+cluster.getQueryEvidenceEnd()+" unique kmers: "+cluster.getNumDifferentKmers());
				if(cluster.getNumDifferentKmers()>=minHits) answer.add(cluster);
				minHits = Math.max(minHits, cluster.getNumDifferentKmers()/2);
			}
		}
		return answer;
	}
	private void printClusters(Map<Integer, List<UngappedSearchHit>> hitsByBin) {
		for(Map.Entry<Integer, List<UngappedSearchHit>> entry:hitsByBin.entrySet()) {
			double [] stats = calculateStatsEstimatedSubjectStart(entry.getValue());
			System.out.println("Start: "+entry.getKey()+" size: "+entry.getValue().size()+" stats: "+stats[0]+" "+stats[1]+" "+stats[3]);
		}
	}
	private void printClusters(List<List<UngappedSearchHit>> clusters) {
		List<List<UngappedSearchHit>> sortedClusters = new ArrayList<>(clusters);
		Collections.sort(sortedClusters,(c1,c2)->(int)(calculateStatsEstimatedSubjectStart(c1)[0])-(int)(calculateStatsEstimatedSubjectStart(c2)[0]));
		for(List<UngappedSearchHit> cluster:sortedClusters) {
			double [] stats = calculateStatsEstimatedSubjectStart(cluster);
			System.out.println("Size: "+cluster.size()+" stats: "+stats[0]+" "+stats[1]+" "+stats[3]);
		}
	}
	private List<List<UngappedSearchHit>> breakByQueryStart(List<UngappedSearchHit> hits) {
		List<List<UngappedSearchHit>> answer = new ArrayList<List<UngappedSearchHit>>();
		Collections.sort(hits, (h1,h2)->h1.getQueryStart()-h2.getQueryStart());
		List<UngappedSearchHit> nextSubcluster = new ArrayList<UngappedSearchHit>();
		int lastQueryStart = 0;
		for(UngappedSearchHit hit:hits) {
			//TODO: define better this parameter
			if(hit.getQueryStart()>lastQueryStart+30000) {
				if(nextSubcluster.size()>0) answer.add(nextSubcluster); 
				nextSubcluster = new ArrayList<UngappedSearchHit>();
			}
			nextSubcluster.add(hit);
			lastQueryStart = hit.getQueryStart();
		}
		if(nextSubcluster.size()>0) answer.add(nextSubcluster);
		return answer;
	}
	private double [] calculateStatsEstimatedSubjectStart(List<UngappedSearchHit> clusterHits) {
		double sum = 0;
		int n=0;
		List<Integer> starts = new ArrayList<>();
		for(UngappedSearchHit hit:clusterHits) {
			int start = estimateSubjectStart(hit);
			starts.add(start);
			sum+=start;
			n++;
		}
		Collections.sort(starts);
		double [] answer = {0,0,1,1,0};
		if(n>0) {
			answer[0] = sum/n;
			answer[1] = starts.get(n/2);
			if(n%2==0) {
				answer[1] += starts.get(n/2-1);
				answer[1] /=2;
			}
		}
		if(n>1) {
			//Better than sum squares because numbers are big
			double sumd2 = 0;
			Distribution dist = new Distribution(-500, 500, 50);
			for(int start:starts) {
				double d = start-answer[0];
				dist.processDatapoint(d);
				sumd2+=d*d;
			}
			answer[2] = sumd2/(n-1);
			//if(answer[2]<0) answer[2]=0;
			answer[3] = Math.sqrt(answer[2]);
			int modeDist = (int) Math.round(dist.getLocalMode(-500, 500));
			answer[4] = modeDist;
		}
		return answer;
	}
	private int estimateSubjectStart(UngappedSearchHit hit) {
		return hit.getSubjectStart() - hit.getQueryStart();
	}
	private int estimateQueryStart(UngappedSearchHit hit) {
		return hit.getQueryStart() - hit.getSubjectStart();
	}
	/**
	 * 
	 * @param queryLength
	 * @param subjectLength
	 * @param hits sorted by queryidx
	 * @return
	 */
	private List<UngappedSearchHit> collapseAndSelectSortedHits(int queryLength, int subjectIdx, int subjectLength, List<UngappedSearchHit> inputHits) {
		
		if(debug) System.out.println("KmerHitsCluster. Clustering "+inputHits.size()+" hits. Subject idx: "+subjectIdx);
		
		
		double [] stats = calculateStatsEstimatedSubjectStart(inputHits);
		
		int median = (int) Math.round(stats[1]);
		//System.out.println("Sum: "+sum+" sum2: "+sum2);
		double variance = stats[2];
		double rawKmerHitsSubjectStartSD = Math.max(1, stats[3]);
		Distribution dist = new Distribution(-500, 500, 50);
		Distribution distAbs = new Distribution(0, 500, 1);
		for(UngappedSearchHit hit:inputHits) {
			int start = estimateSubjectStart(hit);
			int distance = start-median;
			int distanceAbs = Math.abs(distance);
			if (distanceAbs < 2*rawKmerHitsSubjectStartSD) {
				dist.processDatapoint(distance);
				distAbs.processDatapoint(distanceAbs);
			}
		}
		int modeDist = (int) Math.round(dist.getLocalMode(-400, 400));
		if(Math.abs(modeDist)>=100) median+=(modeDist/2);
		
		//System.out.println(subjectStarts);
		if(debug) dist.printDistributionInt(System.out);
		int maxDistance = (int) Math.max(distAbs.getAverage(), Math.sqrt(distAbs.getVariance()+1));
		maxDistance *=5;
		maxDistance=Math.max(queryLength/40,maxDistance);
		maxDistance=Math.min(300,maxDistance);
		maxDistance=Math.max(10,maxDistance);
		//if(maxDistance<0.01*query.length()) maxDistance*=2;
		if(debug) System.out.println("KmerHitsCluster. Num hits: "+inputHits.size()+" median: "+median+" average: "+stats[0]+" variance: "+variance+" stdev: "+rawKmerHitsSubjectStartSD+" abs distance avg: "+distAbs.getAverage()+" stdev "+Math.sqrt(distAbs.getVariance())+" max distance: "+maxDistance);
		//Index hits by query kmer start
		Map<Integer, List<UngappedSearchHit>> hitsMultiMap = indexByQueryStart(inputHits);
		if(debug) System.out.println("KmerHitsCluster. Num different kmers: "+hitsMultiMap.size());
		List<UngappedSearchHit> selectedHits = selectHitsByDistanceWithMedian(hitsMultiMap, median, maxDistance);
		if(selectedHits.size()<1) { 
			if(debug) System.err.println("WARN. Empty list of selected hits for subject: "+subjectIdx);
			return selectedHits;
		}
		Set<Integer> outliersToRemove = replaceHitsByLocalAgreement(selectedHits, hitsMultiMap, median, maxDistance, queryLength);
		List<UngappedSearchHit> filteredHits = new ArrayList<UngappedSearchHit>();
		//if(outliersToRemove.size()>0 && outliersToRemove.size()<0.2*selectedHits.size()) {
			//for(int i=0;i<selectedHits.size();i++) {
				//if(!outliersToRemove.contains(i)) filteredHits.add(selectedHits.get(i));
				//else if (debug) System.out.println("Removing outlier at qpos "+selectedHits.get(i).getQueryStart()+" hit: "+selectedHits.get(i).getSubjectStart());
				
			//}
		//} else filteredHits.addAll(selectedHits);
		filteredHits.addAll(selectedHits);
		filteredHits = removeDisorganized (filteredHits, median);
		if(filteredHits.size()<1) {
			if(debug) System.err.println("WARN. Empty list of sorted hits for subject: "+subjectIdx+" selected hits: "+selectedHits.size()+" query length: "+queryLength);
			
		}
		return filteredHits;
	}
	private Map<Integer, List<UngappedSearchHit>> indexByQueryStart(List<UngappedSearchHit> inputHits) {
		Map<Integer,List<UngappedSearchHit>> hitsMultiMap = new TreeMap<Integer, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:inputHits) {
			//if (subjectIdx==idxSubjectDebug && queryLength == queryLengthDebug) System.out.println("Next qpos "+hit.getQueryIdx()+" hit: "+hit.getStart()+" kmer: "+hit.getQuery()+ " estq: "+estimateQueryStart(hit)+" - "+estimateQueryEnd(hit)+" estS: "+estimateSubjectStart(hit)+" - "+estimateSubjectEnd(hit));
			List<UngappedSearchHit> list = hitsMultiMap.computeIfAbsent(hit.getQueryStart(), l -> new ArrayList<UngappedSearchHit>());
			list.add(hit);
		}
		return hitsMultiMap;
	}
	private List<UngappedSearchHit> selectHitsByDistanceWithMedian(Map<Integer, List<UngappedSearchHit>> hitsMultiMap, int median, int maxDistance) {
		List<UngappedSearchHit> selectedHits = new ArrayList<UngappedSearchHit>(hitsMultiMap.size());
		for(List<UngappedSearchHit> hits:hitsMultiMap.values()) {
			UngappedSearchHit hit = selectHit(hits, median, maxDistance);
			if(hit!=null) {
				if (debug) System.out.println("Selected hits. Next qpos "+hit.getQueryStart()+" hit: "+hit.getSubjectStart()+" estq: "+estimateQueryStart(hit)+" estS: "+estimateSubjectStart(hit)+" all starts: "+calculateHitStarts(hits));
				selectedHits.add(hit);
			}
		}
		return selectedHits;
	}

	private List<Integer> calculateHitStarts(List<UngappedSearchHit> hits) {
		List<Integer> hitStarts = new ArrayList<Integer>();
		for(UngappedSearchHit hit2:hits) hitStarts.add(estimateSubjectStart(hit2));
		return hitStarts;
	}
	
	private UngappedSearchHit selectHit(List<UngappedSearchHit> hits, int median, int maxDistance) {
		UngappedSearchHit answer = null;
		int minCost = 0;
		
		for(UngappedSearchHit hit:hits) {
			int estStart = estimateSubjectStart(hit);
			int distance = Math.abs(estStart-median);
			int cost = distance;
			//if(debug && hit.getQueryIdx()==5314) System.out.println("Selecting hit. Next hit: "+hit.getStart()+" distance: "+distance+" max "+maxDistance+" cost: "+cost);
			if(distance <= maxDistance) {
				if(answer == null || minCost>cost) {
					answer = hit;
					minCost = cost;
				}
			}
		}
		return answer;
	}
	
	private Set<Integer> replaceHitsByLocalAgreement(List<UngappedSearchHit> selectedHits, Map<Integer, List<UngappedSearchHit>> hitsMultiMap, int median, int maxDistance, int queryLength) {
		//Find trustable site
		int minHitPos = -1;
		int minCost = -1;
		int midPoint = queryLength/2;
		for(int i=0;i<selectedHits.size();i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			int estStart = estimateSubjectStart(hit);
			int distance = Math.abs(estStart-median);
			int d2 = Math.abs(hit.getQueryStart()-midPoint);
			int cost = distance + d2/10;
			if(minHitPos==-1 || cost<minCost) {
				minHitPos = i;
				minCost = cost;
			}
		}
		if(minHitPos==-1) return new HashSet<Integer>();
		UngappedSearchHit minHit = selectedHits.get(minHitPos);
		if (debug) System.out.println("Hit closest to median. qpos "+minHit.getQueryStart()+" hit: "+minHit.getSubjectStart()+" cost: "+minCost+" estq: "+estimateQueryStart(minHit)+" estS: "+estimateSubjectStart(minHit));
		int vicinityEstStart = estimateSubjectStart(minHit);
		LinkedList<Integer> localValues = new LinkedList<Integer>();
		localValues.add(vicinityEstStart);
		Set<Integer> posHitsToRemove = new HashSet<Integer>();
		for(int i=minHitPos-1;i>=0;i--) {
			UngappedSearchHit nextHit = selectedHits.get(i);
			int estStart = estimateSubjectStart(nextHit);
			if(Math.abs(estStart-vicinityEstStart)> 5) {
				List<UngappedSearchHit> hits = hitsMultiMap.get(nextHit.getQueryStart());
				vicinityEstStart = getMedian(localValues);
				UngappedSearchHit updatedHit = selectHit(hits, vicinityEstStart, 10);
				if (debug) System.out.println("Trying to replace outlier at qpos "+nextHit.getQueryStart()+" hit: "+nextHit.getSubjectStart()+" ests: "+estStart+" local median: "+vicinityEstStart+" all starts: "+calculateHitStarts(hits));
				if(updatedHit!=null && updatedHit!=nextHit) {
					if (debug) System.out.println("Replacing hit. qpos "+updatedHit.getQueryStart()+" hit: "+updatedHit.getSubjectStart()+" estq: "+estimateQueryStart(updatedHit)+" estS: "+estimateSubjectStart(updatedHit));
					selectedHits.set(i, updatedHit);
					estStart = estimateSubjectStart(updatedHit);
				} else {
					posHitsToRemove.add(i);
				}
			} else {
				vicinityEstStart = estStart;
			}
			localValues.add(estStart);
			if(localValues.size()>10) localValues.removeFirst();
		}
		vicinityEstStart = estimateSubjectStart(minHit);
		localValues.clear();
		localValues.add(vicinityEstStart);
		for(int i=minHitPos+1;i<selectedHits.size();i++) {
			UngappedSearchHit nextHit = selectedHits.get(i);
			int estStart = estimateSubjectStart(nextHit);
			if(Math.abs(estStart-vicinityEstStart)> 5) {
				List<UngappedSearchHit> hits = hitsMultiMap.get(nextHit.getQueryStart());
				vicinityEstStart = getMedian(localValues);
				UngappedSearchHit updatedHit = selectHit(hits, vicinityEstStart, 10);
				if (debug) System.out.println("Trying to replace outlier at qpos "+nextHit.getQueryStart()+" hit: "+nextHit.getSubjectStart()+" ests: "+estStart+" local median: "+vicinityEstStart+" all starts: "+calculateHitStarts(hits));
				if(updatedHit!=null && updatedHit!=nextHit) {
					if (debug) System.out.println("Replacing hit. qpos "+updatedHit.getQueryStart()+" hit: "+updatedHit.getSubjectStart()+" estq: "+estimateQueryStart(updatedHit)+" estS: "+estimateSubjectStart(updatedHit));
					selectedHits.set(i, updatedHit);
					estStart = estimateSubjectStart(updatedHit);
				} else {
					posHitsToRemove.add(i);
				}
			} else {
				vicinityEstStart = estStart;
			}
			localValues.add(estStart);
			if(localValues.size()>10) localValues.removeFirst();
		}
		return posHitsToRemove;
	}

	private int getMedian(LinkedList<Integer> localValues) {
		ArrayList<Integer> copy = new ArrayList<Integer>(localValues);
		Collections.sort(copy);
		return copy.get(copy.size()/2);
	}

	private List<UngappedSearchHit> removeDisorganized(List<UngappedSearchHit> selectedHits, int median) {
		int n = selectedHits.size();
		//if(n<=30) return selectSorted(selectedHits, 0, n-1, median);
		boolean [] vicinityConsistent = new boolean [n];
		Arrays.fill(vicinityConsistent, false);
		for(int i=10;i<n-10;i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			int start = hit.getSubjectStart();
			boolean consistent = true;
			for(int j=i-1;j>=i-10 && consistent;j--) {
				UngappedSearchHit hit2 = selectedHits.get(j);
				if(hit2.getSubjectStart()>=start) consistent = false; 
			}
			for(int j=i+1;j<=i+10 && consistent;j++) {
				UngappedSearchHit hit2 = selectedHits.get(j);
				if(hit2.getSubjectStart()<=start) consistent = false; 
			}
			vicinityConsistent[i] = consistent;
			if(!consistent) {
				for(int j=i-1;j>=i-10;j--) vicinityConsistent[j] = false;
				for(int j=i+1;j<=i+10 ;j++) vicinityConsistent[j] = false;
				i+=10;
			}
		}
		List<UngappedSearchHit> dpSortedHits = new ArrayList<UngappedSearchHit>();
		double sum = 0;
		double sum2 = 0;
		double nS = 0;
		int i=0;
		while(i<n) {
			UngappedSearchHit hit = selectedHits.get(i);
			if(vicinityConsistent[i]) {
				//if(debug) System.out.println("Adding consistent hit at "+hit.getQueryIdx());
				dpSortedHits.add(hit);
				double s = estimateSubjectStart(hit);
				sum+=s;
				sum2+=s*s;
				nS++;
				i++;
			} else {
				int j;
				for(j=i+1;j<n && !vicinityConsistent[j];j++);
				dpSortedHits.addAll(selectSorted(selectedHits,i,j-1,median));
				i=j;
			}
		}
		//Remove final outliers
		int n2= dpSortedHits.size();
		if(debug) System.out.println("Initial hits: "+n+" consistent: "+nS+" dpSorted: "+n2);
		if(n2<30 || nS < 20) return dpSortedHits;
		int averageStart = (int) (sum/nS);
		double variance = Math.max(1, (sum2-sum*sum/nS)/(nS-1));
		double stdev = Math.sqrt(variance);
		double maxDistance = Math.max(100, 3*stdev);
		int countOutliers =0;
		for(UngappedSearchHit hit:dpSortedHits) {
			double distance = Math.abs(estimateSubjectStart(hit)-averageStart);
			if(distance > maxDistance) countOutliers++;
		}
		
		if(debug) System.out.println("Average start: "+averageStart+" stdev: "+stdev+" max distance: "+maxDistance+" outliers: "+countOutliers+" n "+n2+" nconsistent: "+nS);
		List<UngappedSearchHit> answer = new ArrayList<UngappedSearchHit>(n2);
		if(countOutliers<0.02*n2) {
			for(UngappedSearchHit hit:dpSortedHits) {
				double distance = Math.abs(estimateSubjectStart(hit)-averageStart);
				if(distance <= maxDistance) answer.add(hit);
				else if(debug) System.out.println("Removing outlier: "+hit.getQueryStart()+" subject: "+hit.getSubjectStart()+" estimated start: "+estimateSubjectStart(hit));
			}
			return answer;
		}
		Collections.sort(dpSortedHits,(h1,h2)->Math.abs(estimateSubjectStart(h2)-averageStart)-Math.abs(estimateSubjectStart(h1)-averageStart));
		//if(debug) System.out.println("Sorted by distance");
		for (int k=0;k<dpSortedHits.size();k++) {
			UngappedSearchHit hit = dpSortedHits.get(k);
			double distance = Math.abs(estimateSubjectStart(hit)-averageStart);
			if(distance <= maxDistance || k>=0.02*n2) answer.add(hit);
			else if(debug) System.out.println("Removing outlier after sorting: "+hit.getQueryStart()+" subject: "+hit.getSubjectStart()+" estimated start: "+estimateSubjectStart(hit));
		}
		Collections.sort(answer,(h1,h2)->h1.getQueryStart()-h2.getQueryStart());
		return answer;
	}


	private List<UngappedSearchHit> selectSorted(List<UngappedSearchHit> selectedHits, int first, int last, int median) {
		if(debug) System.out.println("DP Sorting hits from "+first+" to "+last+" median start: "+median);
		//Check sorted first
		List<UngappedSearchHit> regionHits = new ArrayList<UngappedSearchHit>(last-first+1);
		int lastStartSubject = -1;
		int lastStartQuery = -1;
		boolean sortedSubject = true;
		boolean sortedQuery = true;
		for(int i=first;i<=last;i++) {
			UngappedSearchHit hit = selectedHits.get(i);
			regionHits.add(hit);
			int startS = hit.getSubjectStart();
			if(lastStartSubject>=startS) sortedSubject = false;
			lastStartSubject = startS;
			int startQ = hit.getQueryStart();
			if(lastStartQuery>=startQ) sortedQuery = false;
			lastStartQuery = startQ;
		}
		if(!sortedQuery) {
			System.err.println("WARN. Disorganized query hits list");
			Collections.sort(regionHits,(h1,h2)->h1.getQueryStart()-h2.getQueryStart());
			
		} else if (sortedSubject) return regionHits;
		List<UngappedSearchHit> regionHitsBySubject = new ArrayList<UngappedSearchHit>();
		regionHitsBySubject.addAll(regionHits);
		Collections.sort(regionHitsBySubject,(h1,h2)->h1.getSubjectStart()-h2.getSubjectStart());
		int n = regionHits.size();
		//Dynamic programming to select the best organized subset
		int [][] scores = new int [n+1][n+1];
		for(int i=0;i<=n;i++) {
			for(int j=0;j<=n;j++) {
				if(i==0 || j==0) scores[i][j] = 0;
				else {
					UngappedSearchHit hit1 = regionHits.get(i-1);
					UngappedSearchHit hit2 = regionHitsBySubject.get(j-1);
					if(hit1==hit2) scores[i][j] = calculateScore(hit1,median) + scores[i-1][j-1];
					else scores[i][j] = Math.max(scores[i-1][j], scores[i][j-1]);
				}
			}
			 
		}
		//if(debug) printMatrix(scores);
		List<UngappedSearchHit> answer = new ArrayList<UngappedSearchHit>();
		int i= n;
		int j = n;
		while(i>0 && j>0) {
			UngappedSearchHit hit1 = regionHits.get(i-1);
			UngappedSearchHit hit2 = regionHitsBySubject.get(j-1);
			if(hit1==hit2) {
				//if(debug) System.out.println("Added DP sorted hit at "+hit1.getQueryIdx());
				answer.add(hit1);
				i--;
				j--;
			} else if (scores[i][j] == scores[i-1][j]) i--;
			else j--;
		}
		Collections.reverse(answer);
		return answer;
	}
	public void printMatrix(int[][] scores) {
		for(int i=0;i<scores.length;i++) {
			for(int j=0;j<scores[i].length;j++) {
				if(j>0) System.out.print("\t");
				System.out.print(""+scores[i][j]);
			}
			System.out.println();
		}
		
	}
	private int calculateScore(UngappedSearchHit hit, int median) {
		/*int start = estimateSubjectStart(hit);
		int distance = Math.abs(start-median);
		int score = (int) (1000*Math.max(0, 1-((double)distance)/1000.0));
		if(debug) System.out.println("Score for hit: "+hit.getQueryIdx()+" "+score+" start. "+start+" distance: "+distance);
		return score;*/
		return 1;
	}
	
}
