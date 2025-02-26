package ngsep.alignments;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import ngsep.sequences.UngappedSearchHit;

public class UngappedSearchHitClusteringAlgorithmKmeansLike implements UngappedSearchHitClusteringAlgorithm {

	@Override
	public List<List<UngappedSearchHit>> clusterLocalSearchHits(List<UngappedSearchHit> hits) {
		boolean debug = false;
		//UngappedSearchHit firstHit = sequenceHits.get(0);
		//int subjectIdx = firstHit.getSequenceIdx();
		//Initial clustering
		Map<Integer,Integer> countsByQueryIdx = new HashMap<Integer, Integer>();
		//if(debug) System.out.println("Clustering hits: "+sequenceHits.size());
		for(UngappedSearchHit hit:hits) {
			//if (debug && hit.getSubjectStart()<10000) System.out.println("Next qpos "+hit.getQueryStart()+" hit: "+hit.getSubjectStart());
			countsByQueryIdx.compute(hit.getQueryStart(), (k,v)->v==null?1:v+1);
		}
		//double estimatedClusters = 0.5*sequenceHits.size()/queryLength;
		double avg = 0;
		for(int count:countsByQueryIdx.values()) avg+=count;
		if(avg>0) avg/=countsByQueryIdx.size();
		double estimatedClusters = avg;
		 
		if(debug) System.out.println("Clustering hits: "+hits.size()+" query starts: "+countsByQueryIdx.size()+" estimatedClusters: "+estimatedClusters);
		Map<Integer,List<UngappedSearchHit>> hitsByBin = new TreeMap<Integer, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:hits) {
			int estStart = hit.estimateSubjectStart();
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
		hitsByBin = clusterByEstimatedStarts(hits, estimatedStarts);
		if(debug) System.out.println("Clustering kmer hits. Clusters after first relocation");
		if(debug) printHitsByBin(hitsByBin);
		estimatedStarts = calculateEstimatedStarts (hitsByBin.values(), maxClusters, true);
		hitsByBin = clusterByEstimatedStarts(hits, estimatedStarts);
		if(debug) System.out.println("Clustering kmer hits. Clusters after second relocation");
		if(debug) printHitsByBin(hitsByBin);
		List<List<UngappedSearchHit>> hitsClusters=new ArrayList<List<UngappedSearchHit>>(hitsByBin.values());
		return hitsClusters;
	}
	
	private Set<Integer> calculateEstimatedStarts(Collection<List<UngappedSearchHit>> clusters, int maxClusters, boolean useModeDiff) {
		List<List<UngappedSearchHit>> sortedClusters = new ArrayList<List<UngappedSearchHit>>(clusters);
		Collections.sort(sortedClusters,(c1,c2)->c2.size()-c1.size());
		Set<Integer> answer = new HashSet<>();
		for(int i=0;i<sortedClusters.size() && i<maxClusters;i++) {
			List<UngappedSearchHit> cluster = sortedClusters.get(i);
			double [] stats = UngappedSearchHitsClusterBuilder.calculateStatsEstimatedSubjectStart(cluster);
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
			int estStart = hit.estimateSubjectStart();
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
	private void printHitsByBin(Map<Integer, List<UngappedSearchHit>> hitsByBin) {
		for(Map.Entry<Integer, List<UngappedSearchHit>> entry:hitsByBin.entrySet()) {
			double [] stats = UngappedSearchHitsClusterBuilder.calculateStatsEstimatedSubjectStart(entry.getValue());
			System.out.println("Start: "+entry.getKey()+" size: "+entry.getValue().size()+" stats: "+stats[0]+" "+stats[1]+" "+stats[3]);
		}
	}
}
