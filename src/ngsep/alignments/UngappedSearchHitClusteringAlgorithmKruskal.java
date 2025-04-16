package ngsep.alignments;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.graphs.DisjointSets;
import ngsep.sequences.UngappedSearchHit;

public class UngappedSearchHitClusteringAlgorithmKruskal implements UngappedSearchHitClusteringAlgorithm {

	private boolean debug = false;
	
	public boolean isDebug() {
		return debug;
	}
	public void setDebug(boolean debug) {
		this.debug = debug;
	}
	@Override
	public List<List<UngappedSearchHit>> clusterLocalSearchHits(List<UngappedSearchHit> hits) {
		int maxDistance = 80;
		List<List<UngappedSearchHit>> clusters = clusterLocalSearchHits(hits,maxDistance);
		List<List<UngappedSearchHit>> answer = new ArrayList<List<UngappedSearchHit>>();
		int rounds = 3;
		for(int i=0;i<rounds && clusters.size()>0;i++) {
			maxDistance/=2;
			List<List<UngappedSearchHit>> roundClusters = new ArrayList<List<UngappedSearchHit>>();
			for(List<UngappedSearchHit> cluster:clusters) {
				if(isRepetitive(cluster)) {
					roundClusters.addAll(clusterLocalSearchHits(cluster,maxDistance));
				} else {
					answer.add(cluster);
				}
			}
			if(debug) System.err.println("Clustering kmer hits. Clusters kruskal end round "+i+" with maxDistance "+maxDistance+". Initial: "+clusters.size()+" nonrep: "+answer.size()+" new: "+roundClusters.size());
			//if(debug) printClusters(clusters);
			clusters = roundClusters;
		}
 		answer.addAll(clusters);
 		Collections.sort(answer,(c1,c2)->c2.size()-c1.size());
 		int n = answer.size();
 		int limit = 0;
 		if(n>=100) {
 			limit = answer.get(0).size()/10;
 			answer = selectBySize(answer,limit);
 		}
 		if(debug && n>0) System.err.println("Clustering kmer hits. Before size selection: "+n+" first: "+answer.get(0).size()+"limit: "+limit+" after: "+answer.size());
		return answer;
	}
	private List<List<UngappedSearchHit>> clusterLocalSearchHits(List<UngappedSearchHit> hits, int maxDistance) {
		int n = hits.size();
		DisjointSets sets = new DisjointSets(n);
		//int [] estimates = new int [n];
		//for(int i=0;i<n;i++) estimates[i] = estimateSubjectStart(sequenceHits.get(i));
		List<UngappedSearchHit> sortedHits = new ArrayList<>(hits);
		Collections.sort(sortedHits,(h1,h2)->h1.estimateSubjectStart()-h2.estimateSubjectStart());
		int [] distanceNext = new int [n];
		int estP = sortedHits.get(0).estimateSubjectStart();
		List<Integer> sortedPos = new ArrayList<>(n-1);
		for(int i=1;i<n;i++) {
			sortedPos.add(i-1);
			UngappedSearchHit hit = sortedHits.get(i);
			int estI = hit.estimateSubjectStart();
			distanceNext[i-1] = estI-estP;
			if(debug && maxDistance==80) System.err.println("Next: "+hit.getQueryStart()+" subjectPos: "+hit.getSubjectStart()+" estimate "+estI+" distance: "+distanceNext[i-1]);
			estP = estI;
		}
 		Collections.sort(sortedPos, (p1,p2)->distanceNext[p1]-distanceNext[p2]);
 		for (int i:sortedPos) {
 			int d = distanceNext[i];
 			if(d>maxDistance) break;
 			if(!sets.sameSubsets(i, i+1)) sets.union(i, i+1);
 			//if(d>0 && sets.getNumSubsets()<=2*estimatedClusters) break;
 		}
 		List<List<UngappedSearchHit>> clusters = buildClustersFromSubsets(sortedHits,sets);
 		return clusters;
	}
	private boolean isRepetitive(List<UngappedSearchHit> cluster) {
		int n = cluster.size();
		if(n<20) return false;
		int d = 0;
		Map<Integer,Integer> qstarts = new HashMap<Integer,Integer>();
		for(UngappedSearchHit hit:cluster) {
			if(qstarts.containsKey(hit.getQueryStart())) d++;
			else qstarts.put(hit.getQueryStart(),hit.getSubjectStart());
		}
		return d>5 && d>0.2*n;
	}
	private List<List<UngappedSearchHit>> selectBySize(List<List<UngappedSearchHit>> clusters, int limit) {
		List<List<UngappedSearchHit>> answer = new ArrayList<List<UngappedSearchHit>>();
		for(List<UngappedSearchHit> cluster:clusters) {
			if(cluster.size()>=limit) answer.add(cluster);
		}
		return answer;
	}
	private void printClusters(List<List<UngappedSearchHit>> initialClusters) {
		System.out.println("PartitionClustering. Printing "+initialClusters.size()+" clusters");
		for(List<UngappedSearchHit> hits:initialClusters) {
			System.out.println("PartitionClustering. Printing cluster with "+hits.size()+" hits");
			long s = 0;
			int d = 0;
			int n=0;
			Map<Integer,Integer> qstarts = new HashMap<Integer,Integer>();
			for(UngappedSearchHit hit:hits) {
				int e = hit.estimateSubjectStart();
				n++;
				s+=e;
				if(qstarts.containsKey(hit.getQueryStart())) {
					d++;
					//System.out.println("PartitionClustering. Two hits for q start: "+hit.getQueryStart()+" : "+hit.getSubjectStart()+" "+qstarts.get(hit.getQueryStart()));
				}
				else qstarts.put(hit.getQueryStart(),hit.getSubjectStart());
			}
			System.out.println("PartitionClustering. Next cluster with "+n+" hits. Distinct starts: "+qstarts.size()+" avg estimate: "+(s/n)+" duplicated: "+d);
		}
		
		
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

}
