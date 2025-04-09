package ngsep.alignments;

import java.util.ArrayList;
import java.util.Collections;
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
			if(debug) System.out.println("Next: "+hit.getQueryStart()+" subjectPos: "+hit.getSubjectStart()+" estimate "+estI+" distance: "+distanceNext[i-1]);
			estP = estI;
		}
 		Collections.sort(sortedPos, (p1,p2)->distanceNext[p1]-distanceNext[p2]);
 		for (int i:sortedPos) {
 			int d = distanceNext[i];
 			if(d>30) break;
 			if(!sets.sameSubsets(i, i+1)) sets.union(i, i+1);
 			//if(d>0 && sets.getNumSubsets()<=2*estimatedClusters) break;
 		}
 		List<List<UngappedSearchHit>> initialClusters = buildClustersFromSubsets(sortedHits,sets);
 		//if(debug) System.out.println("Clustering kmer hits. Clusters kruskal");
		//if(debug) printClusters(initialClusters);
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

}
