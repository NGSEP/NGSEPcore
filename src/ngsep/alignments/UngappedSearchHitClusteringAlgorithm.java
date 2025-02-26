package ngsep.alignments;

import java.util.List;

import ngsep.sequences.UngappedSearchHit;

public interface UngappedSearchHitClusteringAlgorithm {
	/**
	 * Build clusters from the given search hits
	 * @param hits to cluster
	 * @return List of List of hits with the calculated clusters
	 */
	public List<List<UngappedSearchHit>> clusterLocalSearchHits(List<UngappedSearchHit> hits);
}
