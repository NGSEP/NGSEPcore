package ngsep.clustering;

public interface DistanceMatrixClustering {
	/**
	 * Builds a dendrogram representing the clustering given the current distance matrix
	 * @param distances Matrix of distances
	 * @return Dendrogram representing the clusters
	 */
	public Dendrogram buildDendrogram(DistanceMatrix distances);
}
