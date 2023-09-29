package ngsep.clustering;

public class DendrogramEdge {

	private final double weight;
	private final Dendrogram destination;

	public DendrogramEdge(double weight, Dendrogram destination) {
		this.weight = weight;
		this.destination = destination;
	}
	/**
	 * @return the weight
	 */
	public double getWeight() {
		return weight;
	}

	/**
	 * @return the destination
	 */
	public Dendrogram getDestination() {
		return destination;
	}
}
