package ngsep.clustering;

public class DendrogramEdge {

	private double weight;
	private Dendrogram destination;

	public DendrogramEdge(double weight, Dendrogram destination) {
		super();
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
