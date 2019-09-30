package ngsep.clustering;

public class DendrogramEdge {
	private double weight;
	private Dendrogram origin;
	private Dendrogram destination;
	public DendrogramEdge(double weight, Dendrogram origin, Dendrogram destination) {
		super();
		this.weight = weight;
		this.origin = origin;
		this.destination = destination;
	}
	/**
	 * @return the weight
	 */
	public double getWeight() {
		return weight;
	}
	/**
	 * @return the origin
	 */
	public Dendrogram getOrigin() {
		return origin;
	}
	/**
	 * @return the destination
	 */
	public Dendrogram getDestination() {
		return destination;
	}
	
	
	
}
