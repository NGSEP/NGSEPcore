package ngsep.genome;

import java.util.ArrayList;
import java.util.List;

public class SyntenyVertex {
	private HomologyCluster homologyCluster;
	private int weight;
	private List<SyntenyEdge> edges = new ArrayList<SyntenyEdge>();
	public SyntenyVertex(HomologyCluster homologyCluster) {
		super();
		this.homologyCluster = homologyCluster;
	}
	public int getWeight() {
		return weight;
	}
	public void setWeight(int weight) {
		this.weight = weight;
	}
	public HomologyCluster getHomologyCluster() {
		return homologyCluster;
	}
	public void addEdge (SyntenyEdge edge) {
		edges.add(edge);
	}
	
	public List<SyntenyEdge> getEdges() {
		return edges;
	}
	
	
}
