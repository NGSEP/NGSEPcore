package ngsep.genome;

import java.util.ArrayList;
import java.util.List;

public class SyntenyVertex {
	private HomologyEdge homologyRelationship;
	private int weight;
	private List<SyntenyEdge> edges = new ArrayList<SyntenyEdge>();
	public SyntenyVertex(HomologyEdge homologyRelationship, int weight) {
		super();
		this.homologyRelationship = homologyRelationship;
		this.weight = weight;
	}
	public int getWeight() {
		return weight;
	}
	public void setWeight(int weight) {
		this.weight = weight;
	}
	public HomologyEdge getHomologyRelationship() {
		return homologyRelationship;
	}
	public void addEdge (SyntenyEdge edge) {
		edges.add(edge);
	}
	
	public List<SyntenyEdge> getEdges() {
		return edges;
	}
	
	
}
