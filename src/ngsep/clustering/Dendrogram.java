package ngsep.clustering;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;



public class Dendrogram {

	private String label;
	private List<DendrogramEdge> children;

	public Dendrogram(String label) {
	    this.label = label;
	    children = new ArrayList<DendrogramEdge>();
	}
	public Dendrogram(String label, List<DendrogramEdge> children) {
	    this.label = label;
	    this.children = new ArrayList<DendrogramEdge>();
	    this.children.addAll(children);
	}
	
	public void printTree(final PrintStream ps) {    
		ps.print(label);
		for(DendrogramEdge edge:children) {
			Dendrogram child = edge.getDestination();
			child.printTree(ps);
		}
	}
}