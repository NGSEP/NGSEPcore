package ngsep.clustering;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;


public class Dendrogram {

	private String label;
	private ArrayList<DendrogramEdge> children;

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

	public String toNewick(){
		if (children.isEmpty())
			return "();";
		else {
			DendrogramEdge firstL = this.children.get(0);
			DendrogramEdge firstR = this.children.get(1);
			Dendrogram firstLt = firstL.getDestination();
			Dendrogram firstRt = firstR.getDestination();
			double ld = firstL.getWeight();
			double rd = firstR.getWeight();

			return String.format(Locale.ROOT, "(%s:%f,%s:%f);", toNewick(firstLt), ld, toNewick(firstRt), rd);
		}
	}

	private String toNewick(Dendrogram t){
		if (isLeaf(t)) return t.label;
		else {
			ArrayList<DendrogramEdge> currentChildren = t.children;
			String level = "(";
			for (int i = 0; i < currentChildren.size(); i++) {
				DendrogramEdge e = currentChildren.get(i);
				level += String.format(Locale.ROOT, "%s:%f", toNewick(e.getDestination()), e.getWeight());

				if (i != currentChildren.size() - 1)
					level += ",";
			}
			level += ")";
			return level;
		}
	}

	private boolean isLeaf(Dendrogram t){
		return t.children.isEmpty();
	}

	public ArrayList<DendrogramEdge> getChildren(){
		return this.children;
	}

	public void setChildren(ArrayList<DendrogramEdge> children){this.children = children;}
}