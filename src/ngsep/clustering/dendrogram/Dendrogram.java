package ngsep.clustering.dendrogram;

import ngsep.clustering.Pair;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;


public class Dendrogram {

	private final String label;
	private ArrayList<DendrogramEdge> children;

	public Dendrogram(String label) {
		this.label = label;
		children = new ArrayList<>();
	}
	public Dendrogram(String label, List<DendrogramEdge> children) {
		this.label = label;
		this.children = new ArrayList<>();
		this.children.addAll(children);
	}

	public void printTree(PrintStream out) {
		out.println(this.toNewick());
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
			StringBuilder level = new StringBuilder("(");
			for (int i = 0; i < currentChildren.size(); i++) {
				DendrogramEdge e = currentChildren.get(i);
				level.append(String.format(Locale.ROOT, "%s:%f", toNewick(e.getDestination()), e.getWeight()));

				if (i != currentChildren.size() - 1)
					level.append(",");
			}
			level.append(")");
			return level.toString();
		}
	}

	private boolean isLeaf(Dendrogram t){
		return t.children.isEmpty();
	}

	public void setChildren(ArrayList<DendrogramEdge> children){this.children = children;}

	/**
	 * Joins a pair of nodes (u, v) to a new node x. Returns the new
	 * tree with x as the root
	 * @param name - Label of the new tree x
	 * @param distances - A pair (dux, dvx) with the distance between u and x, and the distance between v and x respectively
	 * @param nodes - The pair of trees (u, v) to be joined
	 * @return A new tree with x as the root
	 */
	public static Dendrogram join2 (
			String name,
			Pair<Double, Double> distances,
			Pair<Dendrogram, Dendrogram> nodes
	) {
		double dux = distances.first;
		double dvx = distances.second;
		Dendrogram unode = nodes.first;
		Dendrogram vnode = nodes.second;
		return new Dendrogram(name, List.of(
				new DendrogramEdge(dux, unode),
				new DendrogramEdge(dvx, vnode)
		));
	}
}