package ngsep.clustering;


import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Queue;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public class Dendrogram {

	private int id = -1;
	private final String label;
	private int size = 1;
	private List<DendrogramEdge> children;
	private boolean hasIds = false;

	public Dendrogram(String label) {
		this.label = label;
		children = new ArrayList<>();
	}

	public Dendrogram(String label, List<DendrogramEdge> children) {
		this.label = label;
		this.children = new ArrayList<>();
		this.setChildren(children);
	}

	public Dendrogram (InputStream in) throws IOException {
		Dendrogram tree = fromFile(in);
		this.label = tree.label;
		this.children = tree.children;
		this.size = tree.getSize();
		this.id = tree.id;
		this.hasIds = tree.hasIds;
	}

	public void setChildren(List<DendrogramEdge> children){
		this.children = children;
		for (DendrogramEdge e : children) {
			this.size += e.getDestination().getSize();
		}
	}

	public String getLabel () {
		return this.label;
	}

	public int getId () {
		return this.id;
	}

	public int getSize () {
		return this.size;
	}

	private boolean isLeaf(Dendrogram t){
		return t.children.isEmpty();
	}

	public void printTree(PrintStream out) {
		out.println(this.toNewick());
	}

	public void generateIds () {
		if (!hasIds) {
			Queue<Dendrogram> q = new LinkedList<>();
			q.add(this);
			this.id = 0;
			int id = 1;
			while (!q.isEmpty()) {
				Dendrogram v = q.remove();
				for (DendrogramEdge e : v.children) {
					Dendrogram u = e.getDestination();
					u.id = id++;
					q.add(u);
				}
			}
			this.hasIds = true;
		}
	}

	public List<Dendrogram> getLeaves () {
		List<Dendrogram> leaves = new ArrayList<>();
		Queue<Dendrogram> q = new LinkedList<>();
		q.add(this);
		while (!q.isEmpty()) {
			Dendrogram v = q.remove();
			if (v.children.isEmpty()) leaves.add(v);
			else {
				for (DendrogramEdge e : v.children) {
					Dendrogram u = e.getDestination();
					q.add(u);
				}
			}
		}
		return leaves;
	}

	private static List<String> extractNewickChildren (String newick) {
		// remove outer parenthesis
		String childrenStr = newick.substring(1, newick.length() - 1);
		int parenthesis = 0;
		StringBuilder currChildren = new StringBuilder();
		List<String> children = new ArrayList<>();

		for (int i = 0; i < childrenStr.length(); i++) {
			char c = childrenStr.charAt(i);
			if (c == '(') {
				parenthesis++;
				currChildren.append(c);
			}
			else if (c == ')') {
				parenthesis--;
				currChildren.append(c);
			}
			else if (c == ',' && parenthesis == 0) {
				children.add(currChildren.toString());
				currChildren = new StringBuilder();
			} else {
				currChildren.append(c);
			}
		}
		if (currChildren.length() > 0) children.add(currChildren.toString());

		return children;
	}

	private static List<List<String>> separateNewickChildren (String newick) {
		List<String> children = extractNewickChildren(newick);
		List<String> leaves = new ArrayList<>();
		List<String> trees = new ArrayList<>();

		for (String child : children) {
			String trimmed = child.strip();
			if (trimmed.startsWith("(")) {
				// it's a subtree
				trees.add(trimmed);
			} else {
				// it's a leaf
				leaves.add(trimmed);
			}
		}
		List<List<String>> answer = new ArrayList<>(2);
		answer.add(leaves);
		answer.add(trees);
		return answer;
	}

	private static int currentIndex = 0;
	public static Dendrogram fromNewick (String newick) {
		int n = newick.length();
		int end = n - 1;
		while (newick.charAt(end) != ')') end--;
		currentIndex = 0;
		Dendrogram t = fromNewick(currentIndex, newick.substring(0, end + 1));
		t.generateIds();
		return t;
	}
	private static Dendrogram fromNewick (int index, String newick) {
		List<List<String>> newickChildren = separateNewickChildren(newick);
		List<String> newickLeaves = newickChildren.get(0);
		List<String> newickTrees = newickChildren.get(1);
		List<DendrogramEdge> edges = new ArrayList<>();

		// Process leaves
		for (String leaf : newickLeaves) {
			String[] destAndWeight = leaf.split(":");
			int len = destAndWeight.length;
			String label = String.join(":",
					Arrays.copyOfRange(destAndWeight, 0, len - 1)
			).strip();
			double weight = Double.parseDouble(destAndWeight[len - 1]);
			edges.add(new DendrogramEdge(weight, new Dendrogram(label)));
		}

		// Process subtrees
		for (String tree: newickTrees) {
			int lastSemicolon = tree.lastIndexOf(':');
			String dest = tree.substring(0, lastSemicolon);
			double weight = Double.parseDouble(tree.substring(lastSemicolon + 1));
			currentIndex++;
			Dendrogram subtree = fromNewick(currentIndex, dest);
			edges.add(new DendrogramEdge(weight, subtree));
		}

		return new Dendrogram("T" + index, edges);
	}

	

	private Dendrogram fromFile (InputStream in) {
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String newick = br.lines()
				.map(s -> s.replace("\n", "").strip()
				)
				.collect(Collectors.joining());
		return fromNewick(newick);
	}

	public String toNewick(){
		return children.isEmpty() ? "();" : toNewick(this) + ";";
	}

	private String toNewick(Dendrogram t){
		if (isLeaf(t)) return t.label;
		else {
			List<DendrogramEdge> currentChildren = t.children;
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

	public double getBranchLengthSum () {
		double sum = 0.0;
		Queue<Dendrogram> q = new LinkedList<>();
		q.add(this);
		while (!q.isEmpty()) {
			Dendrogram v = q.remove();
			for (DendrogramEdge e : v.children) {
				Dendrogram u = e.getDestination();
				double weight = e.getWeight();
				sum += weight;
				q.add(u);
			}
		}
		return sum;
	}

	private Dendrogram normalize (Dendrogram original, double totalSum) {
		Dendrogram norm = new Dendrogram(original.label);
		List<DendrogramEdge> children = new ArrayList<>();
		for (DendrogramEdge e : original.children) {
			Dendrogram u = e.getDestination();
			double weight = e.getWeight();
			Dendrogram w = normalize(u, totalSum);
			children.add(new DendrogramEdge(weight / totalSum, w));
		}
		norm.setChildren(children);
		return norm;
	}

	public Dendrogram normalize () {
		Dendrogram t =  this.normalize(this, this.getBranchLengthSum());
		t.generateIds();
		return t;
	}

	/**
	 * Joins a pair of nodes (u, v) to a new node x. Returns the new
	 * tree with x as the root
	 * @param name - Label of the new tree x
	 * @param distances - A pair (dux, dvx) with the distance between u and x, and the distance between v and x respectively
	 * @param nodes - The pair of trees (u, v) to be joined
	 * @return A new tree with x as the root
	 */
	public static Dendrogram join2 (String name, Dendrogram node1, double distance1, Dendrogram node2, double distance2) {
		return new Dendrogram(name, List.of(
				new DendrogramEdge(distance1, node1),
				new DendrogramEdge(distance2, node2)
		));
	}
}