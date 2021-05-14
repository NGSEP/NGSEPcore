package ngsep.clustering.dendrogram;

import ngsep.clustering.Pair;

import java.io.*;
import java.util.*;
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

	public List<List<Pair<Integer, Double>>> toAdjacencyList () {
		int n = this.size;
		List<List<Pair<Integer, Double>>> adj = Stream.generate(ArrayList<Pair<Integer, Double>>::new)
				.limit(n)
				.collect(Collectors.toList());
		Queue<Dendrogram> q = new LinkedList<>();
		q.add(this);
		while (!q.isEmpty()) {
			Dendrogram v = q.remove();
			for (DendrogramEdge e : v.children) {
				Dendrogram u = e.getDestination();
				adj.get(v.getId()).add(new Pair<>(
						u.getId(), e.getWeight()
				));
				q.add(u);
			}
		}
		return adj;
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

	private static Pair<List<String>, List<String>> separateNewickChildren (String newick) {
		List<String> children = extractNewickChildren(newick);
		List<String> leaves = new ArrayList<>();
		List<String> trees = new ArrayList<>();

		for (String child : children) {
			String trimmed = child.strip();
			if (child.startsWith("(")) {
				// it's a subtree
				trees.add(trimmed);
			} else {
				// it's a leaf
				leaves.add(trimmed);
			}
		}
		return new Pair<>(leaves, trees);
	}

	private static Pair<Integer, Dendrogram> fromNewick (int index, String newick) {
		Pair<List<String>, List<String>> newickChildren = separateNewickChildren(newick);
		List<String> newickLeaves = newickChildren.first;
		List<String> newickTrees = newickChildren.second;
		List<DendrogramEdge> edges = new ArrayList<>();

		// Process leaves
		for (String leaf : newickLeaves) {
			String[] destAndWeight = leaf.split(":");
			String label = destAndWeight[0].strip();
			double weight = Double.parseDouble(destAndWeight[1]);
			edges.add(new DendrogramEdge(weight, new Dendrogram(label)));
		}

		// Process subtrees
		int currIndex = index;
		for (String tree: newickTrees) {
			int lastSemicolon = tree.lastIndexOf(':');
			String dest = tree.substring(0, lastSemicolon);
			double weight = Double.parseDouble(
					tree.substring(lastSemicolon + 1)
			);
			Pair<Integer, Dendrogram> indexAndTree = fromNewick(currIndex + 1, dest);
			currIndex = indexAndTree.first;
			Dendrogram subtree = indexAndTree.second;
			edges.add(new DendrogramEdge(weight, subtree));
		}

		return new Pair<>(
				currIndex,
				new Dendrogram("T" + index, edges)
		);
	}

	public static Dendrogram fromNewick (String newick) {
		int n = newick.length();
		int end = n - 1;
		while (newick.charAt(end) != ')') end--;
		Dendrogram t = fromNewick(0, newick.substring(0, end + 1)).second;
		t.generateIds();
		return t;
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

	public static void main(String[] args) throws IOException{
		String s = "(x:1,x:2,(x:3,(x:5,x:6):4):10,x:7,((y:1,y:3):5,y:2):14);";
		Dendrogram t0 = Dendrogram.fromNewick(s);
		System.out.println(t0.toAdjacencyList());
		t0.printTree(System.out);
		String file = args[0];
		Dendrogram t1 = new Dendrogram(new FileInputStream(file));
		t1.printTree(System.out);
	}
}