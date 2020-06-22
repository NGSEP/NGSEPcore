/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;
import java.util.Map.Entry;

/**
 * 
 * @author Juan Camilo Bojaca
 *
 */
public class LayoutBuilderMetricMSTChristofides implements LayoutBuilder {

	@Override
	public void findPaths(AssemblyGraph graph) {

		Collection<Edge> mst = MinimumSpanningTreeGenerator.KRUSKAL.generate(graph);
		System.out.println("Generated minimum spanning tree. Number of edges: "+mst.size());
		Set<Integer> oddDegree = selectVerticesWithOddDegree(graph.getVertices().size(), mst);
		System.out.println("Number of nodes with odd degree in MST: "+oddDegree.size());
		Collection<Edge> match = (new Blossom()).generate(graph, oddDegree);
		System.out.println("Calculated min cost match in subgraph. Edges: "+match.size());
		mst.addAll(match);
			
		LinkedList<Integer> path = (new EulerianCircuitGenerator()).generate(mst);
		System.out.println("Calculated eulerian path. Size: "+path.size());
		
		// to start with a valid sequence
		if (path.get(0) / 2 == path.get(path.size() - 1) / 2) {
			path.addFirst(path.removeLast());
		}
		removeDuplicateNodes(graph, path);
		System.out.println("Removed duplicates. Size: "+path.size());
		addPathsToTheGraph(graph, path);
	}

	/**
	 * Take the path split it and adds to the graph
	 * @param graph
	 * @param path
	 */
	private void addPathsToTheGraph(AssemblyGraph graph, LinkedList<Integer> path) {
		// Create the map of edges in the graph
		Map<Integer, Map<Integer, AssemblyEdge>> m = new HashMap<>();
		for (AssemblyEdge edge : graph.getEdges()) {
			int v = edge.getVertex1().getSequenceIndex(), w = edge.getVertex2().getSequenceIndex();
			m.putIfAbsent(v, new HashMap<>());
			m.putIfAbsent(w, new HashMap<>());
			m.get(v).put(w, edge);
			m.get(w).put(v, edge);
		}

		// Cut the path if there is not an real edge
		List<AssemblyEdge> aux = new LinkedList<>();
		for (int i = 0; i < path.size() - 1; i++) {
			int v = path.get(i), w = path.get(i + 1);
			if (m.containsKey(v) && m.get(v).containsKey(w)) {
				aux.add(m.get(v).get(w));
			} else {
				graph.addPath(aux);
				aux.clear();
			}
		}
	}

	/**
	 * Removes the duplicated nodes of the path
	 *
	 * @param G    the graph
	 * @param path the path
	 */
	private void removeDuplicateNodes(AssemblyGraph graph, LinkedList<Integer> path) {
		Integer[] p = (Integer[]) path.toArray();

		// The list of all the appears of the vertex
		Map<Integer, List<Integer>> added = new HashMap<>();
		for (int i = 1; i < p.length - 1; i++) {
			added.putIfAbsent(p[i], new LinkedList<>());
			added.get(p[i]).add(i);
		}
		List<AssemblyVertex> vertices = graph.getVertices();
		// Keep the vertex with the less edge weight
		Set<Integer> toRemove = new HashSet<>();
		for (Entry<Integer, List<Integer>> entry : added.entrySet())
			if (entry.getValue().size() > 2) {
				int min = Integer.MAX_VALUE;
				AssemblyVertex vertex = vertices.get(entry.getKey());
				int indexMin = -1;
				for (int i : entry.getValue()) {
					AssemblyVertex previous = vertices.get(p[i-1]);
					AssemblyEdge e1 = graph.getEdge(vertex, previous);
					AssemblyVertex next = vertices.get(p[i+1]);
					AssemblyEdge e2 = graph.getEdge(vertex, next);
					int nextMin = Math.min(e1.getCost(), e2.getCost());
					
					if (min > nextMin) {
						min = nextMin;
						indexMin = i;
					}
				}
				for (int i : entry.getValue()) {
					if (i != indexMin) toRemove.add(i);
				}
			}

		// remove in the list
		Iterator<Integer> iter = path.iterator();
		while (iter.hasNext()) {
			Integer integer = (Integer) iter.next();
			if (toRemove.contains(integer))
				iter.remove();

		}
	}
	/**
	 * return the array marking the odd degree vertex in the the minimum spanning
	 * tree
	 *
	 * @param graph the graph
	 * @param mst   the minimum spanning tree
	 */
	private static Set<Integer> selectVerticesWithOddDegree(int n, Collection<Edge> mst) {
		int[] count = new int[n];
		for (Edge edge : mst) {
			count[edge.getV()]++;
			count[edge.getW()]++;
		}
		Set<Integer> ans = new HashSet<Integer>();
		for (int i = 0; i < n; i++) {
			if (count[i]%2==1) {
				ans.add(i);
			}
		}
		return ans;
	}

}

@FunctionalInterface
interface MinimumSpanningTreeGenerator {
	public Collection<Edge> generate(AssemblyGraph graph);

	public static MinimumSpanningTreeGenerator KRUSKAL = (AssemblyGraph graph) -> {
		List<Edge> ans = new ArrayList<>();
		// Add all the sets
		List<AssemblyVertex> vertices = graph.getVertices();
		HashMap<Integer, DisjointSet> sets = new HashMap<>();
		for (int v = 0; v < vertices.size(); v++)
			sets.put(v, new DisjointSet());

		// Sort the edges
		PriorityQueue<Edge> edges = new PriorityQueue<>((Edge a, Edge b) -> a.getWeigth() - b.getWeigth());
		for (int v = 1; v < vertices.size(); v++) {
			AssemblyVertex vertex = vertices.get(v);
			List<AssemblyEdge> edgesG = graph.getEdges(vertex);
			for(AssemblyEdge edge:edgesG) {
				AssemblyVertex v2 = edge.getConnectingVertex(vertex);
				//TODO: Make this more efficient
				int w = vertices.indexOf(v2);
				if (w<v) edges.add(new Edge(v, w, edge.getCost()));
			}
		}	

		// Fuse the sets
		for (Edge edge : edges) {
			DisjointSet a = sets.get(edge.getV()).find(), b = sets.get(edge.getW()).find();
			if (a == b) continue;
			ans.add(edge);
			a.union(b);
			if (ans.size() == sets.size())
				break;
		}

		return ans;
	};
}

/**
 * A element of the disjoint set
 *
 * @author jc.bojaca
 */
class DisjointSet {
	private DisjointSet parent;
	private int rank;

	DisjointSet() {
		parent = this;
		rank = 0;
	}

	/**
	 *
	 * @return the set of the value
	 */
	DisjointSet find() {
		DisjointSet next, x = this;
		while (x.parent != x) {
			next = x.parent;
			x.parent = next.parent;
			x = next;
		}
		return x;
	}

	/**
	 *
	 * @param y the set to merge
	 */
	void union(DisjointSet y) {
		DisjointSet xRoot = find();
		DisjointSet yRoot = y.find();

		// x and y are already in the same set
		if (xRoot == yRoot)
			return;

		// x and y are not in same set, so we merge them
		if (xRoot.rank < yRoot.rank) {
			// swap xRoot and yRoot
			DisjointSet t = xRoot;
			xRoot = yRoot;
			yRoot = t;
		}

		// merge yRoot into xRoot
		yRoot.parent = xRoot;
		if (xRoot.rank == yRoot.rank)
			xRoot.rank++;
	}
}



class Blossom {
	Map<Node, Node> matches = new HashMap<>();
	Map<Node, Tree> forest = new HashMap<>();
	
	Set<Node> G;

	public Collection<Edge> generate(AssemblyGraph graph, Set<Integer> selectedVertices) {
		List<Edge> ans = new ArrayList<>();
		//Build internal subgraph
		G = getNodes(graph, selectedVertices);
		System.out.println("Built subgraph with: "+G.size()+" vertices");
		initializeDualVariables();
		System.out.println("Initialized dual variables");
		refreshEdges();
		System.out.println("refreshed edges");
		Collection<Node> unmatched = primalUpdates();
		int numUnmatched = unmatched.size();
		System.out.println("Matching performed. Number of unmatched: "+unmatched.size());
		while (!unmatched.isEmpty()) {
			updateSlackValues();
			System.out.println("Performed update of slack values");
			unmatched = primalUpdates();
			System.out.println("Matching performed. Number of unmatched: "+unmatched.size());
			if(numUnmatched==unmatched.size()) throw new RuntimeException("Number of unmatched edges could not be reduced");
			numUnmatched = unmatched.size();
		}

		adjustMatches();
		List<AssemblyVertex> vertices = graph.getVertices();

		for (Entry<Node, Node> entry : matches.entrySet()) {
			int id1 = entry.getKey().id;
			int id2 = entry.getValue().id;
			if (id1 < id2) {
				AssemblyVertex v1 = vertices.get(id1);
				AssemblyVertex v2 = vertices.get(id2);
				AssemblyEdge edge = graph.getEdge(v1, v2);
				ans.add(new Edge(entry.getKey().id, entry.getValue().id, edge.getCost()));
			}
		}	
		return ans;
	}

	/**
	 * calculate the real match recursively on blossom nodes
	 */
	private void adjustMatches() {
		Stack<Node> pendingVertices = new Stack<>();
		pendingVertices.addAll(G);
		while (!pendingVertices.isEmpty()) {
			Node node = pendingVertices.pop();
			if (node.isBlossom()) {
				Node match = matches.get(node);
				Node real = node.edgeSuport.get(match);
				blossomMatchAdjust(node);
				matches.remove(node);
				matches.put(match, real);
				matches.put(real, match);

				pendingVertices.addAll(node.contracted);
			}
		}
	}

	/**
	 * Corrects the match of the contracted list of a blossom node
	 *
	 * @param node the node
	 */
	private void blossomMatchAdjust(Node node) {
		LinkedList<Node> queue = node.contracted;
		int ind = queue.indexOf(node.edgeSuport.get(matches.get(node)));
		if (ind != 0) {
			// transform in the base
			for (int i = 0; i < ind; i++)
				queue.addLast(queue.pollFirst());

			// remove the current match
			for (Node innerNode : queue)
				matches.remove(innerNode);

			Iterator<Node> iter = queue.iterator();
			iter.next();// base
			for (int i = 1; i < queue.size() - 1; i += 2) {
				Node a = iter.next(), b = iter.next();
				matches.put(a, b);
				matches.put(b, a);
			}
		}
	}

	/**
	 * make all the primary updates over the current graph
	 *
	 * @param unmatched
	 * @return
	 */
	private void updateSlackValues() {
		double value1 = calculateMinimumNegativeBlossomVertexExpand();
		double value2 = calculateMinimumGrow();
		double value3 = calculateMinimumPositiveVertexShrink();
		double delta = Math.min(value1, Math.min(value2, value3));
		for (Tree t : forest.values())
			t.node.y += ((t.positive) ? 1 : -1) * delta;
		refreshEdges();
	}

	/**
	 * @return the minimum delta for Expand operations u with u.isblosom and y(u)=0
	 */
	private double calculateMinimumNegativeBlossomVertexExpand() {
		double delta = Double.POSITIVE_INFINITY;
		for (Node v : G) {
			Tree vTree = forest.get(v);
			if(vTree!=null && !vTree.positive && vTree.node.isBlossom()) {
				delta = Math.min(delta, v.y);
			}
		}		
		return delta;
	}

	/**
	 * @return the maximum delta for argument and shrink operations (u,v) with u(+)
	 *         v(+)
	 */
	private double calculateMinimumPositiveVertexShrink() {
		double delta = Double.POSITIVE_INFINITY;
		for (Node v : G) {
			Tree vTree = forest.get(v);
			if (vTree!=null && vTree.positive) {
				for (Entry<Node, Double> edge : v.edges.entrySet())
					if (forest.containsKey(edge.getKey()) && forest.get(edge.getKey()).positive)
						delta = Math.min(delta, (edge.getValue() - v.y - edge.getKey().y) / (double) 2);
			}
			
		}
			
		return delta;
	}

	/**
	 *
	 * @return the minimum delta for grow operations (4a) (u,v) with u(+) v(0)
	 */
	private double calculateMinimumGrow() {
		double delta = Double.POSITIVE_INFINITY;
		for (Node v : G) {
			Tree vTree = forest.get(v);
			if (vTree!=null) {
				for (Entry<Node, Double> edge : v.edges.entrySet())
					if (forest.containsKey(edge.getKey()) && forest.get(edge.getKey()).positive)
						delta = Math.min(delta, edge.getValue() - v.y - edge.getKey().y);
			}
				
		}
			
		return delta;
	}

	/**
	 * make all the primary updates over the current graph
	 *
	 * @param unmatched
	 * @return
	 */
	private Collection<Node> primalUpdates() {
		LinkedList<Node> unmatchedNodes = getUnmatchedNodes();
		Stack<Node> pendingVertices = new Stack<>();
		pendingVertices.addAll(unmatchedNodes);
		forest.clear();
		for (Node node : unmatchedNodes) forest.put(node, new Tree(node, null));
		
		while (!pendingVertices.isEmpty()) {
			Node v = pendingVertices.pop();
			if (!G.contains(v)) continue;
			//if(!forest.containsKey(v)) continue;

			for (Node w : v.tights) {
				if (matches.get(w) == v) continue;

				if (!forest.containsKey(w)) {
					if (w.isBlossom() && w.y == 0) {
						Expand(w);
						pendingVertices.add(v);// is added again to Grow latter
						break;
					} else {
						// Grow tree
						Tree vTree = forest.get(v);
						Tree wTree = new Tree(w, vTree);
						forest.put(w, wTree);
						Node match = matches.get(w);
						Tree matchTree = new Tree(match, wTree);
						forest.put(match, matchTree);
						pendingVertices.add(match);
					}		
				} else if (haveSameRoot(v, w)) {
					addToMatching(v, w);
					unmatchedNodes = getUnmatchedNodes();
					pendingVertices.clear();
					pendingVertices.addAll(unmatchedNodes);
					forest.clear();
					for (Node node : unmatchedNodes) forest.put(node, new Tree(node, null));
					break;
				} else {
					pendingVertices.add(Shrink(v, w));
					if (!forest.containsKey(v)) break;
				}
			}
		}
		return getUnmatchedNodes();
	}
	



	/**
	 * Primary operation Shrink: contracts a blossom to a single node argument it
	 * takes a edge (nodeA,nodeB) and contacts using as base the common parent
	 *
	 * @param nodeA node leaf in the tree
	 * @param nodeB node leaf in the tree
	 * @return the contracted node
	 */
	private Node Shrink(Node nodeA, Node nodeB) {
		LinkedList<Node> blossomPath = new LinkedList<>();
		Tree parent = getBlossom(forest.get(nodeA), forest.get(nodeB), blossomPath).parent;
		Set<Node> blossomNodes = new HashSet<>(blossomPath);

		Set<Node> tightEdges = new HashSet<>();
		for (Node node : blossomNodes)
			for (Node adjNode : node.tights)
				if (!blossomNodes.contains(adjNode))
					tightEdges.add(adjNode);

		Map<Node, Double> weigths = new HashMap<>();
		Map<Node, Node> suport = new HashMap<>();
		for (Node node : blossomPath)
			for (Entry<Node, Double> edge : node.edges.entrySet())
				if (!blossomNodes.contains(edge.getKey()) && !tightEdges.contains(edge.getKey())) {
					Node adjNode = edge.getKey();
					double weigth = edge.getValue() - node.y;
					if (!weigths.containsKey(adjNode) || weigths.get(adjNode) > weigth) {
						weigths.put(adjNode, weigth);
						suport.put(adjNode, node);
					}
				}

		Node blossomNode = new Node(-1);
		blossomNode.contracted = blossomPath;
		blossomNode.tights = tightEdges;
		blossomNode.edges = weigths;
		blossomNode.edgeSuport = suport;

		updateBlossomInGraph(blossomNode);
		updateBlossomInForest(parent, blossomNode);
		return blossomNode;
	}

	/**
	 *
	 * @param base
	 * @param blossomNodes
	 * @param blossomNode
	 */
	private void updateBlossomInForest(Tree parent, Node blossomNode) {
		for (Node nodeToRemove : blossomNode.contracted)
			forest.remove(nodeToRemove);
		forest.put(blossomNode, new Tree(blossomNode, parent));
	}

	/**
	 * make the corrections in the graph to add a blossom node
	 *
	 * @param blossomNode the new node
	 */
	private void updateBlossomInGraph(Node blossomNode) {
		G.removeAll(blossomNode.contracted);
		for (Node node : G) {
			for (Node nodeToRemove : blossomNode.contracted)
				node.edges.remove(nodeToRemove);

			if (blossomNode.edges.containsKey(node))
				node.edges.put(blossomNode, blossomNode.edges.get(node));

			if (blossomNode.tights.contains(node)) {
				node.tights.removeAll(blossomNode.contracted);
				node.tights.add(blossomNode);
			}
		}
		G.add(blossomNode);
	}

	/**
	 *
	 * @param a    tree
	 * @param b    tree
	 * @param path list to add the path of the blossom sequence starting by the base
	 * @return the base of the blossom
	 */
	private Tree getBlossom(Tree a, Tree b, LinkedList<Node> path) {
		if (a.heigth < b.heigth) {
			Tree t = a;
			a = b;
			b = t;
		}
		while (a.heigth > b.heigth) {
			path.addFirst(a.node);
			a = a.parent;
		}
		while (a != b) {
			path.addFirst(a.node);
			path.addLast(b.node);
			a = a.parent;
			b = b.parent;
		}
		path.addFirst(a.node);
		return a;
	}

	private boolean haveSameRoot(Node v, Node w) {
		Tree tW = forest.get(w);
		Tree tV = forest.get(v);
		if(tV==null) throw new RuntimeException("v is not in the forest");
		if(tW==null) throw new RuntimeException("W is not in the forest");
		return tW!=null && tV!=null && tW.root() != tV.root();
	}

	/**
	 * Changes the matching to include a direct matching between the given nodes.
	 * Removes matching edges between A and B and adds unmatching edges
	 *
	 * @param nodeA node leaf in the tree
	 * @param nodeB node leaf in the tree
	 */
	private void addToMatching(Node nodeA, Node nodeB) {
		// remove the current match
		Tree a = forest.get(nodeA);
		Tree b = forest.get(nodeB);
		for (Tree tree : new Tree[] { a, b }) {
			while (tree != null && tree.parent != null) {
				matches.remove(tree.node);
				matches.remove(tree.parent.node);
				tree = tree.parent.parent;
			}
		}

		// add the new ones
		for (Tree tree : new Tree[] { a.parent, b.parent }) {
			while (tree != null && tree.parent != null) {
				matches.put(tree.node, tree.parent.node);
				matches.put(tree.parent.node, tree.node);
				tree = tree.parent.parent;
			}
		}
		matches.put(nodeA, nodeB);
		matches.put(nodeB, nodeA);
	}


	/**
	 * Primary operation Expand
	 *
	 * @param v the node that is a blossom and needs to be expanded (is not added to
	 *          the tree)
	 */
	private void Expand(Node v) {
		// remove in the graph
		G.remove(v);
		for (Node node : G) {
			node.tights.remove(v);
			node.edges.remove(v);
		}

		// complete the symmetry
		for (Node inderNode : v.contracted) {
			for (Node adyacentNode : inderNode.tights)
				adyacentNode.tights.add(inderNode);

			for (Entry<Node, Double> edge : inderNode.edges.entrySet())
				edge.getKey().edges.put(inderNode, edge.getValue());
		}

		// adjust the Matches
		blossomMatchAdjust(v);
	}



	/**
	 *
	 * @return the collection of nodes without a match
	 */
	private LinkedList<Node> getUnmatchedNodes() {
		LinkedList<Node> ans = new LinkedList<>();
		for (Node node : G)
			if (!matches.containsKey(node)) ans.add(node);
		return ans;
	}

	/**
	 * Calculate the slack of all the edges to add new tight edges
	 */
	private void refreshEdges() {
		for (Node node : G)
			node.refreshTightEdges();
	}

	/**
	 * Initialize the dual variables of the graph
	 */
	private void initializeDualVariables() {
		for (Node node : G)
			node.y = node.minEdge() / (double) 2;
	}

	/**
	 *
	 * @param graph the graph representation
	 * @return Set<Node> the graph as a set of Node objects
	 */
	private Set<Node> getNodes(AssemblyGraph graph, Set<Integer> selectedVertices) {
		Map<Integer,Node> nodes = new HashMap<Integer, Node>();
		List<AssemblyVertex> vertices = graph.getVertices();
		for (int v:selectedVertices) {
			Node a = new Node(v);
			nodes.put(v,a);
		}
		for (int v:selectedVertices) {
			Node a = nodes.get(v);
			AssemblyVertex vertex = vertices.get(v);
			List<AssemblyEdge> edges = graph.getEdges(vertex);
			for(AssemblyEdge edge:edges) {
				AssemblyVertex v2 = edge.getConnectingVertex(vertex);
				//TODO: Make this more efficient
				int w = vertices.indexOf(v2);
				Node n2 = nodes.get(w);
				if (n2!=null) a.edges.put(n2, (double)edge.getCost());
			}
		}
		return new HashSet<>(nodes.values());
	}

}

/**
 * Alternating Tree
 *
 * @author jc.bojaca
 *
 */
class Tree {
	Node node;
	Tree parent;
	boolean positive;
	int heigth = 0;

	public Tree(Node node, Tree parent) {
		this.node = node;
		this.parent = parent;
		this.positive = (parent == null) || !parent.positive;
		if (parent != null) heigth = parent.heigth + 1;
	}

	/**
	 *
	 * @return the root of the tree
	 */
	Tree root() {
		Tree x = this;
		while (x.parent != null) x = x.parent;
		return x;
	}

}

/**
 * A node of the contracted graph could represent a node or a contracted blossom
 *
 * @author jc.bojaca
 *
 */
class Node {
	LinkedList<Node> contracted = new LinkedList<>();
	Map<Node, Node> edgeSuport = new HashMap<>();
	Set<Node> tights = new HashSet<>();
	Map<Node, Double> edges = new HashMap<>();
	double y;
	int id;

	Node(int id) {
		this.id = id;
	}

	boolean isBlossom() {
		return !contracted.isEmpty();
	}

	/**
	 * @return the minimum weight in the edges
	 */
	double minEdge() {
		double min = Double.POSITIVE_INFINITY;
		for (double edgeW : edges.values())
			if (edgeW < min) min = edgeW;
		return min;
	}

	/**
	 * Calculate the slack for every edge and remove it if is zero
	 */
	void refreshTightEdges() {
		Iterator<Entry<Node, Double>> iter = edges.entrySet().iterator();
		while (iter.hasNext()) {
			Entry<Node, Double> entry = (Entry<Node, Double>) iter.next();
			double slack = entry.getValue() - y - entry.getKey().y;
			if (slack == 0) {
				tights.add(entry.getKey());
				iter.remove();
			}
		}
	}

	@Override
	public String toString() {
		return id + "" + Arrays.toString(edges.keySet().stream().mapToInt((n) -> n.id).toArray());
	}

}

class EulerianCircuitGenerator {
	LinkedList<Integer> generate(Collection<Edge> edges) {
		Iterator<Integer> iter;
		// create the map of pending task
		Map<Integer, List<Integer>> map = new HashMap<>();
		for (Edge edge : edges) {
			map.putIfAbsent(edge.getV(), new LinkedList<>());
			map.get(edge.getV()).add(edge.getW());
			map.putIfAbsent(edge.getW(), new LinkedList<>());
			map.get(edge.getW()).add(edge.getV());
		}
		int t = 1;

		// first paths
		ListNode root = new ListNode(0);
		ListNode act = root;
		iter = map.get(0).iterator();
		Integer v = 0, a = iter.next();
		iter.remove();
		map.get(a).remove(v);
		do {
			act.next = new ListNode(a);
			act = act.next;
			t++;
			iter = map.get(a).iterator();
			int b = iter.next();
			iter.remove();
			map.get(b).remove(a);
			a = b;
		} while (a != v);

		// update
		while (t < edges.size()) {
			act = root;
			while (act != null) {
				ListNode next = act.next;
				v = a = act.value;
				while (!map.get(a).isEmpty()) {
					iter = map.get(a).iterator();
					int b = iter.next();
					iter.remove();
					map.get(b).remove(a);
					a = b;

					act.next = new ListNode(a);
					act = act.next;
					t++;
				}
				act.next = next;
				act = act.next;
			}
		}

		// Transform to list
		act = root;
		LinkedList<Integer> ans = new LinkedList<>();
		while (act != null) {
			ans.add(act.value);
			act = act.next;
		}
		return ans;
	}

	static long key(int a, int b) {
		return ((long) a << 32) | b;
	}
}

class ListNode {
	ListNode next;
	int value;

	public ListNode(int value) {
		this.value = value;
	}
}

class Edge {
	/**
	 * Vertexes
	 */
	private int w, v;
	/**
	 * Weight of the edge
	 */
	private int weigth;

	public Edge(int w, int v, int weigth) {
		super();
		this.w = w;
		this.v = v;
		this.weigth = weigth;
	}

	public int getW() {
		return w;
	}

	public void setW(int w) {
		this.w = w;
	}

	public int getV() {
		return v;
	}

	public void setV(int v) {
		this.v = v;
	}

	public int getWeigth() {
		return weigth;
	}

	public void setWeigth(int weigth) {
		this.weigth = weigth;
	}

	@Override
	public String toString() {
		return "[" + w + " , " + v + "]";
	}
}
