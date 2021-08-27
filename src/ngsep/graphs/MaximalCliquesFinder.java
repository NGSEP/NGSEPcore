package ngsep.graphs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class MaximalCliquesFinder {
	
	private static ArrayList<Node> generateInputGraph(List<Integer> idxs, boolean[][] edges) {
		 ArrayList<Node> inputGraph = new ArrayList<>();
		for (int x = 0; x < idxs.size(); x++) {
			inputGraph.add(new Node(idxs.get(x)));
		}
		for(int i = 0; i < edges.length; i++) {
			for(int j = i; j < edges[0].length; j++) {
				if(edges[i][j]) {
					inputGraph.get(i).addNeighbor(inputGraph.get(j));
				}
			}
		}
		return inputGraph;
	}
	private static ArrayList<Node> getIntersections(ArrayList<Node> nodes1, ArrayList<Node> nodes2) {
		ArrayList<Node> intersections = new ArrayList<>(nodes1);
		intersections.retainAll(nodes2);
		return intersections;
	}
	private static ArrayList<Node> getUnions(ArrayList<Node> nodes1, List<Node> nodes2) {
		ArrayList<Node> intersections = new ArrayList<>(nodes1);
		intersections.addAll(nodes2);
		return intersections;
	}
	private static ArrayList<Node> removeNodeNeighbors(ArrayList<Node> nodeList, Node v) { 
		ArrayList<Node> nodes = new ArrayList<Node>(nodeList); 
		nodes.removeAll(v.getNeighbors()); 
		return nodes; 
	} 
	private static Node getPivotNode(ArrayList<Node> nodeList) {
		Collections.sort(nodeList);
		Node pivot = nodeList.get(nodeList.size()-1);
		return pivot;
	}
	private static ArrayList<Integer> getNodeIdxs(ArrayList<Node> nodeList) {
		ArrayList<Integer> idxs = new ArrayList<>();
		for(Node v:nodeList) {
			idxs.add(v.getIndex());
		}
		return idxs;
	}
	private static void pivotBronKerbosch(ArrayList<Node> R, ArrayList<Node> P, ArrayList<Node> X,
			ArrayList<ArrayList<Integer>> maxCliques) {
		if(P.isEmpty() && X.isEmpty()) {
			ArrayList<Node> Rtemp = new ArrayList<>(R);
			maxCliques.add(getNodeIdxs(Rtemp));
			return;
		}
		ArrayList<Node> Ptemp = new ArrayList<Node>(P); 
		Node pivot = getPivotNode(getUnions(P, X));
		P = removeNodeNeighbors(P, pivot);
		for(Node v:P) {
			R.add(v);
			pivotBronKerbosch(R, getIntersections(Ptemp, v.getNeighbors()),
					getIntersections(X, v.getNeighbors()), maxCliques);
			R.remove(v);
			Ptemp.remove(v);
			X.add(v);
		}
	}
	public static ArrayList<ArrayList<Integer>> callMaxCliqueFinder(List<Integer> idxs, boolean[][] edges){
		ArrayList<Node> inputGraph = generateInputGraph(idxs, edges);
		ArrayList<ArrayList<Integer>> result = new ArrayList<>();
		pivotBronKerbosch(new ArrayList<>(), inputGraph, new ArrayList<>(), result);
		return result;
	}
/**	
	public static void main(String[] args) {
		List<List<Integer>> e = Arrays.asList(Arrays.asList(2,5),
				Arrays.asList(1,3,5),
				Arrays.asList(2,4),
				Arrays.asList(3,5,6),
				Arrays.asList(1,2,4),
				Arrays.asList(4));	
		ArrayList<Integer> idxs = new ArrayList<>(Arrays.asList(1,2,3,4,5,6));
		boolean[][] edges = {{false, true, false, false, true,false},
							{true, false, true, false, true,false},
							{false, true, false, true, false,false},
							{false, false, true, false, true,true},
							{true, true, false, true, false,false},
							{false, false, false, true, false,false}};
		ArrayList<ArrayList<Integer>> cliques = callMaxCliqueFinder(idxs, edges);
		for(ArrayList<Integer> r:cliques) {
			System.out.print(" Maximal Clique: ");
			for(int v:r) {
				System.out.print(v + " ");
			}
			System.out.println();
		}
	}
**/
}
class Node implements Comparable <Node>{
	int idx;
	int degree;
	private ArrayList<Node> neighbors;
	
	public Node(int idx){
		super();
		this.idx = idx;
		neighbors = new ArrayList<>();
	}
	
	@Override
	public int compareTo(Node o2) {
		int cmp = 0; 
		if (this.degree < o2.degree) return -1;
        if (this.degree > o2.degree) return 1;
        return cmp;
	}
	
	public int getIndex() {
		return idx;
	}

	public int getDegree() {
		return degree;
	}
	
	public  ArrayList<Node> getNeighbors(){
		return neighbors;
	}
	
	public void addNeighbor(Node n) {
		this.neighbors.add(n);
		this.degree++;
		if (!n.getNeighbors().contains(n)) { 
            n.getNeighbors().add(this); 
            n.degree++; 
        } 
	}
}
