package ngsep.graphs;

import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class StronglyConnectedComponents {
	
	private int clock = 0;
	
	public List<List<Integer>> computeStronglyConnectedComponents(List<Integer> idxs, List<List<Integer>> adjacencyList) {
		List<List<Integer>> components = new ArrayList<>();
		int n = adjacencyList.size();
		boolean[] visited = new boolean[n];
		int[] order = getOrderByReverseGraph(adjacencyList);
		for(int i = n-1; i >= 0; i--) {
			int u = order[i];
			if(!visited[u]) {
				List<Integer> answerCluster = new ArrayList<>();
				DFSToCluster(u, adjacencyList, visited, idxs, components, answerCluster);
				components.add(answerCluster);
			}
		}
		return components;
	}
	
	private void DFSToCluster(int u, List<List<Integer>> adjacencyList, boolean[] visited, List<Integer> idxs,
			List<List<Integer>> components, List<Integer> answerCluster) {
		// TODO Auto-generated method stub
		visited[u] = true;
		List<Integer> neighbors = adjacencyList.get(u);
		for(int v : neighbors) {
			if(!visited[v]) DFSToCluster(v, adjacencyList, visited, idxs, components, answerCluster);
		}
		int idx = idxs.get(u);
		answerCluster.add(idx);
	}
	
	private int[] getOrderByReverseGraph(List<List<Integer>> adjacencyList) {
		// TODO Auto-generated method stub
		int n = adjacencyList.size();
		List<List<Integer>> reverseGraph = new ArrayList<>();
		int[] order = new int[n];
		for(int l = 0; l < n; l++) reverseGraph.add(new ArrayList<>());
		for(int i = 0; i < n; i++) {
			List<Integer> neighbors = adjacencyList.get(i);
			for(int j:neighbors) {
				reverseGraph.get(j).add(i);
			}
		}
		/**for(List<Integer> i : reverseGraph) {
			for(int j : i) {
				System.out.print(j+1 + " ");
			}
			System.out.println();
		}**/
		boolean[] reverseVisited = new boolean[n];
		for(int s = 0; s < n; s++) {
			if(!reverseVisited[s]) DSFToUpdateClock(reverseGraph, s, reverseVisited, order);
		}
		return order;
	}
	
	private void DSFToUpdateClock(List<List<Integer>> reverseGraph, int s, boolean[] visited, int[] order) {
		// TODO Auto-generated method stub
		visited[s] = true;
		List<Integer> neighbors = reverseGraph.get(s);
		for(int t : neighbors) {
			if(!visited[t]) DSFToUpdateClock(reverseGraph, t, visited, order);
		}
		order[clock] = s;
		clock++;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Scanner scanner = new Scanner(System.in);
        int n = scanner.nextInt();
        int m = scanner.nextInt();
        List<Integer> idxs = new ArrayList<>();
        List<List<Integer>> adj = new ArrayList<>();
        for (int i = 0; i < n; i++) {
        	idxs.add(i);
        	adj.add(new ArrayList<>());
        }
        for (int i = 0; i < m; i++) {
            int x, y;
            x = scanner.nextInt();
            y = scanner.nextInt();
            adj.get(x-1).add(y - 1);
        }        
        List<List<Integer>> components = new StronglyConnectedComponents().computeStronglyConnectedComponents(idxs, adj);
       //System.out.println(components.size());
        for(List<Integer> c : components) {
        	System.out.print("{");
        	for(int i = 0; i < c.size(); i++) {
        		if(i==c.size()-1) System.out.println((c.get(i) + 1) + "}");
        		else System.out.print((c.get(i)+ 1) + ",");
        	}
        }
	}

}
