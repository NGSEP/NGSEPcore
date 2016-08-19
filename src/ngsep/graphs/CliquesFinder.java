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
package ngsep.graphs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class CliquesFinder {
	public static List<List<Integer>> findCliques(boolean[][] consistencyMatrix) {
		int n = consistencyMatrix.length;
		//Choose next largest
		boolean [] visited = new boolean [n];
		Arrays.fill(visited, false);
		List<List<Integer>> answer = new ArrayList<List<Integer>>();
		while (true) {
			List<Integer> nextCluster = findNextLargestClique(consistencyMatrix,visited);
			if(nextCluster.size()<=1) break;
			answer.add(nextCluster);
			for(int i:nextCluster) visited[i] = true;
		}
		return answer;
	}

	private static List<Integer> findNextLargestClique(boolean[][] consistencyMatrix, boolean[] visited) {
		int n = consistencyMatrix.length;
		//Build tree with edges for nodes to process
		Map<Integer,List<Integer>> edges = new TreeMap<Integer, List<Integer>>();
		for(int i=0;i<n;i++) {
			edges.put(i, new ArrayList<Integer>());
		}
		for(int i=0;i<n;i++) {
			if(!visited[i]) {
				edges.get(i).add(i);
				for(int j=i+1;j<n;j++) {
					if(!visited[j] && consistencyMatrix[i][j]) {
						edges.get(i).add(j);
						edges.get(j).add(i);
					}
				}
			}
		}
		//Sort nodes by degree
		List<NodeIndexPlusDegree> nodesToProcess = new ArrayList<NodeIndexPlusDegree>();
		for(int i=0;i<n;i++) {
			if(!visited[i]) {
				nodesToProcess.add(new NodeIndexPlusDegree(i, edges.get(i).size()));
			}
		}
		Collections.sort(nodesToProcess);
		Collections.reverse(nodesToProcess);
		List<Integer> maxCluster = new ArrayList<Integer>();
		
		for(NodeIndexPlusDegree next:nodesToProcess) {
			if(next.getDegree()<maxCluster.size()) break;
			List<Integer> cluster = findLargestClique(next.getIndex(),consistencyMatrix,edges);
			if(cluster.size()>maxCluster.size())maxCluster=cluster;
		}
		return maxCluster;
	}

	private static List<Integer> findLargestClique(int pivot,boolean[][] consistencyMatrix, Map<Integer,List<Integer>> edges) {
		List<Integer> candidateCluster = new ArrayList<Integer>();
		candidateCluster.addAll(edges.get(pivot));
		List<Integer> reducedCluster = new ArrayList<Integer>();
		while(candidateCluster.size()>1) {
			//Keep only edges with degree greater or equal than the size of the candidate cluster
			reducedCluster.clear();
			int n = candidateCluster.size();
			//Sort nodes by degree
			List<NodeIndexPlusDegree> nodesToProcess = new ArrayList<NodeIndexPlusDegree>();
			for(int i:candidateCluster) nodesToProcess.add(new NodeIndexPlusDegree(i, edges.get(i).size()));
			Collections.sort(nodesToProcess);
			for(NodeIndexPlusDegree node:nodesToProcess) {
				if(node.getDegree()<n)n--;
				else reducedCluster.add(node.getIndex());
			}
			
			if(reducedCluster.size()<candidateCluster.size()) {
				candidateCluster.clear();
				candidateCluster.addAll(reducedCluster);
			}
			int inconsistentNode = validateClique(candidateCluster,consistencyMatrix); 
			if(inconsistentNode==-1) break;
			candidateCluster.clear();
			for(int i:reducedCluster) if(i!=inconsistentNode) candidateCluster.add(i);
		}
		Collections.sort(candidateCluster);
		return candidateCluster;
	}

	private static int validateClique(List<Integer> cluster, boolean[][] consistencyMatrix) {
		for(int i:cluster) {
			for(int j:cluster) {
				if(!consistencyMatrix[i][j]) {
					return i;
				}
			}
		}
		return -1;
	}
}
class NodeIndexPlusDegree implements Comparable<NodeIndexPlusDegree> {
	private int index;
	private int degree;
	
	/**
	 * @param index
	 * @param degree
	 */
	public NodeIndexPlusDegree(int index, int degree) {
		super();
		this.index = index;
		this.degree = degree;
	}

	@Override
	public int compareTo(NodeIndexPlusDegree o2) {
		if(this.degree!=o2.degree) return this.degree-o2.degree;
		return this.index-o2.index;
	}

	public int getIndex() {
		return index;
	}

	public int getDegree() {
		return degree;
	}
}
