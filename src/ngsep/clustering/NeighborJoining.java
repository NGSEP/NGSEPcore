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
package ngsep.clustering;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;

/**
 * @author Sebastian Lemus
 * @author Jorge Duitama
 * @author Cristian Loaiza
 *
 */
public class NeighborJoining implements DistanceMatrixClustering {

	// Constants for default values

	// Logging and progress
	private Logger log = Logger.getLogger(NeighborJoining.class.getName());
	private ProgressNotifier progressNotifier=null;

	//Parameters
	private String inputFile = null;
	private String outputFile = null;

	// Get and set methods
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}

	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}

	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}

	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}

 	public static void main (String [ ] args) throws Exception {
		NeighborJoining instance = new NeighborJoining();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}

 	public void run () throws IOException {
		DistanceMatrix dm;
		if(inputFile!=null) {
			log.info("Loading matrix from file "+inputFile);
			dm = new DistanceMatrix(inputFile);
		} else {
			log.info("Loading matrix from standard input");
			dm = new DistanceMatrix(System.in);
		}
		Dendrogram njTree = buildDendrogram(dm);
		if(outputFile == null) System.out.println(njTree.toNewick());
		else {
			try (PrintStream out = new PrintStream(outputFile)) {
				out.println(njTree.toNewick());
			}
		}

 	}

	// Empty constructor
	public NeighborJoining () {

	}

	// List helpers

	/**
	 * Returns a new List without the elements specified by the pair
	 * of indices in the given list
	 * @param xs - List of elements
	 * @param u - First index to remove
	 * @param v - Second index to remove
	 * @param <A> - Parameterized type of the list
	 * @return A new list without the elements specified by the pair of indices
	 */
	private <A> List<A> deleteFromList (List<A> xs, int u, int v) {
		int n = xs.size();
		ArrayList<A> ys = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			if (i != u && i != v) ys.add(xs.get(i));
		}
		return ys;
	}

	// Distance calculations

	/**
	 * @param D - Distance Matrix
	 * @return A vector from a distance matrix D of size n x n, defined as:
	 * rowSumVector_i = 1 / (n - 2) * \sum_{j=1}^n D_{i j} for all 1<= i <= n
	 */
	public static double[] neighborJoiningRowSums (double[][] D) {
		int n = D.length; 
		double[] a = new double[n];
		if(n==2) {
			Arrays.fill(a,0);
			return a;
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				a[i] += D[i][j];
			}
			a[i] /= (n - 2.0);
		}
		return a;
	}

	// Neighbor operations

	/**
	 * @param D - Distance matrix
	 * @param rowSumVector - A vector defined as rowSumVector_i = 1 / (n - 2) * \sum_{j=1}^n D_{i j} for all 1<= i <= n
	 * @return a pair of nodes (u, v) = argmin_{(i, j)} D_{i j} - rowSumVector_i - rowSumVector_j
	 */
	private int [] findNeighbors (double[][] D,double[] rowSumVector) {
		int n = D.length;
		int [] answer = new int[2];
		double minS = Double.MAX_VALUE;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				double S = D[i][j] - rowSumVector[i] - rowSumVector[j];
				if (S < minS) {
					minS = S;
					answer[0] = i;
					answer[1] = j;
				}
			}
		}
		
		return answer;
	}

	/**
	 *
	 * @param D - Distance matrix of size n
	 * @param u - first index to join
	 * @param v second index to join
	 * @return A new distance matrix D' of size n - 1, in which the nodes/indices
	 * (u, v) are replaced are merged into a new node x, and is represented as
	 * the last row (and column) in D'. Distances from x to the other nodes is calculated
	 * using distanceBetweenNewAndOldNode.
	 */
	private double[][] recalculateDistances (double[][] D, int u, int v) {
		int n = D.length;
		double[][] newD = new double[n - 1][n - 1];

		int x = 0;
		int y = 0;
		for (int i = 0; i < n; i++) {
			if (i != u && i != v) {
				for (int j = 0; j < n; j++) {
					if (j != u && j != v) {
						newD[x][y] = D[i][j];
						y++;
					}
				}
				y = 0;
				x++;
			}
		}

		x = 0;
		for (int i = 0; i < n; i++) {
			if (i != u && i != v){
				newD[x][n - 2] = NJDistances.distanceBetweenNewAndOldNode(D, u, v, i);
				newD[n - 2][x] = newD[x][n - 2];
				x++;
			}
		}

		return newD;
	}

	// Neighbor joining algorithm

	/**
	 * Clusters a given set of sequences characterized by a pairwise distance
	 * matrix. Runs the runNeighborJoining function until the resulting list
	 * of trees is reduced to only one tree.
	 * @param distances - Initial distance matrix
	 * @return - A binary tree (dendrogram) that clusters the given sequences.
	 */
	@Override
	public Dendrogram buildDendrogram(DistanceMatrix distances) {
		List<String> names = distances.getIds();
		List<Dendrogram> subtrees = new ArrayList<>(names.size());
		for (String name : names) subtrees.add(new Dendrogram(name));
		DistanceMatrix matrix = distances;
		double[][] D = matrix.getDistances();
		while (subtrees.size() > 1) {
			double[] rowSumVector = neighborJoiningRowSums(D);
			int [] neighbors = findNeighbors(D, rowSumVector);
			int u = neighbors[0];
			int v = neighbors[1];
			Dendrogram d1 = subtrees.get(u);
			Dendrogram d2 = subtrees.get(u);
			String newNodeName = d1.getLabel() + "!" + d2.getLabel();
			double [] neighborDistances = NJDistances.distanceBetweenNeighbors(D, rowSumVector, u,v);
			DendrogramEdge e1 = new DendrogramEdge(neighborDistances[0], subtrees.get(u));
			DendrogramEdge e2 = new DendrogramEdge(neighborDistances[1], subtrees.get(v));
			Dendrogram newNode = new Dendrogram(newNodeName, List.of(e1,e2));
			subtrees = deleteFromList(subtrees, u, v);
			subtrees.add(newNode);
			D = recalculateDistances(D, u, v);
		}
		return subtrees.get(0);
	}
}

