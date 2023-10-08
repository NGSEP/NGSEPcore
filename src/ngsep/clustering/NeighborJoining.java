package ngsep.clustering;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;


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
	 * @param neighbors - Pair of indices to remove
	 * @param <A> - Parameterized type of the list
	 * @return A new list without the elements specified by the pair of indices
	 */
	private <A> List<A> deleteFromList (List<A> xs, Pair<Integer, Integer> neighbors) {
		int n = xs.size();
		int u = neighbors.first;
		int v = neighbors.second;
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
	private double[] rowSums (double[][] D) {
		int n = D.length;
		double[] a = new double[n];
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
	private Pair<Integer, Integer> findNeighbors (
			double[][] D,
			double[] rowSumVector
	) {
		int n = D.length;
		int u = 0;
		int v = 0;
		double minS = Double.MAX_VALUE;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				double S = D[i][j] - rowSumVector[i] - rowSumVector[j];
				if (S < minS) {
					minS = S;
					u = i;
					v = j;
				}
			}
		}
		return new Pair<>(u, v);
	}

	/**
	 *
	 * @param D - Distance matrix of size n
	 * @param neighbors - The pair of neighbors (u, v) joined by x
	 * @return A new distance matrix D' of size n - 1, in which the nodes/indices
	 * (u, v) are replaced are merged into a new node x, and is represented as
	 * the last row (and column) in D'. Distances from x to the other nodes is calculated
	 * using distanceBetweenNewAndOldNode.
	 */
	private double[][] recalculateDistances (
			double[][] D,
			Pair<Integer, Integer> neighbors
	) {
		int n = D.length;
		int u = neighbors.first;
		int v = neighbors.second;
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
				newD[x][n - 2] = NJDistances
						.distanceBetweenNewAndOldNode(D, u, v, i);
				newD[n - 2][x] = newD[x][n - 2];
				x++;
			}
		}

		return newD;
	}

	// Neighbor joining algorithm

	/**
	 * Receives a list of nodes (trees) and the distance matrix that specifies the pairwise distance
	 * between each node.
	 * Runs an iteration of the neighbor joining algorithm, described by the following steps:
	 * 1. Calculates a vector rowSumVector_i = 1 / (n - 2) * \sum_{j=1}^n D_{i j} for all 1<= i <= n to use
	 * when finding neighbors and calculating their distances.
	 * 2. Find a pair of neighbors (u, v) in accordance to the findNeighbors function.
	 * 3. Join nodes (u, v) to a new node x calculating their distances (dux, dvx).
	 * 4. Recalculate the labels of the trees and distance matrix taking into account the new merged node x.
	 * Return the new distance matrix and the new list of trees
	 * @param matrix - Distance matrix
	 * @param subtrees - List of node (trees) characterized by the distance matrix
	 * @return A pair with the new distance matrix and the new list of trees
	 */
	private Pair<DistanceMatrix, List<Dendrogram>> runNeighborJoining (
			DistanceMatrix matrix,
			List<Dendrogram> subtrees
	) {
		double[][] D = matrix.getDistances();
		int n = D.length;
		List<String> names = matrix.getIds();
		double[] rowSumVector = n == 2 ? new double[]{0, 0} : rowSums(D);
		Pair<Integer, Integer> neighbors = findNeighbors(D, rowSumVector);
		int u = neighbors.first;
		int v = neighbors.second;
		String newNodeName = names.get(u) + "!" + names.get(v);
		double distance1 = 0.5 * D[u][v];
		double distance2 = 0.5 * D[u][v];
		if(n>2) {
			double [] distances = NJDistances.distanceBetweenNeighbors(D, rowSumVector, u,v);
			distance1 = distances[0];
			distance2 = distances[1];
		}
				
		Dendrogram newNode = Dendrogram.join2(newNodeName, subtrees.get(u), distance1, subtrees.get(v), distance2);

		List<String> newNames = deleteFromList(names, neighbors);
		newNames.add(newNodeName);
		List<Dendrogram> newSubtrees = deleteFromList(subtrees, neighbors);
		newSubtrees.add(newNode);
		double[][] newD = recalculateDistances(D, neighbors);
		return new Pair<>(
				new DistanceMatrix(newNames, newD),
				newSubtrees
		);
	}

	/**
	 *
	 * @param names - The labels of each leaf in the tree to be constructed
	 * @return The initial set of trees, each one containing only one node
	 * labeled by the taxa to be clustered
	 */
	private List<Dendrogram> initializeSubtrees (List<String> names) {
		List<Dendrogram> subtrees = new ArrayList<>(names.size());
		for (String name : names) subtrees.add(new Dendrogram(name));
		return subtrees;
	}

	/**
	 * Clusters a given set of sequences characterized by a pairwise distance
	 * matrix. Runs the runNeighborJoining function until the resulting list
	 * of trees is reduced to only one tree.
	 * @param distances - Initial distance matrix
	 * @return - A binary tree (dendrogram) that clusters the given sequences.
	 */
	@Override
	public Dendrogram buildDendrogram(DistanceMatrix distances) {
		List<Dendrogram> subtrees = initializeSubtrees(distances.getIds());
		DistanceMatrix matrix = distances;
		while (subtrees.size() > 1) {
			Pair<DistanceMatrix, List<Dendrogram>> njResult = runNeighborJoining(matrix, subtrees);
			matrix = njResult.first;
			subtrees = njResult.second;
		}
		return subtrees.get(0);
	}
}

