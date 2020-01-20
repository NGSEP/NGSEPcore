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
	
	/**
	 * Structure to memorize the subtrees that are created by the algorithm
	 */
	private ArrayList<Dendrogram> subTrees;
	
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

	/**
	 * Constructor
	 * @param distanceMatrix
	 */
	public NeighborJoining(DistanceMatrix distanceMatrix){
		initializeSubTrees(distanceMatrix.getIds());
	}

	public NeighborJoining(){

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
		initializeSubTrees(dm.getIds());
		Dendrogram njTree = buildDendrogram(dm);
		if(outputFile == null) System.out.println(njTree.toNewick());
		else {
			try (PrintStream out = new PrintStream(outputFile)) {
				out.println(njTree.toNewick());
			}
		}
		
 	}
 	
 	/**
	 * Initial set of trees
	 * @param names - Names of the nodes
	 */
	public void initializeSubTrees(List<String> names){
		subTrees = new ArrayList<>(names.size());

		for (int i = 0; i < names.size(); i++) {
			subTrees.add(new Dendrogram(names.get(i)));
		}

	}

	/**
	 * Updates the distance matrix, pairs the closest nodes in the matrix and creates a new tree from these nodes.
	 * Q(i,j) = (n - 2)D(i,j) - \sum_{k = 1}^n D(i,k) - \sum_{k = 1}^n D(i,k)
	 * @param matrix - New distance matrix reduced in one dimension
	 * @return - Updated distance matrix with one less dimension and the new created node
	 */
	private DistanceMatrix recalculateMatrix(DistanceMatrix matrix){

		// Extract distances and names
		double[][] D = matrix.getDistances();
		ArrayList<String> names = new ArrayList<>();
		names.addAll(matrix.getIds());
		int n = D.length;

		// Calculate auxiliary Q matrix
		double[][] Q = calculateQMatrix(D);

		// Find minimum value in Q
		Tuple<Integer, Integer> minCoordinates = findMinimumInMatriz(Q);
		int x = minCoordinates.first;
		int y = minCoordinates.second;

		// make a copy of the matrix with the values not in the x or y rows
		double[][] newD = copyMatrix(D, x, y);

		// Distances from the pair of nodes to a new node u
		Tuple<Double, Double> distances = calculateDistanceFromAPairToNewNode(D, x, y);
		double dx = distances.first;
		double dy = distances.second;
		int j = 0;

		// Distance from all nodes to the new node u
		for (int i = 0; i <= n - 1; i++) {
			if (i != x && i != y){
				if (j == n - 2){
					newD[i][i] = 0;
				}else {
					double newDist = calculateDistanceFromTaxaToNewNode(D, x, y, i);
					newD[n - 2][j] = newDist;
					newD[j][n - 2] = newDist;
				}
				j++;
			}
		}

		// Update names and subtrees in accordance to the new node created
		// The new node remains at the end of the matrix in each call of this method
		ArrayList<String> newNames = updateNames(names, x, y);
		String newNode = newNames.get(newNames.size() - 1);
		updateSubTrees(x, y, dx, dy, newNode);

		if (n == 3){ // Last update

			// Find un-joined node z
			int z = 0;
			int[] lastNodes = new int[3];
			lastNodes[x] = 1;
			lastNodes[y] = 1;
			for (int i = 0; i < 3; i++) {
				if (lastNodes[i] != 1) z = i;
			}

			double dz = calculateDistanceFromTaxaToNewNode(D, x, y, z);
			ArrayList<String> finalNodes = updateNames(newNames, 0, 1);
			String finalNode = finalNodes.get(finalNodes.size() - 1);
			updateSubTrees(0, 1, dz/2, dz/2, finalNode);
		}

		return new DistanceMatrix(newNames, newD);

	}

	public Dendrogram buildDendrogram(DistanceMatrix matrix){
		int n = matrix.getNumSamples();
		DistanceMatrix oldMatrix = matrix;

		for (int i = n - 3; i >= 0; i--) {
			DistanceMatrix newMatrix = recalculateMatrix(oldMatrix);
			oldMatrix = newMatrix;
		}

		Dendrogram dendrogram = subTrees.get(0);
		return dendrogram;
	}

	private Tuple<Integer, Integer> findMinimumInMatriz(double[][] Q){
		int n = Q.length;
		int x = 0;
		int y = 0;
		double min = Double.MAX_VALUE;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j && Q[i][j] < min){
					min = Q[i][j];
					x = i;
					y = j;
				}
			}
		}
		return new Tuple<>(x, y);
	}

	/**
	 * Calculates intermediate distance matrix for neighbor joining algorithm
	 * Q(i,j) = (n - 2)*D(i,j) - \sum_{k = 0}^{n - 1} D(i,k) - \sum_{k = 0}^{n - 1} D(j,k)
	 * @param D - Distance Matrix
	 * @return - New distance matrix
	 */
	private double[][] calculateQMatrix(double[][] D){
		int n = D.length;
		double[][] Q = new double[n][n];

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				Q[i][j] = (n - 2)*D[i][j] - sumOverRow(D,i) - sumOverRow(D,j);
			}
		}
		return Q;
	}

	/**
	 * Calculates the sum over the ith row of the matrix D
	 * @param D - Distance matrix
	 * @param i - row index
	 * @return The resulting sum
	 */
	private double sumOverRow(double[][] D, int i){
		double sum = 0;
		for (int j = 0; j < D[i].length; j++) {
			sum += D[i][j];
		}
		return sum;
	}

	/**
	 * Makes a copy of the matrix without the x and y rows/colums
	 * @param D - Distance matrix
	 * @param x - Row/column to remove
	 * @param y - Row/column to remove
	 * @return A new distance matrix with one less dimension
	 */
	private double[][] copyMatrix(double[][] D, int x, int y){

		int n = D.length;
		double[][] newD = new double[n - 1][n - 1];
		int d = newD.length - 1;
		int c = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != x && j != x && i != y && j != y){
					newD[c / d][c % d] = D[j][i];
					c++;
				}
			}
		}
		return newD;
	}

	/**
	 * Calculates the distance from a pair of nodes f and g to a new node u
	 * @param D - Distance matrix
	 * @param f - Node to be paired
	 * @param g - Node to be paired
	 * @return Distances from both nodes to the new node u
	 */
	private Tuple<Double, Double> calculateDistanceFromAPairToNewNode(double[][] D, int f, int g){
		double d_f_u = 0.5f*D[f][g] + (1f/(2.0f*(D.length - 2)))*(sumOverRow(D,f) - sumOverRow(D,g));
		double d_g_u = D[f][g] - d_f_u;
		return new Tuple<>(d_f_u, d_g_u);
	}

	/**
	 * Calculates the distance from a new given node u (which pairs the nodes f and g)
	 * and the rest of the nodes in the tree
	 * @param D - Distance matrix
	 * @param f - Paired node by u
	 * @param g - Paired node by u
	 * @param k - kth node in the tree (different from f and g)
	 * @return
	 */
	private double calculateDistanceFromTaxaToNewNode(double[][] D, int f, int g,  int k){
		return 0.5f*(D[f][k] + D[g][k] - D[f][g]);
	}

	/**
	 * Updates the array of names without the x and y positions
	 * @param names - Array of names
	 * @param x - Element to omit
	 * @param y - Element to omit
	 * @return Updated array of names
	 */
	private ArrayList<String> updateNames(ArrayList<String> names, int x, int y){
		// update names in accordance to the new matrix
		ArrayList<String> newNames = updateArray(String.class, names, x, y);
		String newNode = names.get(x) + "!" + names.get(y);
		newNames.add(newNode);

		return newNames;
	}

	/**
	 * Creates a new subtree pairing the nodes x and y with a new node u, and the respective distances
	 * dx and dy
	 * @param x - Node to be paired
	 * @param y - Node to be paired
	 * @param dx - Distance from x to u
	 * @param dy - Distance from y to u
	 * @param newNode - The name of the new node u
	 */
	private void updateSubTrees(int x, int y, double dx, double dy, String newNode){

		// create new tree with the nodes
		ArrayList<Dendrogram> newSubTrees = updateArray(Dendrogram.class, subTrees, x, y);
		Dendrogram left = subTrees.get(x);
		Dendrogram right = subTrees.get(y);

		DendrogramEdge arcLeft = new DendrogramEdge(dx, left);
		DendrogramEdge arcRight = new DendrogramEdge(dy, right);
		Dendrogram newTree = new Dendrogram(newNode);

		ArrayList<DendrogramEdge> children = new ArrayList<>();
		children.add(arcLeft);
		children.add(arcRight);
		newTree.setChildren(children);

		newSubTrees.add(newTree);

		subTrees = newSubTrees;
	}


	/**
	 * Updates an array except for the two coordinates x and y
	 * @param type - Type of the array
	 * @param arr - Array of elements
	 * @param x - Element to omit
	 * @param y - Element to omit
	 * @return new copy of the given array
	 */
	private <T> ArrayList<T> updateArray(Class<T> type, ArrayList<T> arr, int x, int y){
		int n = arr.size();
		ArrayList<T> newArr = new ArrayList<>(n - 1);

		for (int i = 0; i < n; i++) {
			if (i != x && i != y){
				newArr.add(arr.get(i));
			}
		}
		return newArr;
	}


	private class Tuple<A, B> {
		public final A first;
		public final B second;

		public Tuple(A first, B second) {
			this.first = first;
			this.second = second;
		}
	}

}

