package ngsep.clustering;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;


public class NeighborJoining implements DistanceMatrixClustering {

	private Logger log = Logger.getLogger(NeighborJoining.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	
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

	/**
	 * Structure to memorize the subtrees that are created by the algorithm
	 */
	private ArrayList<Dendrogram> subTrees;

	/**
	 * Distance Matrix
	 */
	private DistanceMatrix distanceMatrix;

	/**
	 * Constructor
	 * @param distanceMatrix
	 */
	public NeighborJoining(DistanceMatrix distanceMatrix){
		this.distanceMatrix = distanceMatrix;
		initializeSubTrees(distanceMatrix.getIds());
	}

	public NeighborJoining(){

	}

	/**
	 * Initial set of trees
	 * @param names - Names of the nodes
	 */
	public void initializeSubTrees(List<String> names){
		subTrees = new ArrayList<>(names.size());

		for (int i = 0; i < subTrees.size(); i++) {
			subTrees.add(new Dendrogram(names.get(i)));
		}

	}
	 
	public static void main (String [ ] args) throws Exception {
			
	 	NeighborJoining nj = new NeighborJoining();
		//Parameters
		int k=CommandsDescriptor.getInstance().loadOptions(nj, args);
		
		String matrixFile = args[k++];
	 	DistanceMatrix dm = new DistanceMatrix(matrixFile);
	 	nj.initializeSubTrees(dm.getIds());
		Dendrogram njTree = nj.buildDendrogram(dm);
		njTree.printTree(System.out);	
	}
	
	
	public Dendrogram buildDendrogram(DistanceMatrix distances) {
	
		double iterableMatrix[][] = distances.getDistances();
		
		int nSamples = distances.getNumSamples();
		
		int nodesToAssign = nSamples;
		
		Map<String, Dendrogram> nodesMerged = new HashMap<>();
		
		Dendrogram njTree = new Dendrogram("");
		
		ArrayList<String> nodesList = new ArrayList<>(distances.getIds());
		
		
		while(nodesToAssign>2){
			
			double njMatrix[][]= new double[nodesToAssign][nodesToAssign];
			double totalDistance[]= new double[nodesToAssign];
			
			//Calculate total distance array
			for(int i=0;i < nodesToAssign; i++){
				for(int j=0;j < nodesToAssign; j++){
					totalDistance[i] += iterableMatrix[i][j];
				}  
			}
			
			//Compute Neighbor Joining matrix
			for(int i=0;i < nodesToAssign; i++){
				for(int j=0;j < i; j++){
					njMatrix[i][j] = njMatrix[j][i] = ((nodesToAssign - 2) * iterableMatrix[i][j]) - totalDistance[i] - totalDistance[j];
				}  
			}
			
			int rowMin = -1;
			int colMin = -1;
			double minValue = Float.POSITIVE_INFINITY;
			
			//Minimum value in matrix
			//need to be improve searching only in the triangle
			for(int i=0;i < nodesToAssign; i++){
				for(int j=0;j < nodesToAssign; j++){
					if(njMatrix[i][j]<minValue && i != j){
						rowMin = i;
						colMin = j;
						minValue = njMatrix[i][j];
					}
				}  
			}
			
						
			nodesToAssign--;
			
			if(nodesToAssign > 2){
				//branch length estimation
				double leftTreeDistance = (0.5f * iterableMatrix[rowMin][colMin]) + ((totalDistance[rowMin] - totalDistance[colMin] ) / (2 * (nodesToAssign - 1)));
				double rightTreeDistance  = (0.5f * iterableMatrix[rowMin][colMin]) + ((totalDistance[colMin] - totalDistance[rowMin])/ (2 * (nodesToAssign - 1)));

				//Central node join and newick formatting		
				String nodeMergeId = "("+nodesList.get(rowMin) +":"+leftTreeDistance+","+ nodesList.get(colMin)+":"+rightTreeDistance+")";
				
				//Try to get back existing trees
				Dendrogram leftTree = nodesMerged.get(nodesList.get(rowMin)); 
				Dendrogram rightTree = nodesMerged.get(nodesList.get(colMin)); 
				
				//If trees not exist create new ones
				if(leftTree == null){
					leftTree = new Dendrogram(nodesList.get(rowMin));
				}
				
				if(rightTree == null){
					rightTree = new Dendrogram(nodesList.get(colMin));
				}
				
				//Make a new tree with the two joined nodes
				List<DendrogramEdge> children = new ArrayList<>();
				children.add(new DendrogramEdge(1, leftTree));
				children.add(new DendrogramEdge(1, rightTree));
				njTree = new Dendrogram( nodeMergeId, children);			
				nodesMerged.put(nodeMergeId, njTree);
				
				//Next iteration nodes positions
				ArrayList<String> currentNodes = new ArrayList<>();
				currentNodes.add(nodeMergeId);
				for(int n=0;n<nodesList.size();n++){
					if(n != rowMin && n != colMin){ //erase joined nodes from nodes list
						currentNodes.add(nodesList.get(n));
					}

				}
				
				// distance matrix update
				double updateDistanceMatrix[][] = new double[nodesToAssign][nodesToAssign];
				
				//adding new row, col and bringing back old distances
				for(int i=1;i<nodesToAssign;i++){
					for(int j=0;j<i;j++){
						int index_i = nodesList.indexOf(currentNodes.get(i));
						if(j==0){ //if new node, calculate distance
							updateDistanceMatrix[j][i] = updateDistanceMatrix[i][j] = 0.5f * (iterableMatrix[rowMin][index_i]+iterableMatrix[colMin][index_i] - iterableMatrix[rowMin][colMin]);
						} else { //else bring old ones
							int index_j = nodesList.indexOf(currentNodes.get(j));
							updateDistanceMatrix[j][i] = updateDistanceMatrix[i][j] = iterableMatrix[index_i][index_j];
						}
						
					}
					
				}
				
			
				nodesList= currentNodes;
				
				iterableMatrix = updateDistanceMatrix;
				
			} else { // FINAL JOIN ---------------------------------------------------------
				
				
				double leftTreeDistance = (0.5f * iterableMatrix[1][2]) + (0.5f * (totalDistance[1] - totalDistance[2]));
				double rightTreeDistance  = (0.5f * iterableMatrix[1][2]) + (0.5f * (totalDistance[2] - totalDistance[1]));
				
				double centralTreeDistance = 0.0f;
				
				//Calculate final distance between join node and the last two
				if(rowMin!=0){
					if(rowMin<colMin || colMin == 0){
						centralTreeDistance = (0.5f * iterableMatrix[0][rowMin]) + (0.5f * (totalDistance[rowMin] - totalDistance[0]));
					} else{
						centralTreeDistance = (0.5f * iterableMatrix[0][colMin]) + (0.5f * (totalDistance[colMin] - totalDistance[0]));
					}
				} else{
					centralTreeDistance = (0.5f * iterableMatrix[0][colMin]) + (0.5f * (totalDistance[colMin] - totalDistance[0]));
				}
				
							
				String nodeMergeId = nodesList.get(0);
				
				// Newick format
				Dendrogram leftTree = new Dendrogram("("+nodesList.get(2)+":"+rightTreeDistance+","+nodesList.get(1)+":"+leftTreeDistance+",");
				Dendrogram rightTree = new Dendrogram(");\n");
			
				njTree = new Dendrogram(nodeMergeId+":"+centralTreeDistance);			

			}
			
		}

		
		return njTree;
	}

	public Dendrogram execute(DistanceMatrix matrix){
		int n = matrix.getNumSamples();
		DistanceMatrix oldMatrix = matrix;

		for (int i = n - 3; i >= 0; i++) {
			DistanceMatrix newMatrix = recalculateMatrix(oldMatrix);
			oldMatrix = newMatrix;
		}

		Dendrogram dendrogram = subTrees.get(0);
		return dendrogram;
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

		int x = 0;
		int y = 0;
		double min = Double.MAX_VALUE;

		// Find minimum value in Q
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j && Q[i][j] < min){
					min = Q[i][j];
					x = i;
					y = j;
				}
			}
		}

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
			double dz = D[n - 1][x] - dx;
			ArrayList<String> finalNodes = updateNames(names, 0, 1);
			String finalNode = finalNodes.get(finalNodes.size() - 1);
			updateSubTrees(0, 1, dz, 0, finalNode);
		}

		return new DistanceMatrix(newNames, newD);

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

		ArrayList<DendrogramEdge> children = newTree.getChildren();
		children.add(arcLeft);
		children.add(arcRight);

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
		int c = 0;
		ArrayList<T> newArr = new ArrayList<>(n - 1);

		for (int i = 0; i < n; i++) {
			if (i != x && i != y){
				newArr.add(arr.get(i));
				c++;
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

