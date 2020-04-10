package ngsep.clustering;

import ngsep.main.ProgressNotifier;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Sebastian Lemus Cadena
 */
public class UPGMA implements DistanceMatrixClustering{

    private Logger log = Logger.getLogger(UPGMA.class.getName());
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
     * Structure to memorize the accumulated distance for a given cluster
     * (sub-tree)
     */
    private ArrayList<Double> accumulatedClusterDistances;

    /**
     * Initial set of trees and accumulated distances
     * @param names - Names of the nodes
     */
    public void initializeSubTrees(List<String> names){
        int n = names.size();
        subTrees = new ArrayList<>(n);

        for (int i = 0; i < n; i++) {
            subTrees.add(new Dendrogram(names.get(i)));
        }

        accumulatedClusterDistances = new ArrayList<Double>(Collections.nCopies(n, 0.0));
    }

    /**
     * Executes UPGMA algorithm .
     * @param distances Matrix of distances
     * @return the resulting Dendrogram from the distance matrix
     */
    @Override
    public Dendrogram buildDendrogram(DistanceMatrix distances) {

        initializeSubTrees(distances.getIds());

        int n = distances.getNumSamples();
        DistanceMatrix oldMatrix = distances;

        for (int i = 0; i < n - 1; i++) {
            DistanceMatrix newMatrix = recalculateMatrix(oldMatrix);
            oldMatrix = newMatrix;
        }

        return subTrees.get(0);
    }

    /**
     * Picks the cell with the least distance in the matrix, makes a new tree for the
     * corresponding two nodes and recalculates the distances from this tree
     * to the other nodes
     * @param matrix - Distance Matrix
     * @return - New matrix with the updated distances
     */
    private DistanceMatrix recalculateMatrix(DistanceMatrix matrix){

        // Extract distances and names
        double[][] D = matrix.getDistances();
        ArrayList<String> names = new ArrayList<>();
        names.addAll(matrix.getIds());
        int n = D.length;

        // find minimum distance and its coordinates
        int x = 0;
        int y = 0;
        double min = Integer.MAX_VALUE;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j && D[i][j] < min){
                    min = D[i][j];
                    x = i;
                    y = j;
                }
            }
        }

        // Make a copy of the matrix without the x and y coordinates
        double[][] newD = copyMatrix(D, x, y);

        // Calculate distace for the new node between the clusters with coordinates x and y
        int A = names.get(x).split("!").length;
        int B = names.get(y).split("!").length;
        calculateDistancesFromNewNode(D, newD, x, y, A, B);

        // update the names according to the new matrix
        ArrayList<String> newNames = updateNames(names, x, y);
        String newNode = newNames.get(newNames.size() - 1);

        // include new sub-tree with the two joined clusters
        updateSubTrees(x, y, D[x][y]/2, newNode);


        return new DistanceMatrix(newNames, newD);
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
     * Calculates the distances fron the newly formed node (between clusters A and B)
     * to the rest of the nodes in the matrix.
     * @param D - Distance matrix
     * @param A - Cardinality of one of the joined clusters
     * @param B - Cardinality of one of the joined clusters
     */
    private void calculateDistancesFromNewNode(double[][] D, double newD[][], int x, int y, int A, int B){

        int n = D.length;
        int m = n - 1;
        int c = 0;

        for (int i = 0; i < n; i++) {
            if (c == m - 1){
                newD[c][m - 1] = 0;
            }else {
                if (i != x && i != y){
                    double newDistance = (D[i][x]*A + D[i][y]*B) / (A + B);
                    newD[c][m - 1] = newDistance;
                    newD[m - 1][c] = newDistance;
                    c ++;
                }
            }
        }
    }

    /**
     * Creates a new subtree pairing the nodes x and y with a new node u, and the respective distances
     * dx and dy
     * @param x - Node to be paired
     * @param y - Node to be paired
     * @param delta - Corresponds to D[x][y] / 2 by the ultrametricity assumption
     * @param newNode - The name of the new node u
     */
    private void updateSubTrees(int x, int y, double delta, String newNode){

        // create new tree with the nodes
        ArrayList<Dendrogram> newSubTrees = updateArray(Dendrogram.class, subTrees, x, y);
        Dendrogram left = subTrees.get(x);
        double leftAccDistance = accumulatedClusterDistances.get(x);
        Dendrogram right = subTrees.get(y);
        double rightAccDistance = accumulatedClusterDistances.get(y);

        DendrogramEdge arcLeft = new DendrogramEdge(delta - leftAccDistance, left);
        DendrogramEdge arcRight = new DendrogramEdge(delta - rightAccDistance, right);
        Dendrogram newTree = new Dendrogram(newNode);

        // update accumulated cluster distances
        accumulatedClusterDistances = updateArray(Double.class, accumulatedClusterDistances, x, y);
        accumulatedClusterDistances.add(delta);

        ArrayList<DendrogramEdge> children = new ArrayList<>();
        children.add(arcLeft);
        children.add(arcRight);
        newTree.setChildren(children);

        newSubTrees.add(newTree);

        subTrees = newSubTrees;
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
}
