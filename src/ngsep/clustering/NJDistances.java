package ngsep.clustering;

/**
 * Utility class containing methods to calculate distances
 * between nodes during a neighbor joining clustering procedure
 */
public final class NJDistances {

    // This class can't be instantiated
    private NJDistances () {

    }

    /**
     * Given two nodes u and v, calculates the distance of each
     * one to a new node x. Returns a pair (dux, dvx) with the distance between u
     * and x, and the distance between v and x respectively
     * @param D - Distance matrix
     * @param rowSumVector - A vector defined as rowSumVector_i = 1 / (n - 2) * \sum_{j=1}^n D_{i j} for all 1<= i <= n
     * @param neighbors - The pair of neighbors (u, v) to be joined
     * @return double[] array with two entries with the distance between u and x, and the distance between v and x respectively
     */
    public static double [] distanceBetweenNeighbors (double[][] D, double[] rowSumVector, int u, int v) {
        double dux = 0.5 * (D[u][v] + rowSumVector[u] - rowSumVector[v]);
        double dvx = D[u][v] - dux;
        double [] answer = {dux,dvx};
        return answer;
    }

    /**
     * Given a new node x that joins the pair of existing nodes (u, v), calculates the distance
     * from x to any old node (present in the distance matrix).
     * @param D - Distance matrix
     * @param newNode - Node x characterized by the pair of neighbors (u, v)
     * @param oldNode - Index of an old node present in the distance matrix
     * @return the distance from x to any old node (present in the distance matrix)
     */
    public static double distanceBetweenNewAndOldNode (double[][] D,int u, int v, int oldNode) {
        return 0.5 * (D[u][oldNode] + D[v][oldNode] - D[u][v]);
    }
}
