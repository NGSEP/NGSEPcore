package ngsep.clustering.dendrogram.distance;

import ngsep.clustering.Pair;
import ngsep.clustering.dendrogram.Dendrogram;

import java.util.List;

public class LCA {

    private final List<List<Pair<Integer, Double>>> adj;
    private final int[] depth;
    private final int[][] up;
    private final int maxDepth;

    public LCA (Dendrogram t) {
        t.generateIds();
        int n = t.getSize();
        this.maxDepth = log2(n);
        this.adj = t.toAdjacencyList();
        this.depth = new int[n];
        this.up = new int[n][maxDepth];
        saveParents(0);
    }

    public LCA (List<List<Pair<Integer, Double>>> adjacencyList) {
        this.adj = adjacencyList;
        int n = adj.size();
        this.maxDepth = log2(n);
        this.depth = new int[n];
        this.up = new int[n][maxDepth];
        saveParents(0);
    }

    private int log2 (int n) {
        int k = 0;
        while (n > (1 << k)) k++;
        return k;
    }

    private void saveParents (int s) {
        for (Pair<Integer, Double> e : adj.get(s)) {
            int v = e.first;
            depth[v] = 1 + depth[s];
            up[v][0] = s;
            for (int j = 1; j < maxDepth; j++) {
                up[v][j] = up[up[v][j - 1]][j - 1];
            }
            saveParents(v);
        }
    }

    public int lowestCommonAncestor (int u, int v) {
        int[] p = depth[u] < depth[v] ? new int[]{v, u} : new int[]{u, v};
        int x = p[0];
        int y = p[1];

        int k = depth[x] - depth[y];
        for (int j = 0; j < maxDepth; j++) {
            if ((k & (1 << j)) != 0) {
                x = up[x][j];
            }
        }

        if (x == y) return x;

        for (int j = maxDepth - 1; j >= 0; j--) {
            if (up[x][j] != up[y][j]) {
                x = up[x][j];
                y = up[y][j];
            }
        }

        return up[x][0];
    }

    public int hopsFromRootToLCA (int u, int v) {
        int w = lowestCommonAncestor(u, v);
        return depth[w];
    }

    private double findDistance (int s, int dest, double acc) {
        if (s == dest) return acc;
        if (adj.get(s).isEmpty()) return -1.0;
        else {
            double ans = -1.0;
            for (Pair<Integer, Double> e : adj.get(s)) {
                int v = e.first;
                double weight = e.second;
                double res = findDistance(v, dest, acc + weight);
                if (res != -1.0) ans = res;
            }
            return ans;
        }
    }

    public double distanceFromRootToLCA (int u, int v) {
        int w = lowestCommonAncestor(u, v);
        return findDistance(0, w, 0.0);
    }

    public int findKthAncestorOf(int v, int k) {
        int x = v;
        for (int j = 0; j < maxDepth; j++) {
            if ((k & (1 << j)) != 0) x = up[x][j];
        }
        return x;
    }

    public int hopsToKthAncestor(int v, int k) {
        int ancestor = findKthAncestorOf(v, k);
        return depth[v] - depth[ancestor];
    }

    public double distanceToKthAncestor(int v, int k) {
        int ancestor = findKthAncestorOf(v, k);
        return findDistance(ancestor, v, 0.0);
    }
}
