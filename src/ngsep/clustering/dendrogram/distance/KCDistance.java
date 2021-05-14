package ngsep.clustering.dendrogram.distance;

import ngsep.clustering.Pair;
import ngsep.clustering.dendrogram.Dendrogram;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;

public class KCDistance implements BiFunction<Dendrogram, Dendrogram, Double> {

    private final double lambda;

    public KCDistance (double lambda) throws Exception {
        if (lambda > 1.0 || lambda < 0.0) {
            throw new Exception ("Lambda parameter has to be between 0 and 1");
        }
        this.lambda = lambda;
    }

    private List<Pair<Integer, Integer>> combine2 (List<Integer> leaves) {
        int k = leaves.size();
        ArrayList<Pair<Integer, Integer>> cmb = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            for (int j = i + 1; j < k; j++) {
                cmb.add(new Pair<>(
                        leaves.get(i), leaves.get(j)
                ));
            }
        }
        return cmb;
    }

    private Vector getTopologyVector (LCA lcaFinder, List<Integer> leaves) {
        int k = leaves.size();
        int n = k * (k + 1) / 2;
        List<Pair<Integer, Integer>> cmb = combine2(leaves);
        return new Vector(n, i -> {
            if (i < cmb.size()) {
                int u = cmb.get(i).first;
                int v = cmb.get(i).second;
                return lcaFinder.hopsFromRootToLCA(u, v) + 0.0;
            } else return  1.0;
        });
    }

    private Vector getDistanceVector (LCA lcaFinder, List<Integer> leaves) {
        int k = leaves.size();
        int n = k * (k + 1) / 2;
        List<Pair<Integer, Integer>> cmb = combine2(leaves);
        return new Vector(n, i -> {
            if (i < cmb.size()) {
                int u = cmb.get(i).first;
                int v = cmb.get(i).second;
                return lcaFinder.distanceFromRootToLCA(u, v);
            } else {
                int u = leaves.get(i - cmb.size());
                return lcaFinder.distanceToKthAncestor(u, 1);
            }
        });
    }

    private Vector getTreeVector (
            List<List<Pair<Integer, Double>>> adj,
            List<Integer> leaves
    ) throws Exception {
        LCA lcaFinder = new LCA(adj);
        Vector m = getTopologyVector(lcaFinder, leaves);
        Vector M = getDistanceVector(lcaFinder, leaves);
        return m.multiply(1.0 - lambda).add(M.multiply(lambda));
    }

    @Override
    public Double apply(Dendrogram t1, Dendrogram t2) {
        List<List<Pair<Integer, Double>>> adj1 = t1.toAdjacencyList();
        List<List<Pair<Integer, Double>>> adj2 = t2.toAdjacencyList();
        List<Integer> leaves1 = t1.getLeaves().stream()
                .sorted(Comparator.comparing(Dendrogram::getLabel))
                .map(Dendrogram::getId)
                .collect(Collectors.toList());
        List<Integer> leaves2 = t2.getLeaves().stream()
                .sorted(Comparator.comparing(Dendrogram::getLabel))
                .map(Dendrogram::getId)
                .collect(Collectors.toList());
        try {
            Vector v1 = getTreeVector(adj1, leaves1);
            Vector v2 = getTreeVector(adj2, leaves2);
            return Vector.distance(v1, v2);
        } catch (Exception e) {
            e.printStackTrace();
            return -1.0;
        }
    }

    @Override
    public <V> BiFunction<Dendrogram, Dendrogram, V> andThen(Function<? super Double, ? extends V> after) {
        return (t1, t2) -> after.apply(this.apply(t1, t2));
    }

    public static void main(String[] args) throws Exception{
        KCDistance kc = new KCDistance(0.5);
        Dendrogram t1 = Dendrogram.fromNewick("((A: 1.2,B: 0.8):0.5,(C:0.8,D:1.0):1.1);");
        Dendrogram t2 = Dendrogram.fromNewick("(((A:0.8,B:1.4):0.3,C:0.7):0.9,D:1.0);");
        t1.printTree(System.out);
        t2.printTree(System.out);
        System.out.println(kc.apply(t1, t2));
    }
}
