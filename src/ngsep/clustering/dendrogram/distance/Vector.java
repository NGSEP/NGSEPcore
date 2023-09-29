package ngsep.clustering.dendrogram.distance;

import java.util.Arrays;
import java.util.function.BiFunction;
import java.util.function.Function;

public class Vector {

    private final double[] v;

    public Vector (double[] v) {
        this.v = v;
    }

    public Vector (int size, double fillValue) {
        v = new double[size];
        for (int i = 0; i < size; i++) {
            v[i] = fillValue;
        }
    }

    public Vector (int size, Function<Integer, Double> init) {
        v = new double[size];
        for (int i = 0; i < size; i++) {
            v[i] = init.apply(i);
        }
    }

    public int size () {
        return v.length;
    }

    public double get (int i) {
        return v[i];
    }

    public Vector map (Function<Double, Double> f) {
        double[] u = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            u[i] = f.apply(v[i]);
        }
        return new Vector(u);
    }

    public Vector zipWith (BiFunction<Double, Double, Double> f, Vector u) throws Exception{
        if (this.size() != u.size()) throw new Exception("Vector should be of the same size");
        double[] w = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            w[i] = f.apply(v[i], u.get(i));
        }
        return new Vector(w);
    }

    public <A> A fold (BiFunction<A, Double, A> f, A seed) {
        A res = seed;
        for (double vi : v) {
            res = f.apply(res, vi);
        }
        return res;
    }

    public Vector multiply (double x) {
        return this.map(vi -> vi * x);
    }

    public Vector negative () {
        return this.multiply(-1.0);
    }

    public Vector add (Vector u) throws Exception {
        return this.zipWith(Double::sum, u);
    }

    public Vector subtract(Vector u) throws Exception {
        return this.add(u.negative());
    }

    public double dot (Vector u) throws Exception {
        return this.zipWith((vi, ui) -> vi * ui, u)
                .fold(Double::sum, 0.0);
    }

    public double norm () throws Exception {
        return Math.sqrt(this.dot(this));
    }

    public static double distance (Vector v, Vector u) throws Exception {
        return v.subtract(u).norm();
    }

    @Override
    public String toString() {
        return Arrays.toString(v);
    }
}
