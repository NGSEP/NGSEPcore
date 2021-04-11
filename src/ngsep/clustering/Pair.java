package ngsep.clustering;

/**
 * A container representing a 2-tuple
 * @param <A> - Type of the first element
 * @param <B> - Type of the second element
 */
public class Pair <A, B> {
    public final A first;
    public final B second;
    public Pair (A first, B second) {
        this.first = first;
        this.second = second;
    }

    @Override
    public String toString() {
        return "(" + this.first + ", " + this.second + ")";
    }
}
