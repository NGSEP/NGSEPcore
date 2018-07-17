package ngsep.assembly;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class MaxPath extends GraphSimplificator {
    /** List with the answer */
    private List<Integer> path = new LinkedList<>();
    /** Array for maxLentghPath (Dynamic algorithm) */
    private int[] maxLentghPath = new int[N];
    /** Array for the decisions in maxLentghPath */
    private int[] sig = new int[N];

    public MaxPath(int N, Map<Integer, Map<Integer, Integer>> Edges) {
	super(N, Edges);
    }

    public List<Integer> maxLentghPath() {
	int[] top = getTopologicalOrder();
	for (int i = top.length - 1; i >= 0; i--)
	    maxLentghPath(top[i]);

	int pos = getMaxPosition();
	while (pos != -1) {
	    path.add(pos);
	    pos = sig[pos];
	}
	return path;
    }

    private void maxLentghPath(int n) {
	int peso, max = 0, posMax = -1;
	for (Entry<Integer, Integer> x : Edges.get(n).entrySet()) {
	    peso = maxLentghPath[x.getKey()] + x.getValue();
	    if (max < peso) {
		posMax = x.getKey();
		max = peso;
	    }
	}
	maxLentghPath[n] = max;
	sig[n] = posMax;
    }

    private int getMaxPosition() {
	int temp = 0, pos = 0;
	for (int i = 0; i < maxLentghPath.length; i++) {
	    if (maxLentghPath[i] > temp) {
		pos = i;
		temp = maxLentghPath[i];
	    }
	}
	return pos;
    }

    public List<Integer> getPath() {
	return path;
    }
}
