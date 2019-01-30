package ngsep.assembly;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

@FunctionalInterface
public interface PathsFinder {
	public List<List<Integer>> findPaths(Map<Integer, List<Edge>> edges);

	public static PathsFinder NONE = (Map<Integer, List<Edge>> edges) -> {
		return new LinkedList<List<Integer>>();
	};

}
