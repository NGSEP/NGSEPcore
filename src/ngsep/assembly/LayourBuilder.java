package ngsep.assembly;

import java.util.LinkedList;
import java.util.List;

@FunctionalInterface
public interface LayourBuilder {
	public List<List<Integer>> findPaths(AssemblyGraph graph);

	public static LayourBuilder NONE = (AssemblyGraph graph) -> {
		return new LinkedList<List<Integer>>();
	};
}
