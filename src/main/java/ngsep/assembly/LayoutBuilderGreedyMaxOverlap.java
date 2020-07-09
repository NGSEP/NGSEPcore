package ngsep.assembly;

import java.util.Comparator;

public class LayoutBuilderGreedyMaxOverlap extends LayoutBuilderGreedy {
	@Override
	public void findPaths(AssemblyGraph graph) {
		Comparator<AssemblyEdge> comparator = new Comparator<AssemblyEdge>() {
			@Override
			public int compare(AssemblyEdge edge1, AssemblyEdge edge2) {
				return edge2.getOverlap()-edge1.getOverlap();
			}
		};
		setComparator(comparator);
		super.findPaths(graph);
	}
}
