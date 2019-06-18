package ngsep.assembly;

import java.util.AbstractMap.SimpleEntry;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Set;

public class LayoutBuilderImplementation implements LayourBuilder {
	private Map<AssemblyVertex, Map<AssemblyVertex, AssemblyEdge>> VertEdge = new HashMap<>();
	private Map<AssemblyVertex, AssemblyVertex> post = new HashMap<>();
	private Map<AssemblyVertex, Integer> maxOvp = new HashMap<>();
	private Set<AssemblyVertex> poseccing = new HashSet<AssemblyVertex>();

	@Override
	public void findPaths(AssemblyGraph graph) {
		// new LayoutGraphFilter(graph);
		
		for (AssemblyEdge assemblyEdge : graph.getEdges()) {
			AssemblyVertex v1 = assemblyEdge.getVertex1();
			AssemblyVertex v2 = assemblyEdge.getVertex2();
			VertEdge.computeIfAbsent(v1, (AssemblyVertex x) -> new HashMap<>()).put(v2, assemblyEdge);
			VertEdge.computeIfAbsent(v2, (AssemblyVertex x) -> new HashMap<>()).put(v1, assemblyEdge);
		}
	
		//Order vertices by number of edges
		PriorityQueue<AssemblyVertex> quee = new PriorityQueue<AssemblyVertex>(new Comparator<AssemblyVertex>() {
			public int compare(AssemblyVertex o1, AssemblyVertex o2) {
				return VertEdge.get(o1).size() - VertEdge.get(o2).size();
			}
		});
		quee.addAll(VertEdge.keySet());

		
		while (!quee.isEmpty()) {
			AssemblyVertex assemblyVertex = quee.poll();
			if (!post.containsKey(assemblyVertex)) {
				calculateMaxOverlap(assemblyVertex);
			}
		}

		int max = 0;
		AssemblyVertex i = null;
		for (Entry<AssemblyVertex, Integer> entry : maxOvp.entrySet()) {
			if (max < entry.getValue()) {
				max = entry.getValue();
				i = entry.getKey();
			}
		}

		List<AssemblyEdge> ans = new LinkedList<AssemblyEdge>();
		AssemblyVertex next = null;
		while (i != null) {
			next = post.get(i);
			if (next == null)
				break;
			ans.add(VertEdge.get(i).get(next));
			i = next;
		}
		graph.addPath(ans);
//		graph.paths();
	}

	private int calculateMaxOverlap(AssemblyVertex assemblyVertex) {
		int max = 0;
		AssemblyVertex next = null;
		poseccing.add(assemblyVertex);

		List<Entry<AssemblyVertex, Integer>> verts = new LinkedList<>();
		for (Entry<AssemblyVertex, AssemblyEdge> entry : VertEdge.get(assemblyVertex).entrySet()) {
			AssemblyVertex to = entry.getKey();
			AssemblyEdge assemblyEdge = entry.getValue();
			if (poseccing.contains(to))
				continue;

			if (to.getIndex() / 2 == assemblyVertex.getIndex() / 2) {
				max = assemblyEdge.getOverlap() + calculateMaxOverlap(to);
				post.put(assemblyVertex, to);
				maxOvp.put(assemblyVertex, max);
				poseccing.remove(assemblyVertex);
				return max;
			}

			if (post.containsKey(to) && (to.getIndex() / 2 == post.get(to).getIndex() / 2)) {
				verts.add(new SimpleEntry<>(to, maxOvp.get(to)));
			}

			if (max < assemblyEdge.getOverlap()) {
				max = assemblyEdge.getOverlap();
				next = to;
			}
		}

		if (next != null) {
			max += calculateMaxOverlap(next);
			verts.add(new SimpleEntry<>(next, max));
		}

		max = 0;
		for (Entry<AssemblyVertex, Integer> en : verts) {
			if (max < en.getValue()) {
				max = en.getValue();
				next = en.getKey();
			}
		}
		post.put(assemblyVertex, next);
		maxOvp.put(assemblyVertex, max);
		poseccing.remove(assemblyVertex);
		return max;
	}

}
