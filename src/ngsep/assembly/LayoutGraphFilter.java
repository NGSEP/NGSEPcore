package ngsep.assembly;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.PriorityQueue;

public class LayoutGraphFilter {
	private AssemblyGraph Graph;
	private HashSet<AssemblyVertex> procesing;
	private HashSet<AssemblyVertex> finished;
	private HashMap<AssemblyVertex, AssemblyVertex> prev;
	private HashMap<AssemblyVertex, HashSet<AssemblyVertex>> transitivePrevs;
	int globalcount = 0;

	public LayoutGraphFilter(AssemblyGraph Graph) {
		System.out.println("	graph filter");
		long ini = System.currentTimeMillis();
		this.Graph = Graph;
		procesing = new HashSet<>();
		finished = new HashSet<>();
		prev = new HashMap<>();
		transitivePrevs = new HashMap<>();

		System.out.println("		# edges: " + Graph.getEdges().size());

		for (AssemblyVertex vert : Graph.getVertices()) {
			if (!finished.contains(vert)) {
				transitivePrevs.put(vert, new HashSet<>());
				findTrasition("", vert);
			}
		}
		System.out.println("		# edges: " + (Graph.getEdges().size() - globalcount));
		System.out.println("	graph filter: " + (System.currentTimeMillis() - ini) / (double) 1000 + " s");
	}

	public void findTrasition(String s, AssemblyVertex vert) {
//		System.out.println("		::" + s + vert.getIndex());
//		System.out.println("		::" + s + "prev: " + ((prev.get(vert) != null) ? prev.get(vert).getIndex() : null));
//
//		System.out.println("		::" + s + "Tprevs: "
//				+ Arrays.toString(transitivePrevs.get(vert).stream().mapToInt((x) -> x.getIndex()).toArray()));
		procesing.add(vert);

		AssemblyVertex comp = Graph.getVertex(vert.getIndex() / 2, !vert.isStart());
		Collection<AssemblyVertex> ady = getAdyacentSortedVertex(comp, vert);
		for (AssemblyVertex next : ady) {
			if (prev.containsKey(next)) {
				AssemblyVertex pre = prev.get(next);
				if (transitivePrevs.get(vert).contains(pre)) {
					// IS Transition
					// System.out.println(" transacion detectada");
					System.out.println("		->" + comp.getIndex() + " ---- " + next.getIndex());
					globalcount++;
				}
			}
			prev.put(next, vert);
			HashSet<AssemblyVertex> transitive = transitivePrevs.computeIfAbsent(next,
					(x) -> new HashSet<AssemblyVertex>());
			transitive.addAll(transitivePrevs.get(vert));
			if (prev.get(vert) != null) {
				transitive.add(prev.get(vert));
			}
		}
		for (AssemblyVertex next : ady) {
			if (!finished.contains(next))
				findTrasition(s + "---", next);
		}

		finished.add(vert);
	}

	public Collection<AssemblyVertex> getAdyacentSortedVertex(AssemblyVertex vert, AssemblyVertex non) {
		ArrayList<AssemblyVertex> ans = new ArrayList<AssemblyVertex>(vert.getEdges().size());
		PriorityQueue<Entry<AssemblyVertex, Integer>> queue = new PriorityQueue<>(
				(a, b) -> b.getValue() - a.getValue());
		for (AssemblyEdge edge : vert.getEdges()) {
			AssemblyVertex next = (vert == edge.getVertex1()) ? edge.getVertex2() : edge.getVertex1();
			if (!procesing.contains(next) && next != non)
				queue.add(new AbstractMap.SimpleEntry<>(next, edge.getOverlap()));
		}
		while (!queue.isEmpty()) {
			ans.add(queue.poll().getKey());
		}
		return ans;
	}

}
