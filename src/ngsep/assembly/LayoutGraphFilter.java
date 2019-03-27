package ngsep.assembly;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.PriorityQueue;

import com.sun.javafx.geom.Edge;

public class LayoutGraphFilter {
	private AssemblyGraph Graph;
	private HashSet<AssemblyVertex> procesing;
	private HashSet<AssemblyVertex> finished;
	private HashMap<AssemblyVertex, AssemblyVertex> prev;
	private HashMap<AssemblyVertex, HashSet<AssemblyVertex>> transitivePrevs;
	int globalcount = 0;

	public LayoutGraphFilter(AssemblyGraph Graph) {
		List<AssemblyVertex> list = Graph.getVertices();
		int ccc = 0;
		for (AssemblyVertex ve : list)
			if (ve.getEdges().size() == 1)
				ccc++;
		System.out.println(ccc);
		
		for(AssemblyEdge edg:Graph.getEdges()) {
			System.out.println(edg.getVertex1().getIndex()+":"+edg.getOverlap()+":"+edg.getVertex2().getIndex());
		}
		
		
//		for (int i = 0; i < list.size(); i++) {
//			AssemblyVertex vert = list.get(i);
//			int[] array = new int[list.size()];
//			for (AssemblyEdge edg : vert.getEdges()) {
//				AssemblyVertex vert2 = (edg.getVertex1().getIndex() == vert.getIndex()) ? edg.getVertex2()
//						: edg.getVertex1();
//				array[vert2.getIndex()] = 1;
//			}
//			StringBuilder str = new StringBuilder();
//			for (int j = 0; j < array.length - 1; j++) {
//				str.append(array[j] + ",");
//			}
//			str.append(array[array.length - 1]);
//			System.out.println(str);
//		}

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

		list = Graph.getVertices();
		ccc = 0;
		for (AssemblyVertex ve : list)
			if (ve.getEdges().size() == 1)
				ccc++;
		System.out.println(ccc);
//		for (int i = 0; i < list.size(); i++) {
//			AssemblyVertex vert = list.get(i);
//			int[] array = new int[list.size()];
//			for (AssemblyEdge edg : vert.getEdges()) {
//				AssemblyVertex vert2 = (edg.getVertex1().getIndex() == vert.getIndex()) ? edg.getVertex2()
//						: edg.getVertex1();
//				array[vert2.getIndex()] = 1;
//			}
//			StringBuilder str = new StringBuilder();
//			for (int j = 0; j < array.length - 1; j++) {
//				str.append(array[j] + ",");
//			}
//			str.append(array[array.length - 1]);
//			System.out.println(str);
//		}
	}

	public void findTrasition(String s, AssemblyVertex vert) {
//		System.out.println("		::" + s + vert.getIndex());
//		System.out.println("		::" + s + "prev: " + ((prev.get(vert) != null) ? prev.get(vert).getIndex() : null));
//
//		System.out.println("		::" + s + "Tprevs: "+ Arrays.toString(transitivePrevs.get(vert).stream().mapToInt((x) -> x.getIndex()).toArray()));
		AssemblyVertex comp = Graph.getVertex(vert.getIndex() / 2, !vert.isStart());

		procesing.add(vert);
		Collection<AssemblyVertex> ady = getAdyacentSortedVertex(comp, vert);
//		System.out.println(
//				"		::" + s + "next" + Arrays.toString(ady.stream().mapToInt((x) -> x.getIndex()).toArray()));
		for (AssemblyVertex next : ady) {
			if (prev.containsKey(next)) {
				AssemblyVertex pre = prev.get(next);
//				System.out.println("		::" + s + "-----: "+next.getIndex()+"()"+pre.getIndex());
				if (transitivePrevs.get(vert).contains(pre) || prev.get(vert) == pre) {
					// IS Transition
					// System.out.println(" transacion detectada");
					AssemblyVertex precom = Graph.getVertex(pre.getIndex() / 2, !pre.isStart());
//					System.out.println("		->" + precom.getIndex() + " ---- " + next.getIndex());
					AssemblyEdge edge = precom.getEdge(next);
					if (edge != null) { // aun no estoy seguro de porque pasa
						precom.removeEdge(edge);
						next.removeEdge(edge);
						globalcount++;
					}
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

		procesing.remove(vert);
		// finished.add(vert);
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
