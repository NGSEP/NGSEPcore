package ngsep.assembly;

import java.io.PrintStream;
import java.util.AbstractMap.SimpleEntry;
import java.util.Arrays;
import java.util.Deque;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

public class GraphSimplificator {
    /** #Vertices Original Graph **/
    protected int N;
    /** Edges with weight **/
    protected Map<Integer, Map<Integer, Integer>> Edges;

    public GraphSimplificator(int N, Map<Integer, Map<Integer, Integer>> Edges) {
	this.N = N;
	this.Edges = Edges;
	removeCicles();
	removeTransitiveEdges();
    }

    public void printWhitFormat(PrintStream outputStream) {
	// http://graphonline.ru/en/
	int sum = 0;
	for (Map<Integer, Integer> map : Edges.values())
	    sum += map.size();
	outputStream.print(N + "," + sum);

	outputStream.println("-------");
	int[] intrs = new int[N];
	for (int i = 0; i < N; i++) {
	    Arrays.fill(intrs, 0);
	    Set<Entry<Integer, Integer>> lista = Edges.get(i).entrySet();
	    if (lista != null) {
		for (Entry<Integer, Integer> e : lista)
		    intrs[e.getKey()] = e.getValue();
	    }
	    for (int j = 0; j < N - 1; j++)
		outputStream.print(intrs[j] + ",");
	    outputStream.println(intrs[N - 1]);
	}
	outputStream.println("-------");
    }

    private void removeCicles() {
	Stack<Integer> aux = new Stack<>();
	Stack<Integer> P = new Stack<>();
	Stack<Integer> A = new Stack<>();
	Stack<Integer> Ahat = new Stack<>();
	List<Entry<Integer, Integer>> cycle;

	boolean[] procesing = new boolean[N];
	boolean[] Done = new boolean[N];

	int to, prev, Node = 0;
	do {
	    A.push(Node);
	    Ahat.push(-1);
	    Stack: while (!A.isEmpty()) {
		Node = A.pop();
		Ahat.pop();

		boolean nothing = true;
		if (!Done[Node]) {
		    P.push(Node);
		    procesing[Node] = true;

		    for (Integer nextNode : Edges.get(Node).keySet()) {
			if (!Done[nextNode]) {
			    if (procesing[nextNode]) {
				cycle = new LinkedList<>();
				aux.clear();

				to = nextNode;
				prev = P.pop();
				aux.push(prev);
				do {
				    cycle.add(0, new SimpleEntry<>(prev, Edges.get(prev).get(to)));
				    to = prev;
				    if (P.isEmpty())
					break;
				    prev = P.pop();
				    aux.push(prev);
				} while (to != nextNode);

				int removed = removeCicle(cycle);
				while (aux.peek() != removed)
				    P.push(aux.pop());

				for (int integer : aux) {
				    procesing[integer] = false;
				    while (!Ahat.isEmpty() && Ahat.peek() == integer) {
					Ahat.pop();
					A.pop();

				    }
				}

				A.push(removed);
				Ahat.push(-1);

				continue Stack;
			    } else {
				A.push(nextNode);
				Ahat.push(Node);
			    }
			    nothing = false;
			}
		    }

		}
		if (nothing) {
		    if (!Ahat.isEmpty()) {
			to = Ahat.peek();
			while (!P.isEmpty() && P.peek() != to)
			    Done[P.pop()] = true;
		    } else
			while (!P.isEmpty())
			    Done[P.pop()] = true;
		}
	    }

	    Node = NextNode(Done);
	} while (Node != -1);

    }

    private int removeCicle(List<Entry<Integer, Integer>> cycle) {
	int Min = Integer.MAX_VALUE;
	int MinNde = 0;
	int posMin = 0;
	for (int i = 0; i < cycle.size(); i++) {
	    Entry<Integer, Integer> entry = cycle.get(i);
	    if (Min > entry.getValue()) {
		Min = entry.getValue();
		MinNde = entry.getKey();
		posMin = i;
	    }
	}
	int pos = (posMin + 1) % cycle.size();
	pos = cycle.get(pos).getKey();
	Edges.get(MinNde).remove(pos);
	return MinNde;
    }

    private int NextNode(boolean[] Done) {
	for (int i = 0; i < Done.length; i++)
	    if (!Done[i])
		return i;
	return -1;
    }

    private Map<Integer, Set<Integer>> removeTransitiveEdges() {
	int[] top = getTopologicalOrder();
	int[] invTop = new int[top.length];
	for (int i = 0; i < top.length; i++)
	    invTop[top[i]] = i;
	Set<Integer> prevs;

	Map<Integer, Set<Integer>> prevMap = new Hashtable<>(N);
	for (Integer integer : Edges.keySet())
	    prevMap.put(integer, new TreeSet<>());

	for (Integer Vert : top) {
	    for (Integer toVert : Edges.get(Vert).keySet()) {
		prevs = prevMap.get(toVert);
		if (prevs.isEmpty())
		    prevs.add(Vert);
		else {
		    Iterator<Integer> iter = prevs.iterator();
		    while (iter.hasNext()) {
			Integer otherEdge = iter.next();
			if (search(otherEdge, Vert, prevMap, invTop, invTop[otherEdge])) {
			    Edges.get(otherEdge).remove(toVert);
			    iter.remove();
			}
		    }
		    prevs.add(Vert);
		}
	    }
	}
	return prevMap;
    }

    private boolean search(Integer searched, int from, Map<Integer, Set<Integer>> prevMap, int[] top, int topLimit) {
	if (top[from] < topLimit)
	    return false;

	Set<Integer> prevs = prevMap.get(from);
	for (Integer integer : prevs)
	    if (integer == searched)
		return true;

	for (Integer integer : prevs)
	    if (search(searched, integer, prevMap, top, topLimit))
		return true;

	return false;
    }

    protected int[] getTopologicalOrder() {
	boolean[] mark = new boolean[N];

	Deque<Integer> topologicalOrder = new LinkedList<>();
	int nodo = 0;
	do {
	    dfsPosOrden(topologicalOrder, mark, nodo);
	    nodo = NextNode(mark);
	} while (nodo != -1);

	int[] ans = new int[topologicalOrder.size()];
	int i = 0;
	while (!topologicalOrder.isEmpty())
	    ans[i++] = topologicalOrder.poll();
	return ans;
    }

    private void dfsPosOrden(Deque<Integer> order, boolean[] mark, int node) {
	mark[node] = true;
	for (Integer to : Edges.get(node).keySet())
	    if (!mark[to])
		dfsPosOrden(order, mark, to);
	order.push(node);
    }
}
