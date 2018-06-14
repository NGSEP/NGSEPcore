package ngsep.assembly;

import java.io.InputStream;
import java.io.PrintStream;
import java.util.AbstractMap.SimpleEntry;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.stream.Stream;
import java.util.Scanner;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

public class GraphSimplificator {
    /** #Vertices Original Graph **/
    protected int N;
    /** Edges with weight **/
    protected Map<Integer, Map<Integer, Integer>> Edges;
    /** Vertices of reducition **/
    private Map<Integer, List<Integer>> Verts;
    // ** Weight of the reduction **/
    private Map<Integer, Integer> weigth;
    
  

    public GraphSimplificator(int N, Map<Integer, Map<Integer, Integer>> Edges) {
	this.N = N;
	this.Edges = Edges;
	Initialize();
    }

    public GraphSimplificator(InputStream inputStream) {
	try (Scanner sc = new Scanner(inputStream)) {
	    N = Integer.parseInt(sc.nextLine());
	    Edges = new HashMap<>();

	    for (int i = 0; i < N; i++)
		Edges.put(i, new HashMap<>());

	    int[] aux;
	    while (true) {
		String line = sc.nextLine();
		if (line.isEmpty())
		    break;
		aux = Stream.of(line.split(",")).mapToInt((x) -> Integer.parseInt(x)).toArray();
		Edges.get(aux[0]).put(aux[2], aux[1]);
	    }
	    Initialize();
	} catch (Exception e) {
	    e.printStackTrace();
	}
    }

    private void Initialize() {
	//printWhitFormat(System.out);
	System.out.println(roots());
	removeCicles();
	int[] top = getTopologicalOrder();
	removeTransitiveEdges(top);
	printWhitFormat(System.out);
	System.out.println(roots());
    }

    private List<Integer> roots() {
	int[] top = getTopologicalOrder();
	boolean[] nop = new boolean[N];
	for (Integer i : top) {
	    if (Edges.get(i).isEmpty())
		nop[i] = true;
	    for (Integer j : Edges.get(i).keySet())
		nop[j] = true;
	}
	List<Integer> ans = new LinkedList<>();

	for (int i = 0; i < N; i++) {
	    if (!nop[i])
		ans.add(i);
	}

	return ans;
    }
    

    private void reduce(int[] top, Map<Integer, Set<Integer>> foward) {
	boolean[] moreThanOne = new boolean[N];
	for (int i = 0; i < N; i++)
	    moreThanOne[i] = Edges.get(i).size() > 1 || foward.get(i).size() > 1;

	int[] map = new int[N];
	Arrays.fill(map, -1);
	boolean[] mark = new boolean[N];
	Verts = new Hashtable<>();
	weigth = new Hashtable<>();
	List<Integer> aux;
	int prev, sum, Node = top[N - 1];
	int N = 0;

	do {
	    int firts = Node;
	    mark[Node] = true;
	    aux = new LinkedList<>();
	    aux.add(Node);
	    sum = 0;
	    if (!moreThanOne[Node] && foward.get(Node).iterator().hasNext()) {
		prev = Node;
		Node = foward.get(Node).iterator().next();
		while (!moreThanOne[Node] && foward.get(Node).iterator().hasNext()) {
		    mark[Node] = true;
		    aux.add(0, Node);
		    sum += Edges.get(Node).get(prev);
		    prev = Node;
		    Node = foward.get(Node).iterator().next();
		}
		if (firts != prev) {
		    Edges.get(prev).clear();
		    Edges.get(prev).putAll(Edges.get(firts));
		    Iterator<Integer> iter = aux.iterator();
		    iter.next();
		    while (iter.hasNext())
			Edges.put(iter.next(), new HashMap<>());
		}
		Node = prev;
	    }
	    Verts.put(Node, aux);
	    weigth.put(Node, sum);
	    map[Node] = N++;
	    Node = NextNodeInv(mark);
	} while (Node != -1);

	Map<Integer, Map<Integer, Integer>> EdgesAux = new Hashtable<>();
	Map<Integer, List<Integer>> VertsAux = new Hashtable<>();
	Map<Integer, Integer> weigthAux = new Hashtable<>();

	for (Entry<Integer, List<Integer>> entry : Verts.entrySet())
	    VertsAux.put(map[entry.getKey()], entry.getValue());
	Verts = VertsAux;
	for (Entry<Integer, Integer> entry : weigth.entrySet())
	    weigthAux.put(map[entry.getKey()], entry.getValue());
	weigth = weigthAux;
	for (Entry<Integer, Map<Integer, Integer>> entry : Edges.entrySet()) {
	    if (map[entry.getKey()] != -1) {
		weigthAux = new Hashtable<>();
		for (Entry<Integer, Integer> entry2 : entry.getValue().entrySet())
		    weigthAux.put(map[entry2.getKey()], entry2.getValue());
		EdgesAux.put(map[entry.getKey()], weigthAux);
	    }
	}
	this.N = N;
	Edges = EdgesAux;

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

    private int NextNodeInv(boolean[] Done) {
	for (int i = Done.length - 1; i >= 0; i--)
	    if (!Done[i])
		return i;
	return -1;
    }

    private Map<Integer, Set<Integer>> removeTransitiveEdges(int[] top) {
	int[] invTop = new int[top.length];
	for (int i = 0; i < top.length; i++)
	    invTop[top[i]] = i;
	Set<Integer> prevs;

	Map<Integer, Set<Integer>> prevMap = new HashMap<>(N);
	for (Integer integer : Edges.keySet())
	    prevMap.put(integer, new TreeSet<>());

	for (Integer Vert : top) {
	    Map<Integer, Integer> edges = Edges.get(Vert);
	    for (Integer toVert : edges.keySet()) {
		prevs = prevMap.get(toVert);
		if (prevs.isEmpty())
		    prevs.add(Vert);
		else {
		    Iterator<Integer> iter = prevs.iterator();
		    while (iter.hasNext()) {
			Integer otherEdge = iter.next();
			if (search(otherEdge, Vert, prevMap, 15, invTop, invTop[otherEdge])) {
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

    private boolean search(Integer searched, int from, Map<Integer, Set<Integer>> prevMap, int max, int[] top,
	    int topLimit) {
	if (max == 0 || top[from] < topLimit)
	    return false;

	Set<Integer> prevs = prevMap.get(from);
	for (Integer integer : prevs)
	    if (integer == searched)
		return true;

	for (Integer integer : prevs)
	    if (search(searched, integer, prevMap, max - 1, top, topLimit))
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

    public static void main(String[] args) {
	new GraphSimplificator(System.in);
    }
}
