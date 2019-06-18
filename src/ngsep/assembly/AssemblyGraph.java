package ngsep.assembly;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.TreeMap;
import ngsep.sequences.DNAMaskedSequence;

public class AssemblyGraph implements Serializable {

	private List<CharSequence> sequences;
	private List<AssemblyVertex> vertices;

	private List<AssemblyEdge> edges = new ArrayList<>();
	private Map<Integer, List<AssemblyEmbedded>> embeddedSequences = new HashMap<>();

	// Indexes in the vertices list
	private List<List<AssemblyEdge>> paths = new ArrayList<List<AssemblyEdge>>();

	public AssemblyGraph(List<CharSequence> sequences) {
		this.sequences = Collections.unmodifiableList(sequences);
		vertices = new ArrayList<>();
		edges = new ArrayList<>();
		for (CharSequence seq : sequences) {
			AssemblyVertex vS = new AssemblyVertex(seq, true, vertices.size());
			vertices.add(vS);
			AssemblyVertex vE = new AssemblyVertex(seq, false, vertices.size());
			vertices.add(vE);
			addEdge(vS, vE, seq.length());
		}
	}

	/**
	 * @return the sequences
	 */
	public List<CharSequence> getSequences() {
		return sequences;
	}

	public void addEdge(AssemblyVertex v1, AssemblyVertex v2, int overlap) {
		AssemblyEdge edge = new AssemblyEdge(v1, v2, overlap); 
		edges.add(edge);
		v1.addEdge(edge);
		v2.addEdge(edge);
	}

	public AssemblyVertex getVertex(int indexSequence, boolean start) {
		return vertices.get(2 * indexSequence + (start ? 0 : 1));
	}

	public void addEmbedded(int ind, AssemblyEmbedded embedded) {
		List<AssemblyEmbedded> list = embeddedSequences.computeIfAbsent(ind, key -> new LinkedList<>());
		list.add(embedded);
	}

	/**
	 * Return the list of embedded sequences for the given read
	 * 
	 * @param index of the read
	 * @return list of embedded sequences
	 */
	public List<AssemblyEmbedded> getEmbedded(int index) {
		return embeddedSequences.get(index);
	}

	public void addPath(List<AssemblyEdge> path) {
		paths.add(path);
	}

	/**
	 * @return the vertices
	 */
	public List<AssemblyVertex> getVertices() {
		return vertices;
	}

	/**
	 * @return the edges
	 */
	public List<AssemblyEdge> getEdges() {
		return edges;
	}

	/**
	 * @return the paths
	 */
	public List<List<AssemblyEdge>> getPaths() {
		return paths;
	}

	public List<Map<Integer, CharSequence>> paths() {
		// calcular los vertices correspondientes
		List<Map<Integer, CharSequence>>  ans = new ArrayList<>();
		for (int i = 0; i < paths.size(); i++) {
			List<AssemblyVertex> vertx = new ArrayList<AssemblyVertex>();
			List<AssemblyEdge> list = paths.get(i);

			AssemblyVertex v1 = list.get(0).getVertex1();
			AssemblyVertex v2 = list.get(0).getVertex2();
			AssemblyVertex v3 = list.get(1).getVertex1();
			AssemblyVertex v4 = list.get(1).getVertex2();

			AssemblyVertex vact;
			if (v1 == v3) {
				vertx.add(v2);
				vertx.add(v1);
				vact = v4;
			} else if (v1 == v4) {
				vertx.add(v2);
				vertx.add(v1);
				vact = v3;
			} else if (v2 == v3) {
				vertx.add(v1);
				vertx.add(v2);
				vact = v4;
			} else {
				vertx.add(v1);
				vertx.add(v2);
				vact = v3;
			}

			for (int j = 2; j < list.size(); j++) {
				vertx.add(vact);
				v3 = list.get(j).getVertex1();
				v4 = list.get(j).getVertex2();
				vact = (vact == v3) ? v4 : v3;
			}
			vertx.add(vact);
			
			System.out.println( Arrays.toString(vertx.stream().mapToInt((x) -> x.getIndex()).toArray()));

	
			TreeMap<Integer,CharSequence> ord = new TreeMap<>();
			int llll = 0;
			for (int j = 0; j < vertx.size() - 2; j += 2) {
				CharSequence str = vertx.get(j).getRead();
				boolean rev = false;
				if (vertx.get(j).getIndex() > vertx.get(j + 1).getIndex()) {
					str = DNAMaskedSequence.getReverseComplement(str);
					rev = true;
				}
				ord.put(llll,str);

				if (embeddedSequences.get(vertx.get(j).getIndex() / 2) != null) {
					for (AssemblyEmbedded emb : embeddedSequences.get(vertx.get(j).getIndex() / 2)) {
						str = emb.getRead();
						int ill = emb.getStartPosition();
						if (rev ^ emb.isReverse()) {
							str = DNAMaskedSequence.getReverseComplement(str);
							ill = list.get(j).getOverlap() - ill - emb.getRead().length();
						}
						ord.put( llll + ill,str);
					}
				}

				llll += list.get(j).getOverlap() - list.get(j + 1).getOverlap();
			}

			int j = vertx.size() - 2;
			CharSequence str = vertx.get(vertx.size() - 1).getRead();
			boolean rev = false;
			if (vertx.get(vertx.size() - 2).getIndex() > vertx.get(vertx.size() - 1).getIndex()) {
				str = DNAMaskedSequence.getReverseComplement(str);
				rev = true;
			
			}
			ord.put(llll,str);

			if (embeddedSequences.get(vertx.get(j).getIndex() / 2) != null) {
				for (AssemblyEmbedded emb : embeddedSequences.get(vertx.get(j).getIndex() / 2)) {
					str = emb.getRead();
					int ill = emb.getStartPosition();
					if (rev ^ emb.isReverse()) {
						str = DNAMaskedSequence.getReverseComplement(str);
						ill = list.get(j).getOverlap() - ill - emb.getRead().length();
					}
					ord.put( llll + ill,str);
				}
			}
			System.out.println(ord.size());
			
			ans.add(ord);
			print(ord);
		}
		return ans;
	}
	
	private void print(Map<Integer, CharSequence> ord) {
		PriorityQueue<StringBuilder> prints = new PriorityQueue<>(new Comparator<StringBuilder>() {
			public int compare(StringBuilder o1, StringBuilder o2) {
				return o1.length()-o2.length();
			}
		});
		prints.add(new StringBuilder());
		StringBuilder trs = null;
		for(Entry<Integer,CharSequence> ent: ord.entrySet()) {
			if(prints.peek().length() <= ent.getKey()) {
				trs = prints.poll();
			}else {
				trs = new StringBuilder();
			}
			int id = ent.getKey() - trs.length();
			for(int k =0;k<id;k++)
				trs.append('-');
			trs.append(ent.getValue());
			prints.add(trs);
			
		}
		
		while(!prints.isEmpty())
			System.out.println(prints.poll());
	}
}
