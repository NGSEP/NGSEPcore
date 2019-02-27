package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class AssemblyGraph {

	private List<CharSequence> sequences;
	private List<AssemblyVertex> vertices;
	private List<AssemblyEdge> edges = new ArrayList<>();
	private Map<Integer, List<AssembyEmbedded>> embeddedSequences = new HashMap<>();
	
	//Indexes in the vertices list
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
			edges.add(new AssemblyEdge(vS, vE, seq.length()));
		}
	}

	/**
	 * @return the sequences
	 */
	public List<CharSequence> getSequences() {
		return sequences;
	}

	public void addEdge(AssemblyVertex v1, AssemblyVertex v2, int overlap) {
		edges.add(new AssemblyEdge(v1, v2, overlap));
	}

	public AssemblyVertex getVertex(int indexSequence, boolean start) {
		return vertices.get(2 * indexSequence + (start ? 0 : 1));
	}

	public void addEmbedded(int ind, AssembyEmbedded embedded) {
		List<AssembyEmbedded> list = embeddedSequences.computeIfAbsent(ind, key -> new LinkedList<>());
		list.add(embedded);
	}
	
	/**
	 * Return the list of embedded sequences for the given read
	 * @param index of the read
	 * @return list of embedded sequences
	 */
	public List<AssembyEmbedded> getEmbedded(int index)
	{
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
}
