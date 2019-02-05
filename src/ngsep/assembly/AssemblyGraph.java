package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class AssemblyGraph {

	private List<CharSequence> sequences;
	private List<AssemblyVertex> vertices;
	private List<AssemblyEdge> edges = new ArrayList<>();
	
	public AssemblyGraph (List<CharSequence> sequences) {
		this.sequences = Collections.unmodifiableList(sequences);
		vertices = new ArrayList<>();
		edges = new ArrayList<>();
		for(CharSequence seq:sequences) {
			AssemblyVertex vS = new AssemblyVertex(seq, true); 
			vertices.add(vS);
			AssemblyVertex vE = new AssemblyVertex(seq, false);
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

	public void addEdge (AssemblyVertex v1, AssemblyVertex v2, int overlap) {
		edges.add(new AssemblyEdge(v1, v2, overlap));
	}
	
	public AssemblyVertex getVertex (int indexSequence, boolean start) {
		return vertices.get(2*indexSequence+(start?0:1));
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
	
	
}
