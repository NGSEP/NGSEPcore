package ngsep.assembly;

import java.io.Serializable;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class AssemblyVertex implements Serializable {
	private CharSequence read;
	private boolean start;
	private int index;
	private Map<Integer, AssemblyEdge> connectedVertices = new HashMap<>();

	public AssemblyVertex(CharSequence read, boolean start, int index) {
		this.read = read;
		this.start = start;
		this.index = index;
	}

	/**
	 * @return the read
	 */
	public CharSequence getRead() {
		return read;
	}

	/**
	 * @return the start
	 */
	public boolean isStart() {
		return start;
	}
	public int getIndex () {
		return index;
	}
	
	public Collection<AssemblyEdge> getEdges () {
		return connectedVertices.values();
	}
	public boolean isConnected (AssemblyVertex vertex) {
		return connectedVertices.containsKey(vertex.getIndex());
	}
	
	public void removeEdge(AssemblyEdge edge) {
		AssemblyVertex v1 = edge.getVertex1();
		AssemblyVertex v2 = edge.getVertex2();
		if(v1.getIndex()==index) connectedVertices.remove(v2.getIndex());
		else connectedVertices.remove(v1.getIndex());
	}
	public void addEdge (AssemblyEdge edge) {
		AssemblyVertex v1 = edge.getVertex1();
		AssemblyVertex v2 = edge.getVertex2();
		if(v1.getIndex()==index) connectedVertices.put(v2.getIndex(),edge);
		else connectedVertices.put(v1.getIndex(),edge);
	}
	public AssemblyEdge getEdge (AssemblyVertex vertex) {
		return connectedVertices.get(vertex.getIndex());
	}
	public AssemblyVertex getConnectedVertex (int index) {
		AssemblyEdge edge = connectedVertices.get(index);
		if(edge == null) return null;
		AssemblyVertex v1 = edge.getVertex1();
		AssemblyVertex v2 = edge.getVertex2();
		if(v1.getIndex()==index) return v1;
		return v2;
	}
	
	
	
	
}

