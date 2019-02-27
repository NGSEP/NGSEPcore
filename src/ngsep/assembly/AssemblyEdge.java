package ngsep.assembly;

public class AssemblyEdge {
	private AssemblyVertex vertex1;
	private AssemblyVertex vertex2;
	private int overlap;

	public AssemblyEdge(AssemblyVertex vertex1, AssemblyVertex vertex2, int overlap) {
		this.vertex1 = vertex1;
		this.vertex2 = vertex2;
		this.overlap = overlap;
	}

	/**
	 * @return the vertex1
	 */
	public AssemblyVertex getVertex1() {
		return vertex1;
	}

	/**
	 * @return the vertex2
	 */
	public AssemblyVertex getVertex2() {
		return vertex2;
	}

	/**
	 * @return the overlap
	 */
	public int getOverlap() {
		return overlap;
	}

}
