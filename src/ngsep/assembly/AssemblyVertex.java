package ngsep.assembly;

public class AssemblyVertex {
	private CharSequence read;
	private boolean start;
	private int index;

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
}
