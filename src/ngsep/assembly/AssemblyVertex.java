package ngsep.assembly;

public class AssemblyVertex {
	private CharSequence read;
	private boolean start;
	public AssemblyVertex(CharSequence read, boolean start) {
		super();
		this.read = read;
		this.start = start;
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
}
