package ngsep.assembly;

public class AssembyEmbedded {
	private CharSequence read;
	private int startPosition;
	private boolean isReverse;

	public AssembyEmbedded(CharSequence read, int startPosition, boolean isReverse) {
		this.read = read;
		this.startPosition = startPosition;
		this.isReverse = isReverse;
	}

	public CharSequence getRead() {
		return read;
	}

	public void setRead(CharSequence read) {
		this.read = read;
	}

	public int getStartPosition() {
		return startPosition;
	}

	public void setStartPosition(int startPosition) {
		this.startPosition = startPosition;
	}

	public boolean isReverse() {
		return isReverse;
	}

	public void setReverse(boolean isReverse) {
		this.isReverse = isReverse;
	}
}
