package ngsep.sequences;

public class MinimizersTableEntry {
	private int minimizer;
	private int sequenceId;
	private int start;
	
	
	public MinimizersTableEntry(int minimizer, int sequenceId, int start) {
		this.minimizer = minimizer;
		this.sequenceId = sequenceId;
		this.start = start;
	}

	public MinimizersTableEntry (int minimizer, long entryCode) {
		this.minimizer = minimizer;
		start = (int) (entryCode & 0xFFFFFFFF);
		sequenceId = (int) (entryCode >> 32);
	}

	/**
	 * Encodes the information of this object in a long number. The minimizer number is not saved
	 * @return long number with the sequence id in the upper 32 bits and the start position in the lower 32 bits
	 */
	public long encode() {
		long code = sequenceId;
		code = code << 32;
		code+=start;
		return code;
	}

	public int getMinimizer() {
		return minimizer;
	}

	public int getSequenceId() {
		return sequenceId;
	}
	public int getStart() {
		return start;
	}
}
