package ngsep.sequences;

public class KmerCodesTableEntry {
	private long kmerCode;
	private int sequenceId;
	private int start;
	
	
	public KmerCodesTableEntry(long kmerCode, int sequenceId, int start) {
		this.kmerCode = kmerCode;
		this.sequenceId = sequenceId;
		this.start = start;
	}

	public KmerCodesTableEntry (long kmerCode, long entryCode) {
		this.kmerCode = kmerCode;
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

	public long getKmerCode() {
		return kmerCode;
	}

	public int getSequenceId() {
		return sequenceId;
	}
	public int getStart() {
		return start;
	}
}
