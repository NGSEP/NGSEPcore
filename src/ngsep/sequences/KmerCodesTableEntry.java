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
	public KmerCodesTableEntry(long kmerCode, long entryCode) {
		this.kmerCode = kmerCode;
		int [] dec = decode(entryCode);
		sequenceId = dec[0];
		start = dec[1];
	}
	/**
	 * Encodes the information of this object in a long number. The minimizer number is not saved
	 * @return long number with the sequence id in the upper 32 bits and the start position in the lower 32 bits
	 */
	public long encode() {
		return encode(sequenceId,start);
	}

	public static int [] decode(long entryCode) {
		int start = (int) (entryCode & 0xFFFFFFFF);
		int sequenceId = (int) (entryCode >> 32);
		int [] answer = {sequenceId,start};
		return answer;
	}

	/**
	 * Encodes the information of this object in a long number. The minimizer number is not saved
	 * @return long number with the sequence id in the upper 32 bits and the start position in the lower 32 bits
	 */
	public static long encode(int sequenceId, int start) {
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
