package ngsep.sequences;

import java.io.Serializable;

public class UngappedSearchHit implements Serializable {
	private static final long serialVersionUID = 5920964399280119643L;
	private int queryIdx=-1;
	private CharSequence query;
	private int sequenceIdx=-1;
	private String sequenceName;
	//Zero based start of this hit within the sequence
	private int start;
	//Weight of this hit encoded as a byte
	private byte weight = 100;
	
	/**
	 * Creates a new ungapped hit to a subject sequence
	 * @param query sequence
	 * @param sequenceIdx id of the subject sequence
	 * @param start zero based first position of the hit within the subject
	 */
	public UngappedSearchHit(CharSequence query, int sequenceIdx, int start) {
		this.query = query;
		this.sequenceIdx = sequenceIdx;
		this.start = start;
	}
	public CharSequence getQuery() {
		return query;
	}
	public int getQueryIdx() {
		return queryIdx;
	}
	
	public void setQueryIdx(int queryIdx) {
		this.queryIdx = queryIdx;
	}
	public int getSequenceIdx() {
		return sequenceIdx;
	}
	public String getSequenceName() {
		return sequenceName;
	}
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
	
	public int getStart() {
		return start;
	}
	public void setWeight (double weight) {
		if(weight>1) weight=1;
		if(weight<0) weight = 0;
		this.weight = (byte)Math.round(100.0*weight);
	}
	public double getWeight () {
		return 0.01*weight;
	}
}
