package ngsep.sequences;

import java.io.Serializable;

public class UngappedSearchHit implements Serializable {
	private static final long serialVersionUID = 5920964399280119643L;
	private int queryIdx=-1;
	private CharSequence query;
	private int sequenceIdx=-1;
	private String sequenceName;
	private int sequenceLength;
	//Zero based start of this hit within the sequence
	private int start;
	private int totalHitsQuery = 0;
	/**
	 * Creates a new ungapped hit to an FMIndex sequence
	 * @param query sequence
	 * @param sequenceIdx id of the subject sequence
	 * @param sequenceName Name of the subject sequence
	 * @param start zero based first position of the hit within the subject
	 */
	public UngappedSearchHit(CharSequence query, String sequenceName, int start) {
		this.query = query;
		this.sequenceName = sequenceName;
		this.start = start;
	}
	/**
	 * Creates a new ungapped hit to an FMIndex sequence
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
	public int getSequenceLength() {
		return sequenceLength;
	}
	public void setSequenceLength(int sequenceLength) {
		this.sequenceLength = sequenceLength;
	}
	public int getStart() {
		return start;
	}
	public int getTotalHitsQuery() {
		return totalHitsQuery;
	}
	public void setTotalHitsQuery(int totalHitsQuery) {
		this.totalHitsQuery = totalHitsQuery;
	}
}
