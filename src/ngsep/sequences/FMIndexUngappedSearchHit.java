package ngsep.sequences;

public class FMIndexUngappedSearchHit {
	private int queryIdx=-1;
	private String query;
	private int sequenceIdx=-1;
	private String sequenceName;
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
	public FMIndexUngappedSearchHit(String query, String sequenceName, int start) {
		this.query = query;
		this.sequenceName = sequenceName;
		this.start = start;
	}
	/**
	 * Creates a new ungapped hit to an FMIndex sequence
	 * @param query sequence
	 * @param sequenceIdx id of the subject sequence
	 * @param sequenceName Name of the subject sequence
	 * @param start zero based first position of the hit within the subject
	 */
	public FMIndexUngappedSearchHit(String query, int sequenceIdx, String sequenceName, int start) {
		this.query = query;
		this.sequenceIdx = sequenceIdx;
		this.sequenceName = sequenceName;
		this.start = start;
	}
	public String getQuery() {
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
