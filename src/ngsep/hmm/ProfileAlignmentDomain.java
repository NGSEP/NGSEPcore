package ngsep.hmm;

public class ProfileAlignmentDomain {
	private final String query;
	private int start;
	private int end;
	private String hmmID;
	private double evalue;
	private String alignmentQuery;
	private String alignmentProfile;
	public ProfileAlignmentDomain(String query, int start, int end, String hmmID) {
		super();
		this.query = query;
		this.start = start;
		this.end = end;
		this.hmmID = hmmID;
	}
	public double getEvalue() {
		return evalue;
	}
	public void setEvalue(double evalue) {
		this.evalue = evalue;
	}
	public String getQuery() {
		return query;
	}
	public int getStart() {
		return start;
	}
	public int getEnd() {
		return end;
	}
	public String getHmmID() {
		return hmmID;
	}
	
	public String getAlignmentQuery() {
		return alignmentQuery;
	}
	public String getAlignmentProfile() {
		return alignmentProfile;
	}
	public void setAlignment(String alignmentQuery, String alignmentProfile) {
		this.alignmentQuery = alignmentQuery;
		this.alignmentProfile = alignmentProfile;
	}
}
