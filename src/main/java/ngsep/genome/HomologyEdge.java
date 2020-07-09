package ngsep.genome;

public class HomologyEdge {
	private HomologyUnit queryUnit;
	private HomologyUnit subjectUnit;
	private double score;
	private double ks;
	private double ka;
	
	public HomologyEdge(HomologyUnit queryUnit, HomologyUnit subjectUnit, double score) {
		super();
		this.queryUnit = queryUnit;
		this.subjectUnit = subjectUnit;
		this.score = score;
	}

	public HomologyUnit getQueryUnit() {
		return queryUnit;
	}

	public HomologyUnit getSubjectUnit() {
		return subjectUnit;
	}

	public double getScore() {
		return score;
	}
	
	public double getKs() {
		return ks;
	}

	public void setKs(double ks) {
		this.ks = ks;
	}

	public double getKa() {
		return ka;
	}

	public void setKa(double ka) {
		this.ka = ka;
	}

}
