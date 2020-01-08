package ngsep.genome;

public class Transposon extends GenomicRegionImpl{
	
	private String type;
	
	private int score;
	
	public Transposon(String sequenceName, int first, int last, String ptype, int pscore) {
		super(sequenceName, first, last);
		type = ptype;
		score = pscore;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}
	
	public int getScore() {
		return score;
	}
	
	public void setScore(int score) {
		this.score = score;
	}
	
}
