package ngsep.genome;

public class TransposableElementAnnotation implements GenomicRegion{
	private String sequenceName;
	private int first;
	private int last;
	private String taxonomy;
	
	
	public TransposableElementAnnotation(String sequenceName, int first, int last) {
		super();
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
	}
	public String getTaxonomy() {
		return taxonomy;
	}
	public void setTaxonomy(String taxonomy) {
		this.taxonomy = taxonomy;
	}
	public String getSequenceName() {
		return sequenceName;
	}
	public int getFirst() {
		return first;
	}
	public int getLast() {
		return last;
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return true;
	}
	@Override
	public boolean isNegativeStrand() {
		return false;
	}
}
