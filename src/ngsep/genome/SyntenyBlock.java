package ngsep.genome;

import java.util.List;

public class SyntenyBlock {

	private List<SyntenyEdge> homologies;
	
	private AnnotatedReferenceGenome genome1;
	
	private AnnotatedReferenceGenome genome2;
	
	private int first;
	
	private int last;
	
	private String sequenceName;
	
	public SyntenyBlock(String sequenceName, int first, int last, AnnotatedReferenceGenome g1, AnnotatedReferenceGenome g2, List<SyntenyEdge> homologies) {
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		this.homologies = homologies;
		this.genome1 = g1;
		this.genome2 = g2;
	}

	public List<SyntenyEdge> getHomologies() {
		return homologies;
	}

	public AnnotatedReferenceGenome getGenome1() {
		return genome1;
	}

	public AnnotatedReferenceGenome getGenome2() {
		return genome2;
	}

	public int getFirst() {
		return first;
	}

	public int getLast() {
		return last;
	}

	public String getSequenceName() {
		return sequenceName;
	}
	
	public int length() {
		return last - first + 1;
	}
	
	
}
