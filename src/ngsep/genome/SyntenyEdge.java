package ngsep.genome;

public class SyntenyEdge {
	
	private SyntenyVertex source;
	
	private SyntenyVertex target;
	
	private int weight;

	public SyntenyEdge(SyntenyVertex source, SyntenyVertex target, int weight) {
		super();
		this.source = source;
		this.target = target;
		this.weight = weight;
	}

	public SyntenyVertex getSource() {
		return source;
	}

	public SyntenyVertex getTarget() {
		return target;
	}

	public int getWeight() {
		return weight;
	}
	

}
