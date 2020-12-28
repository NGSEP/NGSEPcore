package ngsep.genome;

public class SyntenyEdge {
	
	private HomologyEdge source;
	
	private HomologyEdge target;
	
	private int weight;

	public SyntenyEdge(HomologyEdge source, HomologyEdge target, int weight) {
		super();
		this.source = source;
		this.target = target;
		this.weight = weight;
	}

	public HomologyEdge getSource() {
		return source;
	}

	public HomologyEdge getTarget() {
		return target;
	}

	public int getWeight() {
		return weight;
	}
	

}
