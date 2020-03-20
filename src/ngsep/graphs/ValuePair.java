package ngsep.graphs;

public class ValuePair {
	public int index;
	public double value;
	
	public ValuePair(int index, double value) {
		super();
		this.index = index;
		this.value = value;
	}
	
	@Override
	public String toString() {
		return String.format("(%d, %f)", this.index, this.value);
	}
}