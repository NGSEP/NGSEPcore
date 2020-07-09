package ngsep.graphs;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class SparseVector {
	private int length = 0;  
	private HashMap<Integer, Double> vector;
	
	public SparseVector(int length) {
		vector = new HashMap<>();
		this.length = length;
	}
	
	public void set(int pos, double k) {
		if (pos < 0 || pos >= length) {
			throw new IndexOutOfBoundsException();
		} else {
			vector.put(pos, k);
			length++;
		}
	}
	
	public double get(int pos) {
		if (pos < 0 || pos >= length) {
			throw new IndexOutOfBoundsException();
		} else {
			if (vector.containsKey(pos)) {
				return vector.get(pos);
			} else {
				return 0;
			}
		}
	}
	
	public void append(double k) {
		vector.put(length++, k);
	}
	
	public List<ValuePair> valuePairs() {
		ArrayList<ValuePair> tuples = new ArrayList<ValuePair>();
		for (int index : vector.keySet()) {
			tuples.add(new ValuePair(index, vector.get(index)));
		}
		return tuples;
	}
	
	class ValuePair {
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
}