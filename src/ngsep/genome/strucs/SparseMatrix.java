package ngsep.genome.strucs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class SparseMatrix {
	private SparseVector[] vectors;
	private int length;
	private int height;
	
	public SparseMatrix(int length, int height) {
		vectors = new SparseVector[length];
		for (int i = 0; i < length; i++) vectors[i] = new SparseVector(height);
		this.length = length;
		this.height = height;
	}
	
	public double get(int i, int j) {
		return vectors[i].get(j);
	}
	
	public void set(int i, int j, double k) {
		vectors[i].set(j, k);
	}
	
	public double sumOfRow(int i) {
		double sum = 0;
		List<ValuePair> values = vectors[i].valuePairs();
		for(ValuePair p : values) {
			sum += p.value;
		}
		return sum;
	}
	
	public List<ValuePair> getRowAsTuples(int i) {
		return vectors[i].valuePairs();
	}
	
	public double sumOfColumn(int j) {
		double sum = 0;
		for(int i = 0; i < length; i++) {
			sum += this.get(i, j);
		}
		return sum;
	}
	
	public int length() {
		return length;
	}
	
	public int height() {
		return height;
	}
	
	public String getLineAsString(int i) {
		String line = "[";
		for (int j = 0; j < height; j++) {
			line = line.concat(String.format(" %f,", this.get(i, j)));
		}
		line = line.substring(0,line.length() - 1);
		line = line.concat("]");
		return line;
	}
	
	public Collection<String> getMatrixAsString() {
		ArrayList<String> matrix = new ArrayList<String>();
		for (int i = 0; i < length; i++) {
			List<ValuePair> tuples = vectors[i].valuePairs();
			if (tuples.size() != 0) {
				String line = String.format("Row %d: %s", i, Arrays.toString(tuples.toArray()));
				matrix.add(line);
			}
		}
		return matrix;
	}
	
	public static void main(String[] args) {
		SparseMatrix m = new SparseMatrix(5, 5);
		m.set(1, 4, 10.47);
		System.out.println(m.getMatrixAsString());
	}
}
