package ngsep.genome.strucs;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;

import com.sun.xml.internal.bind.v2.schemagen.xmlschema.List; 

public class SparseMatrix {
	private SparseVector[] vectors;
	
	public SparseMatrix(int length, int height) {
		vectors = new SparseVector[length];
		for (int i = 0; i < length; i++) vectors[i] = new SparseVector(height);
	}
	
	public double get(int i, int j) {
		return vectors[i].get(j);
	}
	
	public void set(int i, int j, double k) {
		vectors[i].set(j, k);
	}
	
	public String getLineAsString(int i) {
		String line = "[";
		for (int j = 0; j < vectors.length; j++) {
			line = line.concat(String.format(" %f,", this.get(i, j)));
		}
		line = line.substring(0,line.length() - 1);
		line = line.concat("]");
		return line;
	}
	
	public Collection<String> getMatrixAsString() {
		ArrayList<String> matrix = new ArrayList();
		for (int i = 0; i < vectors.length; i++) {
			matrix.add(this.getLineAsString(i));
		}
		return matrix;
	}
	
	class SparseVector {
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
	}
	
	public static void main(String[] args) {
		SparseMatrix m = new SparseMatrix(5, 5);
		m.set(1, 4, 10.47);
		System.out.println(m.getMatrixAsString());
	}
}
