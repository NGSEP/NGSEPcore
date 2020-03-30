package ngsep.graphs;

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
}
