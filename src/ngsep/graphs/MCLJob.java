package ngsep.graphs;

import java.util.Arrays;
import java.util.List;

public class MCLJob extends Thread {
	private static final float FAKE_LINK_VALUE = 20;
	private static final float ZERO_THRESHOLD = 0.0001f;
	private static final float INFLATION_COEFFICIENT = 2;
	
	private float[][] similarityMatrix;
	
	public MCLJob(float[][] similarityMatrix) {
		super();
		this.similarityMatrix = similarityMatrix;
	}

	@Override
	public void run() {
		similarityMatrix = completeMatrix(similarityMatrix);
		similarityMatrix = normalizeMatrix(similarityMatrix);
	}
	
	public float[][] getSimilarityMatrix() {
		return similarityMatrix;
	}
	
	/**
	 * Normalizes the matrix given as input, also adds self loops equal to 1/n, where n is the quantity of viable neighbours for that node
	 * @param matrix matrix to be normalized
	 * @return
	 */
	private float[][] normalizeMatrix(float[][] matrix) {
		for(int i = 0; i < matrix.length; i++) {
			float sum = 0;
			int valid = 0;
			for(int j = 0; j < matrix.length; j++) {
				if (matrix[i][j] > 0) {
					sum += matrix[i][j];
					valid++;
				}
			}
			
			if(sum > 0) {
				matrix[i][i] = sum/valid;
				sum += (sum/valid);
				for(int j = 0; j < matrix.length; j++) matrix[i][j] = matrix[i][j]/sum;
			} else {
				matrix[i][i] = 1;
			}
		}
		return matrix;
	}
	
	/**
	 * Returns a matrix where if node_i has a link to node_j, then node_j will have a link to node_i. The added links have the value given by constant FAKE_LINK_VALUE
	 * @param matrix matrix to be completed
	 * @return
	 */
	private float[][] completeMatrix(float[][] matrix) {
		for(int i = 0; i < matrix.length; i++) {
			for(int j = 0; j < matrix.length - i; j++) {
				if(i == j) continue;
				
				if(matrix[i][j] > 0 && matrix[j][i] == 0) {
					matrix[j][i] = FAKE_LINK_VALUE;
				} else if (matrix[j][i] > 0 && matrix[i][j] == 0) {
					matrix[i][j] = FAKE_LINK_VALUE;
				}
			}
		}
		return matrix;
	}
	
	/**
	 * Squares the provided matrix a given number of times
	 * @param matrix matrix to be squared, this matrix is square (n x n)
	 * @param times iterations to square the matrix, has to be > 0. e.g. times = 4 returns M^16
	 * @return
	 */
	private float[][] squareMatrixTimes(float[][] matrix, int times) {
		if(times < 1) {
			return matrix;
		}
		
		for(int k = 0; k < times; k++) {
			float[][] base = new float[matrix.length][];
			for(int i = 0; i < matrix.length; i++) base[i] = new float[matrix.length];
			
			for(int i = 0; i < matrix.length; i++) {
				for(int j = 0; j < matrix.length; j++) {
					float val = 0;
					for(int z = 0; z < matrix.length; z++) {
						val += matrix[i][z]*matrix[z][j];
					}
					base[i][j] = val;
				}
			}
			
			matrix = base;
		}
		
		return matrix;
	}
	
	/**
	 * Takes the provided matrix and multiplies every value by the given coefficient.
	 * @param matrix matrix to be inflated.
	 * @param coefficient coefficient to be used.
	 * @return
	 */
	private float[][] inflateMatrix(float[][] matrix, float coefficient) {
		for(int i = 0; i < matrix.length; i++) {
			for(int j = 0; j < matrix.length; j++) {
				matrix[i][j] = matrix[i][j]*coefficient;
			}
		}
		return matrix;
	}
	
	/**
	 * Generates a copy of the given matrix.
	 * @param matrix matrix to be copied.
	 * @return
	 */
	private float[][] copyOf(float[][] matrix) {
		float[][] copy = new float[matrix.length][];
		for(int i = 0; i < matrix.length; i++) copy[i] = Arrays.copyOf(matrix[i], matrix[i].length);
		return copy;
	}
	
	public static void main(String[] args) {
		float[][] test = {{1,2,3},{4,5,6},{7,8,9}};
		MCLJob job = new MCLJob(test);
		test = job.squareMatrixTimes(test, 3);
		
		for(float[] arr : test) {
			System.out.println(Arrays.toString(arr));
		}
	}
}
