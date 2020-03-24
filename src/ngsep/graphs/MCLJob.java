package ngsep.graphs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MCLJob extends Thread {
	private static final double FAKE_LINK_VALUE = 20;
	private static final double ZERO_THRESHOLD = 0.001;
	private static final int E_POWER = 2;
	private static final double INFLATION_COEFFICIENT = 2;
	private static final double DEVIATION_THRESHOLD = 0.000001;
	
	private double[][] matrixCopy;
	private double[][] similarityMatrix;
	private List<List<Integer>> clusters;
	
	public MCLJob(double[][] providedMatrix) {
		super();
		this.similarityMatrix = providedMatrix;
		this.matrixCopy = copyOf(similarityMatrix);
		this.clusters = new ArrayList<>();
	}

	@Override
	public void run() {
		similarityMatrix = completeMatrix(similarityMatrix);
		similarityMatrix = fitMatrix(similarityMatrix);
		
		int runs = 0;
		boolean convergenceState = false;
		double[][] backup;
		while(!convergenceState) {
			//System.out.println(String.format("Runs count %d", ++runs));
			
			backup = copyOf(similarityMatrix);
			backup = squareMatrixTimes(backup, E_POWER);
			backup = inflateMatrix(backup, INFLATION_COEFFICIENT);
			convergenceState = verifyConvergence(backup, similarityMatrix, DEVIATION_THRESHOLD);
			similarityMatrix = backup;
		}
		
		similarityMatrix = consolidateAttractors(similarityMatrix);
		clusters = extractResults(similarityMatrix);
	}

	public double[][] getSimilarityMatrix() {
		return similarityMatrix;
	}
	
	public void printMatrix(double[][] matrix) {
		for(double[] arr : matrix) {
			double sum = 0;
			for(int i = 0; i < arr.length; i++) sum += arr[i];
			System.out.println(String.format("%s || %f", Arrays.toString(arr), sum));
		}
	}
	
	private List<List<Integer>> extractResults(double[][] matrix) {
		List<List<Integer>> clusters = new ArrayList<>();
		
		for(int j = 0; j < matrix.length; j++) {
			ArrayList<Integer> cluster = new ArrayList<>();
			for(int i = 0; i < matrix.length; i++) {
				if(matrix[i][j] > 0) {
					cluster.add(i);
				}
			}
			if(cluster.size() > 0) clusters.add(cluster);
		}
		
		return clusters;
	}
	
	/**
	 * Verifies if the standard deviation the matrix is below the desired threshold. Each matrix needs to have the same dimensions.
	 * @param oldState one state to be compared
	 * @param newState other state to be compared
	 * @return true if the std is below the desired threshold, false otherwise.
	 */
	private boolean verifyConvergence(double[][] oldState, double[][] newState, double threshold) {
		double squaredSum = 0;
		for(int i = 0; i < oldState.length; i++) {
			for(int j = 0; j < oldState[0].length; j++) {
				double errVal = oldState[i][j] - newState[i][j];
				squaredSum += (errVal)*(errVal);
			}
		}
		double count = (oldState.length*oldState[0].length) - 1;
		double std = Math.sqrt(squaredSum/(count));
		//System.out.println(String.format("STD WAS %f", std));
		return std <= threshold;
	}
	
	private double[][] consolidateAttractors(double[][] matrix) {
		for(int i = 0; i < matrix.length; i++) {
			int index = 0;
			double max = 0;
			for(int j = 0; j < matrix.length; j++) {
				if(max < matrix[i][j]) {
					index = j;
					max = matrix[i][j];
				}
				matrix[i][j] = 0;
			}
			matrix[i][index] = 1;
		}
		return matrix;
	}
	
	/**
	 * Normalizes the matrix given as input, also adds self loops equal to 1/n, where n is the quantity of viable neighbours for that node
	 * @param matrix matrix to be normalized
	 * @return
	 */
	private double[][] fitMatrix(double[][] matrix) {
		for(int i = 0; i < matrix.length; i++) {
			double sum = 0;
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
	 * Normalizes a given row, taking the weight of each edge and returning a distribution probability over all edges.
	 * @param row row to be normalized
	 * @return normalized row
	 */
	private double[] normalizeRow(double[] row) {
		double sum = 0;
		for(int j = 0; j < row.length; j++) {
			if (row[j] > 0) {
				sum += row[j];
			}
		}
		
		if(sum > 0) {
			for(int j = 0; j < row.length; j++) row[j] = row[j]/sum;
		}
		
		return row;
	}
	
	/**
	 * Returns a matrix where if node_i has a link to node_j, then node_j will have a link to node_i. The added links have the value given by constant FAKE_LINK_VALUE
	 * @param matrix matrix to be completed
	 * @return
	 */
	private double[][] completeMatrix(double[][] matrix) {
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
	private double[][] squareMatrixTimes(double[][] matrix, int times) {
		if(times < 1) {
			return matrix;
		}
		
		for(int k = 0; k < times; k++) {
			double[][] base = new double[matrix.length][];
			for(int i = 0; i < matrix.length; i++) base[i] = new double[matrix.length];
			
			for(int i = 0; i < matrix.length; i++) {
				for(int j = 0; j < matrix.length; j++) {
					double val = 0;
					for(int z = 0; z < matrix.length; z++) {
						val += matrix[i][z]*matrix[z][j];
					}
					if(val <= ZERO_THRESHOLD) val = 0;
					base[i][j] = val;
				}
			}
			
			matrix = base;
		}
		
		return matrix;
	}
	
	/**
	 * Takes the provided matrix and elevates a random row to the given coefficient, then normalizes all rows.
	 * @param matrix matrix to be inflated.
	 * @param coefficient coefficient to be used.
	 * @return
	 */
	private double[][] inflateMatrix(double[][] matrix, double coefficient) {
		//Select random row
		int k = (int) (Math.random()*matrix.length);
		//Elevate to the coefficient
		for(int j = 0; j < matrix.length; j++) {
			matrix[k][j] = (double) Math.pow(matrix[k][j], coefficient);
		}
		//Normalize row
		for(int i = 0; i < matrix.length; i++) matrix[i] = this.normalizeRow(matrix[i]);
		return matrix;
	}
	
	/**
	 * Generates a copy of the given matrix.
	 * @param matrix matrix to be copied.
	 * @return
	 */
	private double[][] copyOf(double[][] matrix) {
		double[][] copy = new double[matrix.length][];
		for(int i = 0; i < matrix.length; i++) copy[i] = Arrays.copyOf(matrix[i], matrix[i].length);
		return copy;
	}
	
	public static void main(String[] args) {
		double[][] basic = {
				{},
				{},
				{},
				{},
				{},
				{}
		};
		
		double[][] complex = {
				{0, 0.33f, 0.34f, 0.33f, 0, 0, 0, 0},
				{0.5f, 0, 0.5f, 0, 0, 0, 0, 0},
				{0.399f, 0.3f, 0, 0.3f, 0.001f, 0, 0, 0},
				{0.5f, 0, 0.5f, 0, 0, 0, 0, 0},
				{0, 0, 0.001f, 0, 0, 0.3f, 0.3f, 0.399f},
				{0, 0, 0, 0, 0.5f, 0, 0, 0.5f},
				{0, 0, 0, 0, 0.5f, 0, 0, 0.5f},
				{0, 0, 0, 0, 0.34f, 0.33f, 0.33f, 0},
			};
		
		double[][] test = complex;
		for(double[] arr : test) {
			double sum = 0;
			for(int i = 0; i < arr.length; i++) sum += arr[i];
			System.out.println(String.format("%s || %f", Arrays.toString(arr), sum));
		}
		
		MCLJob job = new MCLJob(test);
		job.run();
		List<List<Integer>>results = job.getResults();
		for(List<Integer> arr : results) {
			System.out.println(String.format("%s", Arrays.toString(arr.toArray())));
		}
	}

	public List<List<Integer>> getResults() {
		return clusters;
	}
}