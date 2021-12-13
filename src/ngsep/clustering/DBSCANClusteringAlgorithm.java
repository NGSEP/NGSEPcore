package ngsep.clustering;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import ngsep.discovery.DBSCANClusteringDetectionAlgorithm;

public class DBSCANClusteringAlgorithm {
	
	public static final int UNDEFINED_LABEL = -1;
	public static final int NOISE_LABEL = 0;
	
	public static List<List<Integer>> runDBSCANClustering(List<Integer> idxs,
			double [][] distanceMatrix, int minPts, double epsilon) {
		// TODO Auto-generated method stub
		int [] clusterLabels = DBSCAN(idxs, distanceMatrix, minPts, epsilon);
		List<List<Integer>> clusters = assignClustersByLabel(idxs, clusterLabels);
		return clusters;
	}
	
	private static int [] DBSCAN(List<Integer> idxs, double [][] distanceMatrix, int minPts, double epsilon) {
		int [] labels = new int[idxs.size()];
		Arrays.fill(labels, UNDEFINED_LABEL);
		int c = 0;
		for(int i = 0; i < idxs.size(); i++) {
			if(labels[i] != UNDEFINED_LABEL) continue;
			//This list must return the idxs.indexOf(neighbor) as the values, which correspond to i and labels' indexes
			List<Integer> neighbors = rangeQuery(distanceMatrix, epsilon, i);
			if (neighbors.size() < minPts) {
				labels[i] = NOISE_LABEL;
				continue;
			}
			c++;
			labels[i] = c;
			Queue<Integer> seedQueue = new LinkedList<>(neighbors);
			while(!seedQueue.isEmpty()) {
				int j = seedQueue.remove();
				//if(labels[j] == NOISE_LABEL) labels[j] = c;
				if(labels[j] != UNDEFINED_LABEL) continue;
				neighbors = rangeQuery(distanceMatrix, epsilon, j);
				labels[j] = c;
				if(neighbors.size() < minPts) continue;
				seedQueue.addAll(neighbors);
			}
		}
		return labels;
	}

	private static List<List<Integer>> assignClustersByLabel(List<Integer> idxs, int[] labels) {
		// TODO Auto-generated method stub
		Map <Integer, List<Integer>> clustersMap = new HashMap<>();
		for(int i = 0; i < labels.length; i++) {
			int c = labels[i];
			if(c==NOISE_LABEL) continue;
			int idx = idxs.get(i);
			List<Integer> cluster = clustersMap.computeIfAbsent(c, v -> new ArrayList<>());
			cluster.add(idx);
		}
		List<List<Integer>> clusters = new ArrayList<>(clustersMap.values());
		return clusters;
	}

	private static List<Integer> rangeQuery(double[][] distanceMatrix, double epsilon,
			int j) {
		// TODO Auto-generated method stub
		List<Integer> neighbors = new ArrayList<>();
		for(int i = 0; i < distanceMatrix.length; i++) {
			double distance = distanceMatrix[i][j];
			if(distance <= epsilon && i != j) {
				neighbors.add(i);
			}
		}
		return neighbors;
	}
	
	//Main method used for testing the implementation outside the ngsep context. Not supposed to be used in application.
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		//Test a dataset in contrast with pythons sklearn DBSCAN implementation
		List<Double> xDimNums = new ArrayList<>();
		List<Double> yDimNums = new ArrayList<>();
		try(BufferedReader r = new BufferedReader(new FileReader(args[0]))){
			String line = r.readLine();
			while(line != null) {
				String[] nums = line.split(",");
				xDimNums.add(Double.parseDouble(nums[0]));
				yDimNums.add(Double.parseDouble(nums[1]));
				line = r.readLine();
			}
			if(xDimNums.size() != yDimNums.size()) throw new Exception("x and y dimensions have different sizes");
		}
		int n = xDimNums.size();
		double [][] distanceMatrix = new double [n][n];
		for(int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double x1 = xDimNums.get(i);
				double x2 = xDimNums.get(j);
				double y1 = yDimNums.get(i);
				double y2 = yDimNums.get(j);
				double xDistance = Math.pow((x2-x1), 2.0);
				double yDistance = Math.pow((y2-y1), 2.0);
				double eDistance = xDistance + yDistance;
				eDistance = Math.sqrt((double) eDistance);
				distanceMatrix[i][j] = eDistance;
			}
		}
		List<Integer> idxs = new ArrayList<>();
		for (int idx = 0; idx < n; idx++) {
			idxs.add(idx);
			
		}
		/**
		int[] testClusters = DBScan(idxs, distanceMatrix, 13, 1.0);
		if(xDimNums.size() != testClusters.length) throw new Exception("x and labels dimensions have different sizes");
		try(PrintWriter writer = new PrintWriter("ngsep_labels.txt")){
			for(int label : testClusters) writer.write(label + "\n");
		}**/
		List<List<Integer>> clusters = runDBSCANClustering(idxs, distanceMatrix, 13, 1.0);
		try(PrintWriter writer = new PrintWriter("ngsep_clusters.txt")){
			int c = 0;
			for(List<Integer> cluster : clusters) {
				for(int idx : cluster) {
					writer.write(xDimNums.get(idx) + "," + yDimNums.get(idx) + "," + c + "\n");
				}
				c++;
			}
		}
	}
}
