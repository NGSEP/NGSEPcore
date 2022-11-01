package ngsep.discovery;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;

import ngsep.clustering.DBSCANClusteringAlgorithm;
import ngsep.clustering.Pair;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.MaximalCliquesFinder;
import ngsep.math.Distribution;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;

public class DBSCANClusteringDetectionAlgorithm implements LongReadVariantDetectorClusteringAlgorithm {
	
	public final static int MIN_DEFAULT_POINTS = 4;
	public final static double DEFAULT_EPSILON = 250;
	public static final double MAX_DOWNSTREAM_CONSENSUS_LENGTH = 4000;
	
	private int minPoints = MIN_DEFAULT_POINTS;
	private double epsilon = DEFAULT_EPSILON;
	
	private GenomicRegionSortedCollection<GenomicVariant> signatures;
	
	public DBSCANClusteringDetectionAlgorithm(ReferenceGenome genome, 
			GenomicRegionSortedCollection<GenomicVariant>  signatures, int SVEventLength){
		this.signatures = signatures;
	}

	@Override
	public Map<String, List<List<Integer>>> callVariantClusters() {
		// TODO Auto-generated method stub
		Map<String, List<List<Integer>>> clusters = new LinkedHashMap<>();
		List<String> keys = signatures.getSequenceNames().getNamesStringList();
		for(String k:keys){
			List<GenomicVariant> signList = signatures.getSequenceRegions(k).asList();
			List<List<Integer>> chrClusters = new ArrayList<>();
			List<GenomicVariant> indel = new ArrayList<>();
			List<GenomicVariant> del = new ArrayList<>();
			List<GenomicVariant> ins = new ArrayList<>();
			List<GenomicVariant> inv = new ArrayList<>();
			List<GenomicVariant> dup = new ArrayList<>();
			List<GenomicVariant> tan = new ArrayList<>();
			List<GenomicVariant> und = new ArrayList<>();
			List<List<GenomicVariant>> signsPerType = new ArrayList<>();
			for(GenomicVariant sign:signList) {
				byte type = sign.getType();

				if(GenomicVariant.TYPE_INDEL == type) indel.add(sign);
				else if(GenomicVariant.TYPE_LARGEDEL == type) del.add(sign);
				else if(GenomicVariant.TYPE_LARGEINS == type) ins.add(sign);
				else if(GenomicVariant.TYPE_DUPLICATION == type) dup.add(sign);
				else if(GenomicVariant.TYPE_INVERSION == type) inv.add(sign);
				else if(GenomicVariant.TYPE_CNV == type) tan.add(sign);
				else if(GenomicVariant.TYPE_UNDETERMINED == type) und.add(sign);
			}
			signsPerType.add(indel);
			signsPerType.add(del);
			signsPerType.add(ins);
			signsPerType.add(dup);
			signsPerType.add(inv);
			signsPerType.add(tan);
			signsPerType.add(und);
			for(List<GenomicVariant> typePartition:signsPerType) {
				if(typePartition.isEmpty()) continue;
				boolean[] visited = new boolean[typePartition.size() - 1];
				List<GenomicVariant> toCluster = new ArrayList<>();
				for(int i = 0; i < typePartition.size() - 1; i++) {
					GenomicVariant current = typePartition.get(i);
					visited[i] = true;
					GenomicVariant next = typePartition.get(i+1);
					toCluster.add(current);
					if(!testDownstreamSignatureCompatibility(current,next) || i==typePartition.size() - 2 || 
							toCluster.size() > 10000) {
						List<Integer> idxs = new ArrayList<>();
						for(GenomicVariant sign:toCluster) {
							idxs.add(signList.indexOf(sign));
						}
							//double [][] distanceMatrix = calculateDistanceMatrix(toCluster);
							List<List<Integer>> adjacencyList = constructSignatureGraph(toCluster);
							DBSCANClusteringAlgorithm instance = new DBSCANClusteringAlgorithm();
							//List<List<Integer>> partClusters = instance
							//		.runDBSCANClustering(idxs, distanceMatrix, minPoints);
							List<List<Integer>> partClusters = instance
									.runDBSCANClustering(idxs, adjacencyList, minPoints);
							chrClusters.addAll(partClusters);
							List<Integer> noisePoints = instance.getNoisePoints();
							if(noisePoints.size() > 1) {
								List<GenomicVariant> noiseSignsToCluster = new ArrayList<>();
								for(int np:noisePoints) noiseSignsToCluster.add(signList.get(np));
								//double [][] weakDistanceMatrix = calculateDistanceMatrix(noiseSignsToCluster);
								List<List<Integer>> weakAdjacencyList = constructSignatureGraph(noiseSignsToCluster);
								DBSCANClusteringAlgorithm inst = new DBSCANClusteringAlgorithm();
								List<List<Integer>> weakClusters = inst
										.runDBSCANClustering(noisePoints, weakAdjacencyList, 0);
								chrClusters.addAll(weakClusters);
							}
							toCluster.clear();
					}
				}
			}
			clusters.put(k, chrClusters);
		}
		
		return clusters;
	}

	private List<List<Integer>> constructSignatureGraph(List<GenomicVariant> candidateSignatures) {
		List<List<Integer>> graph = new ArrayList<>();
		int n = candidateSignatures.size();
		for(int i = 0; i < n; i++) graph.add(new ArrayList<>());
		for(int i = 0; i < n; i++){
			GenomicVariant iSign = candidateSignatures.get(i);
			List<Integer> iNeighbors = graph.get(i);
			for(int j = i; j < n; j++){
				GenomicVariant jSign = candidateSignatures.get(j);
				List<Integer> jNeighbors = graph.get(j);
				if(i!=j){
					double distance = calculateThreeDimEuclideanDistance(iSign, jSign);
					if(distance < epsilon) {
						iNeighbors.add(j);
						jNeighbors.add(i);
					}
				}
			}
		}
		return graph;
	}

	public static double calculateThreeDimEuclideanDistance(GenomicVariant s1, GenomicVariant s2) {
		return calculateThreeDimEuclideanDistance(s1.getFirst(), s2.getFirst(), s1.getLast(), s2.getLast(), 
				s1.length(), s2.length());
	}
	
	public static double calculateThreeDimEuclideanDistance(int x1, int x2, int y1, int y2, int z1, int z2) {
		double eDistance;
		double xDistance = Math.pow(x2-x1 ,2);
		double yDistance = Math.pow(y2-y1, 2);
		double zDistance = Math.pow(z2-z1, 2);
		eDistance = xDistance + yDistance + zDistance;
		eDistance = Math.sqrt(eDistance);
		return eDistance;
	}
	
	public static boolean testDownstreamSignatureCompatibility(GenomicVariant s1, GenomicVariant s2) {
		return Math.abs(s2.getFirst() - s1.getFirst()) < MAX_DOWNSTREAM_CONSENSUS_LENGTH;
	}
	private double [][] calculateDistanceMatrix(List<GenomicVariant> candidateSignatures) {
		// TODO Auto-generated method stub
		int n = candidateSignatures.size();
		double [][] distanceMatrix = new double[n][n];
		for(int i = 0; i < n; i++) {
			GenomicVariant si = candidateSignatures.get(i);
			for(int j = 0; j < n; j++) {
				GenomicVariant sj = candidateSignatures.get(j);
				if(!si.equals(sj)) {
					double distance = calculateThreeDimEuclideanDistance(si, sj);
					distanceMatrix[i][j] = distance;
				}
			}
		}
		return distanceMatrix;
	}
	private List<List<Integer>> clusterConnectedComponents(List<Integer> noisePoints, double[][] weakDistanceMatrix,
														   double epsilon) {
		// TODO Auto-generated method stub
		DBSCANClusteringAlgorithm inst = new DBSCANClusteringAlgorithm();
		List<List<Integer>> weakClusters = inst
				.runDBSCANClustering(noisePoints, weakDistanceMatrix, 0, epsilon);
		return weakClusters;
	}
}
