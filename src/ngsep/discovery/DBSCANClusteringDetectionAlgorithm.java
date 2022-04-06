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
	
	public final static int MIN_DEFAULT_POINTS = 3;
	public final static double DEFAULT_EPSILON = 250;
	public static final double MAX_DOWNSTREAM_CONSENSUS_LENGTH = 4000;
	
	private int minPoints = MIN_DEFAULT_POINTS;
	private double epsilon = DEFAULT_EPSILON;
	//just for test
	//private boolean ranDistTest = false;
	
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
				//System.out.println(" Sign in type list: " + sign.getFirst() + " last: " + sign.getLast() + " chr: "
					//	+ sign.getSequenceName() + " with length: " + sign.length());
				//System.out.println(GenomicVariantImpl.getVariantTypeName(type));
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
				//toCluster.addAll(typePartition);
				for(int i = 0; i < typePartition.size() - 1; i++) {
					//to remove after testing
					//int r = new Random().nextInt(typePartition.size()-1);
					GenomicVariant current = typePartition.get(i);
					visited[i] = true;
					GenomicVariant next = typePartition.get(i+1);
					toCluster.add(current);
					//if(i==typePartition.size() - 2) System.out.println(" Last sign in partition: " + next.getFirst() + " last: " + next.getLast() + " chr: "
						//	+ next.getSequenceName() + " with length: " + next.length());
					if(!testDownstreamSignatureCompatibility(current,next) || i==typePartition.size() - 2) {
						List<Integer> idxs = new ArrayList<>();
						//System.out.println("DSC: " + testDownstreamSignatureCompatibility(current,next));
						for(GenomicVariant sign:toCluster) {
							idxs.add(signList.indexOf(sign));
						}
							//only for testing
							/**
							if(r == i && !ranDistTest) {
								try {
									testDistributions(GenomicVariantImpl.getVariantTypeName(signList.get(0).getType()), 
											toCluster);
									ranDistTest = true;
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
							}**/
							double [][] distanceMatrix = calculateDistanceMatrix(toCluster);
							//double [][] distanceMatrix = calculateNormalizedDistanceMatrix(toCluster);
							//System.out.println("#partition size: " + toCluster.size());
							//for(GenomicVariant sign:toCluster) {
								//System.out.println(" Sign in partition list: " + sign.getFirst() + " last: " + sign.getLast() + " chr: "
									//	+ sign.getSequenceName() + " with length: " + sign.length());
							//}
							//System.out.println("#Processing partition with types: " + getSignatureType(typePartition.get(0)) 
								//	+ " from: " + idxs.get(0) + " to: " + idxs.get(idxs.size()-1));
							DBSCANClusteringAlgorithm instance = new DBSCANClusteringAlgorithm();
							List<List<Integer>> partClusters = instance
									.runDBSCANClustering(idxs, distanceMatrix, minPoints, epsilon);
							//System.out.println("Cliques found");
							//for(List<Integer> cluster : partClusters) if(cluster.size() > 3) chrClusters.add(cluster);
							List<Integer> noisePoints = instance.getNoisePoints();
							chrClusters.addAll(partClusters);
							if(noisePoints.size() > 2) {
								List<GenomicVariant> noiseSigns = new ArrayList<>();
								for(int np:noisePoints) noiseSigns.add(signList.get(np));
								boolean[][] adjMatrix = MaxCliqueClusteringDetectionAlgorithm.
										calculateAdjacencyMatrix(noiseSigns, 
												MaxCliqueClusteringDetectionAlgorithm.DEFAULT_PD_NORM_FACTOR, 
												0.8);
								ArrayList<ArrayList<Integer>> weakClusters = MaximalCliquesFinder
										.callMaxCliqueFinder(noisePoints, adjMatrix);
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
	/**
	private boolean[][] calculateAdjacencyMatrix(List<GenomicVariant> noiseSigns) {
		// TODO Auto-generated method stub
		int n = noiseSigns.size();
		boolean [][] adjMatrix = new boolean[n][n];
		for(int i = 0; i < n; i++) {
			GenomicVariant si = noiseSigns.get(i);
			for(int j = 0; j < n; j++) {
				GenomicVariant sj = noiseSigns.get(j);
				if(!si.equals(sj)) {
					double distance = calculateThreeDimEuclideanDistance(si, sj);
					adjMatrix[i][j] = distance < epsilon;
				}
			}
		}
		return adjMatrix;
	}
	**/
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
						/**System.out.println("distance between " + si.getFirst() + " " + si.getLast() + " " 
								+ si.length() + " and " + sj.getFirst() + " " + sj.getLast() + " " +
								sj.length() + " is " + distance);**/
					distanceMatrix[i][j] = distance;
				}
			}
		}
		return distanceMatrix;
	}

	private double calculateThreeDimEuclideanDistance(GenomicVariant s1, GenomicVariant s2) {
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
	
	private boolean testDownstreamSignatureCompatibility(GenomicVariant s1, GenomicVariant s2) {
		//System.out.println(s1.getLast() + " " + s2.getFirst() + " " + Math.abs(s2.getFirst() - s1.getFirst()));
		return Math.abs(s2.getFirst() - s1.getFirst()) < MAX_DOWNSTREAM_CONSENSUS_LENGTH;
	}

	/**
	//Only for testing
	public void testDistributions(String listType, List<GenomicVariant> variants) throws IOException {
		try(PrintWriter writer = new PrintWriter("distributions.txt")){
			List<Integer> firstDistValues = new ArrayList<>();
			List<Integer> lastDistValues = new ArrayList<>();
			List<Integer> lengthDistValues = new ArrayList<>();
			SortedMap<Integer, Integer> firstDist = new TreeMap<>();
			SortedMap<Integer, Integer> lastDist = new TreeMap<>();
			SortedMap<Integer, Integer> lengthDist = new TreeMap<>();
			int initJ = 1;
			int j = 1;
			for(int i = 0; i < variants.size(); i++){
				GenomicVariant var = variants.get(i);
				for (; j < variants.size(); j++) {
					GenomicVariant var2 = variants.get(j);
					firstDistValues.add(Math.abs(var2.getFirst() - var.getFirst()));
					lastDistValues.add(Math.abs(var.getLast() - var2.getLast()));
					lengthDistValues.add(Math.abs(var.length() - var2.length()));
				}
				initJ++;
				j = initJ;
			}
			for(int x = 0; x < firstDistValues.size(); x++) {
				firstDist.compute(firstDistValues.get(x), (k,v) -> (v == null) ? 1:v+1);
				lastDist.compute(lastDistValues.get(x), (k,v) -> (v == null) ? 1:v+1);
				lengthDist.compute(lengthDistValues.get(x), (k,v) -> (v == null) ? 1:v+1);
			}
			for(int i = 0; i < firstDist.lastKey(); i++) {
				if(!firstDist.containsKey(i)) firstDist.put(i, 0); 
			}
			for(int i = 0; i < lastDist.lastKey(); i++) {
				if(!lastDist.containsKey(i)) lastDist.put(i, 0); 
			}
			for(int i = 0; i < lengthDist.lastKey(); i++) {
				if(!lengthDist.containsKey(i)) lengthDist.put(i, 0); 
			}
			writer.println("Distributions for random types: " + listType);
			writer.println("------------------------------------------------------------------------------------------ ");
			writer.println("FIRST_DISTRIBUTION");
			 for(Map.Entry<Integer, Integer> entry : firstDist.entrySet()) {
				 writer.println(entry.getKey() + "\t" + entry.getValue());
			 	}
			 writer.println("LAST_DISTRIBUTION");
			 for(Map.Entry<Integer, Integer> entry : lastDist.entrySet()) {
				 writer.println(entry.getKey() + "\t" + entry.getValue());
				 }
			 writer.println("LENGTH_DISTRIBUTION");
			 for(Map.Entry<Integer, Integer> entry : lengthDist.entrySet()) {
				 writer.println(entry.getKey() + "\t" + entry.getValue());
				 }
		}
	}

	private double calculateNormalizedThreeDimEuclideanDistance(GenomicVariant s1, GenomicVariant s2, 
			List<Integer> firstList, List<Integer> lastList, List<Integer> lengthList) {
		int firstDistance = Math.abs(s2.getFirst() - s1.getFirst());
		int lastDistance = Math.abs(s2.getLast() - s1.getLast());
		int lengthDistance = Math.abs(s2.length() - s1.length());
		int minFirst = firstList.get(0);
		int minLast = lastList.get(0);
		int minLength = lengthList.get(0);
		int n = firstList.size();
		int maxFirst = firstList.get(n-1);
		int maxLast = lastList.get(n-1);
		int maxLength = lengthList.get(n-1);
		return calculateNormalizedThreeDimEuclideanDistance(firstDistance, lastDistance, lengthDistance,
				minFirst, maxFirst, minLast, maxLast, minLength, maxLength);
	}

	public double calculateNormalizedThreeDimEuclideanDistance(int x, int y, int z, int minFirst, int maxFirst,
			int minLast, int maxLast,int minLength, int maxLength) {
		double firstNorm = minMaxScaler(x, minFirst, maxFirst);
		double lastNorm = minMaxScaler(y, minLast, maxLast);
		double lengthNorm = minMaxScaler(z, minLength, maxLength);
		double xDistance = Math.pow(firstNorm ,2);
		double yDistance = Math.pow(lastNorm, 2);
		double zDistance = Math.pow(lengthNorm, 2);
		double eDistance = xDistance + yDistance + zDistance;
		eDistance = Math.sqrt(eDistance);
		return eDistance;
	}
	private double [][] calculateNormalizedDistanceMatrix(List<GenomicVariant> candidateSignatures) {
		// TODO Auto-generated method stub
		List<Integer> firstDistValues = new ArrayList<>();
		List<Integer> lastDistValues = new ArrayList<>();
		List<Integer> lengthDistValues = new ArrayList<>();
		int init = 1;
		int y = 1;
		for(int x = 0; x < candidateSignatures.size(); x++){
			GenomicVariant var1 = candidateSignatures.get(x);
			for (; y < candidateSignatures.size(); y++) {
				GenomicVariant var2 = candidateSignatures.get(y);
				firstDistValues.add(Math.abs(var2.getFirst() - var1.getFirst()));
				lastDistValues.add(Math.abs(var2.getLast() - var1.getLast()));
				lengthDistValues.add(Math.abs(var2.length() - var1.length()));
			}
			init++;
			y = init;
		}
		Collections.sort(firstDistValues);
		Collections.sort(lastDistValues);
		Collections.sort(lengthDistValues);
		//Distribution firstDist = LongReadStructuralVariantDetector.calculateDistribution(firstDistValues);
		//Distribution lastDist = LongReadStructuralVariantDetector.calculateDistribution(lastDistValues);
		//Distribution lengthDist = LongReadStructuralVariantDetector.calculateDistribution(lengthDistValues);
		int n = candidateSignatures.size();
		double [][] distanceMatrix = new double[n][n];
		for(int i = 0; i < n; i++) {
			GenomicVariant si = candidateSignatures.get(i);
			for(int j = 0; j < n; j++) {
				GenomicVariant sj = candidateSignatures.get(j);
				if(!si.equals(sj)) {
					double distance = calculateNormalizedThreeDimEuclideanDistance(si, sj,
							firstDistValues, lastDistValues, lengthDistValues);
					if(si.length()>=10000) {
						System.out.println("distance between " + si.getFirst() + " " + si.getLast() + " " 
								+ si.length() + " and " + sj.getFirst() + " " + sj.getLast() + " " +
								sj.length() + " is " + distance);
					}
					distanceMatrix[i][j] = distance;
				}
			}
		}
		return distanceMatrix;
	}
	public static void main (String[] args) {
		double dist = calculateThreeDimEuclideanDistance(28216662, 28219678, 28220049, 28219724, 3388, 45);
		System.out.println(dist);
	}
**/
	

}
