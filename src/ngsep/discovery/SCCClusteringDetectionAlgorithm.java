package ngsep.discovery;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.StronglyConnectedComponents;
import ngsep.variants.GenomicVariant;

public class SCCClusteringDetectionAlgorithm implements LongReadVariantDetectorClusteringAlgorithm {
	
	public final static double DEFAULT_MIN_ED = 100;
	public static final double MAX_DOWNSTREAM_CONSENSUS_LENGTH = 4000;

	private GenomicRegionSortedCollection<GenomicVariant> signatures;
	
	private double minED = DEFAULT_MIN_ED;
	/**
	private int testGraphNumber = 0;
	private int maxTests = 10;
	**/
	public SCCClusteringDetectionAlgorithm(ReferenceGenome genome, 
			GenomicRegionSortedCollection<GenomicVariant>  signatures, int SVEventLength){
		this.signatures = signatures;
	}
	
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
					if(!DBSCANClusteringDetectionAlgorithm.
							testDownstreamSignatureCompatibility(current,next) || i==typePartition.size() - 2 || 
							toCluster.size() > 10000) {
						List<Integer> idxs = new ArrayList<>();
						for(GenomicVariant sign:toCluster) {
							idxs.add(signList.indexOf(sign));
						}
							List<List<Integer>> directedGraph = createGraphAsAdjacencyList(toCluster);
							StronglyConnectedComponents instance = new StronglyConnectedComponents();
							List<List<Integer>> partClusters = instance.
									computeStronglyConnectedComponents(idxs, directedGraph);
							chrClusters.addAll(partClusters);
							/**if(r==1 && maxTests>0 && toCluster.size()>1) {
								int r = new Random().nextInt(100);
								try {
									printDirectedGraphForTesting(directedGraph, testGraphNumber);
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
								testGraphNumber++;
								maxTests--;
							}
							**/
							/**DBSCANClusteringAlgorithm instance = new DBSCANClusteringAlgorithm();
							List<List<Integer>> partClusters = instance
									.runDBSCANClustering(idxs, distanceMatrix, minPoints, epsilon);
							chrClusters.addAll(partClusters);
							List<Integer> noisePoints = instance.getNoisePoints();
							if(noisePoints.size() > 1) {
								List<GenomicVariant> noiseSignsToCluster = new ArrayList<>();
								for(int np:noisePoints) noiseSignsToCluster.add(signList.get(np));
								double [][] weakDistanceMatrix = calculateDistanceMatrix(noiseSignsToCluster);
								List<List<Integer>> weakClusters = 
										clusterConnectedComponents(noisePoints, weakDistanceMatrix, epsilon);
								chrClusters.addAll(weakClusters);
							}**/
							toCluster.clear();
					}
				}
			}
			clusters.put(k, chrClusters);
		}
		
		return clusters;
	}
	/**
	private void printDirectedGraphForTesting(List<List<Integer>> directedGraph, int testGraphNumber) throws IOException {
		// TODO Auto-generated method stub
		try(PrintWriter w = new PrintWriter("graphTesting" + testGraphNumber + ".tsv")){
			w.println("source" + "\t" + "target");
			for(int i = 0; i < directedGraph.size(); i++) {
				List<Integer> iNeighbors = directedGraph.get(i);
				for(int j : iNeighbors) {
					w.println(i + "\t" + j);
				}
			}
		}
	}**/

	private List<List<Integer>> createGraphAsAdjacencyList(List<GenomicVariant> signatures) {
		// TODO Auto-generated method stub
		List<List<Integer>> adjacencyList = new ArrayList<>();
		for(int l = 0; l < signatures.size(); l++) adjacencyList.add(new ArrayList<>());
		for(int i = 0; i < signatures.size(); i++) {
			GenomicVariant si = signatures.get(i);
			for(int j = 0; j < signatures.size(); j++) {
				GenomicVariant sj = signatures.get(j);
				if(i==j) continue;
				if(i<j) {
					double distance = DBSCANClusteringDetectionAlgorithm.calculateThreeDimEuclideanDistance(si, sj);
					if(distance < minED) adjacencyList.get(i).add(j);
				}
				else {
					if(sj.length() == si.length()) adjacencyList.get(i).add(j);
				}
			}
		}
		return adjacencyList;
	}
}
