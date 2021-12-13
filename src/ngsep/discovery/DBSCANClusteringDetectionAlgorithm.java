package ngsep.discovery;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import ngsep.clustering.DBSCANClusteringAlgorithm;
import ngsep.clustering.Pair;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.MaximalCliquesFinder;
import ngsep.variants.GenomicVariant;

public class DBSCANClusteringDetectionAlgorithm implements LongReadVariantDetectorClusteringAlgorithm {
	
	public final static int MIN_DEFAULT_POINTS = 4;
	public final static double DEFAULT_EPSILON = 50;
	public static final double MAX_DOWNSTREAM_CONSENSUS_LENGTH = 1000;
	
	public int minPoints = MIN_DEFAULT_POINTS;
	public  double epsilon = DEFAULT_EPSILON;
	
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
				//System.out.println(" Sign first: " + sign.getFirst() + " last: " + sign.getLast() + " chr: "
					//	+ sign.getSequenceName() + " with length: " + sign.length());
				//System.out.println(type);
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
				List<GenomicVariant> toCluster = new ArrayList<>();
				//toCluster.addAll(typePartition);
				for(int i = 0; i < typePartition.size() - 1; i++) {
					GenomicVariant current = typePartition.get(i);
					GenomicVariant next = typePartition.get(i+1);
					toCluster.add(current);
					if(!testDownstreamSignatureCompatibility(current,next)) {
						List<Integer> idxs = new ArrayList<>();
						for(GenomicVariant sign:toCluster) idxs.add(signList.indexOf(sign));
						if(toCluster.size() < 4) {
							toCluster.clear();
							continue;
						}
						else{
							double [][] distanceMatrix = calculateDistanceMatrix(toCluster);
							//System.out.println("#partition size: " + toCluster.size());
							//System.out.println("#Processing partition with types: " + getSignatureType(typePartition.get(0)) 
								//	+ " from: " + idxs.get(0) + " to: " + idxs.get(idxs.size()-1));
							List<List<Integer>> partClusters = DBSCANClusteringAlgorithm
									.runDBSCANClustering(idxs, distanceMatrix, minPoints, epsilon);
							//System.out.println("Cliques found");
							for(List<Integer> cluster : partClusters) if(cluster.size() > 3) chrClusters.add(cluster);;
							toCluster.clear();
						}
					}
				}
			}
			clusters.put(k, chrClusters);
		}
		return clusters;
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
	
	public boolean testDownstreamSignatureCompatibility(GenomicVariant s1, GenomicVariant s2) {
		return s2.getFirst() - s1.getLast() < MAX_DOWNSTREAM_CONSENSUS_LENGTH;
	}
	/**
	public static void main (String[] args) {
		double dist = calculateThreeDimEuclideanDistance(28216662, 28219678, 28220049, 28219724, 3388, 45);
		System.out.println(dist);
	}
	**/
}
