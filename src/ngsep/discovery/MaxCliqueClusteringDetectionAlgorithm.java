package ngsep.discovery;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import JSci.maths.statistics.SampleStatistics;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.MaximalCliquesFinder;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledInversion;
import ngsep.variants.CalledLargeIndel;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;

public class MaxCliqueClusteringDetectionAlgorithm implements LongReadVariantDetectorClusteringAlgorithm {
	
	public static final double DEFAULT_PD_NORM_FACTOR = 900;
	public static final double DEFAULT_EDGE_TRESHOLD = 0.7;
	public static final double MAX_DOWNSTREAM_CONSENSUS_LENGTH = 50;
	
	private double pdNormFactor = DEFAULT_PD_NORM_FACTOR;
	private double edgeTreshold = DEFAULT_EDGE_TRESHOLD;
	private GenomicRegionSortedCollection<GenomicVariant> signatures;

	public MaxCliqueClusteringDetectionAlgorithm(ReferenceGenome genome, 
			GenomicRegionSortedCollection<GenomicVariant>  signatures, int SVEventLength){
		this.signatures = signatures;
	}
	
	public Map<String, List<List<Integer>>> callVariantClusters() {
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
				for(int i = 0; i < typePartition.size() - 1; i++) {
					GenomicVariant current = typePartition.get(i);
					GenomicVariant next = typePartition.get(i+1);
					toCluster.add(current);
					if(!testDownstreamSignatureCompatibility(current,next) || toCluster.size() >= 300 ||
							i==typePartition.size() - 2) {
						List<Integer> idxs = new ArrayList<>();
						for(GenomicVariant sign:toCluster) idxs.add(signList.indexOf(sign));
						if(toCluster.size() < 4) {
							toCluster.clear();
							continue;
						}
						else{
							boolean[][] adjMatrix = calculateAdjacencyMatrix(toCluster, pdNormFactor, edgeTreshold);
							//System.out.println("#partition size: " + toCluster.size());
							//System.out.println("#Processing partition with types: " + getSignatureType(typePartition.get(0)) 
								//	+ " from: " + idxs.get(0) + " to: " + idxs.get(idxs.size()-1));
							ArrayList<ArrayList<Integer>> maxCliques = MaximalCliquesFinder
									.callMaxCliqueFinder(idxs, adjMatrix);
							//System.out.println("Cliques found");
							chrClusters.addAll(maxCliques);
							toCluster.clear();
						}
					}
				}
			}
			clusters.put(k, chrClusters);
		}
		return clusters;
	}
	
	public static double calculateSPD(GenomicVariant sign1, GenomicVariant sign2, double PDNormFactor) {
		double SPD = 0;
		double SD = 0;
		int PD = 0;
		int span1 = sign1.length();
		int span2 = sign2.length();
		int first1 = sign1.getFirst();
		int last1 = sign1.getLast();
		int first2 = sign2.getFirst();
		int last2 = sign2.getLast();
		if(last1 - first1 < 2) last1 = first1 + span1 - 1;
		if(last2 - first2 < 2) last2 = first2 + span2 - 1;
		SD = calculateSD(span1, span2);
		PD = calculatePD(first1, first2, last1, last2);
		SPD = SD + (double) PD/PDNormFactor;
		return SPD;
	}
	
	public static double calculateSD(int span1, int span2) {
		double SD = Math.abs((span1)-(span2));
		SD = (double) SD/Math.max((span1), (span2));
		return SD;
	}
	
	public static int calculatePD(int first1, int first2, int last1, int last2) {
		int PD = Math.min(Math.abs(first1-first2), Math.abs(last1-last2));
		PD = Math.min(PD, Math.abs(((first1-last1)/2)) - ((first2-last2)/2));
		return PD;
	}
	
	public static boolean[][] calculateAdjacencyMatrix(List<GenomicVariant> candidateSignatures, double PDNormFactor, 
			double edgeTreshold) {
		int n = candidateSignatures.size();
		boolean[][] adjacencyMatrix = new boolean[n][n];
		for(int i = 0; i < n; i++) {
			GenomicVariant si = candidateSignatures.get(i);
			for(int j = 0; j < n; j++) {
				GenomicVariant sj = candidateSignatures.get(j);
				double spd = calculateSPD(si, sj, PDNormFactor);
				if(spd < edgeTreshold && !si.equals(sj)) {
					adjacencyMatrix[i][j] = true;
				}
			}
		}
		return adjacencyMatrix;
	}	
	
	public boolean testDownstreamSignatureCompatibility(GenomicVariant s1, GenomicVariant s2) {
		return s2.getFirst() - s1.getLast() < MAX_DOWNSTREAM_CONSENSUS_LENGTH;
	}
	/**
	public byte getSignatureType(GenomicRegion signature) {
		byte type = 0;
		if(signature instanceof CalledLargeIndel) {
			CalledLargeIndel indelSign = (CalledLargeIndel) signature;
			type = indelSign.getType();
		}
		if(signature instanceof CalledInversion) {
			type = GenomicVariant.TYPE_INVERSION;
		}
		return type;
	}**/
	
	/**
	@Override
	public void setRefGenome(ReferenceGenome genome) {
		this.refGenome = genome;
	}
	public void setSignatures(GenomicRegionSortedCollection<GenomicRegion> signatures) {
		this.signatures = signatures;
	}
	**/
}
