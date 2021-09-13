package ngsep.discovery;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import JSci.maths.statistics.SampleStatistics;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.MaximalCliquesFinder;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;

public class MaxCliqueClusteringDetectionAlgorithm implements LongReadVariantDetectorAlgorithm{
	
	public static final double DEFAULT_PD_NORM_FACTOR = 900;
	public static final double DEFAULT_EDGE_TRESHOLD = 0.7;
	public static final double MAX_CONSENSUS_LENGTH = 100;
	public static final double MAX_DOWNSTREAM_CONSENSUS_LENGTH = 100;
	
	private double pdNormFactor = DEFAULT_PD_NORM_FACTOR;
	private double edgeTreshold = DEFAULT_EDGE_TRESHOLD;
	private GenomicRegionSortedCollection<GenomicVariant> signatures;
	
	@Override
	//public void setSignatures(Map<String, List<GenomicVariant>> signatures) {
		// TODO Auto-generated method stub
		//this.signatures = signatures;
	//}
	
	public void setSignatures(GenomicRegionSortedCollection<GenomicVariant> signatures) {
		// TODO Auto-generated method stub
		this.signatures = signatures;
	}
	public double calculateSPD(GenomicVariant sign1, GenomicVariant sign2) {
		double SPD = 0;
		double SD = 0;
		int PD = 0;
		int first1 = sign1.getFirst();
		int last1 = sign1.getLast();
		int first2 = sign2.getFirst();
		int last2 = sign2.getLast();
		SD = Math.abs((last1-first1)-(last2-first2));
		SD = (double) SD/Math.max((last1-first1), (last2-first2));
		PD = Math.min(Math.abs(first1-first2), Math.abs(last1-last2));
		PD = Math.min(PD, Math.abs(((first1-last1)/2)) - ((first2-last2)/2));
		SPD = SD + (double) PD/pdNormFactor;
		return SPD;
	}
	public boolean[][] calculateAdjacencyMatrix(List<GenomicVariant> candidateSignatures) {
		int n = candidateSignatures.size();
		boolean[][] adjacencyMatrix = new boolean[n][n];
		for(int i = 0; i < n; i++) {
			GenomicVariant si = candidateSignatures.get(i);
			for(int j = 0; j < n; j++) {
				GenomicVariant sj = candidateSignatures.get(j);
				double spd = calculateSPD(si, sj);
				if(spd < edgeTreshold && !si.equals(sj)) {
					adjacencyMatrix[i][j] = true;
				}
			}
		}
		return adjacencyMatrix;
	}	
	public short calculateVariantScore(List<GenomicVariant> candidates) {
		double score = 0;
		SampleStatistics calcSpan = new SampleStatistics();
		SampleStatistics calcPos = new SampleStatistics();
		for(GenomicVariant sign:candidates) {
			calcPos.update(sign.getFirst());
			calcSpan.update(sign.length());
		}
		double numSign = Math.min(80, candidates.size());
		double spanMean = calcSpan.getMean();
		double normPos = (numSign/8)*(1-Math.min(1, 
						(Math.sqrt(calcPos.getVariance()))/spanMean));
		double normSpan = (numSign/8)*(1-Math.min(1, 
						(Math.sqrt(calcSpan.getVariance()))/spanMean));
		score = normSpan + normPos + numSign;
		return (short) score;
	}
	public boolean testCandidateCompatibility(GenomicVariant v1, GenomicVariant v2) {
		if(!v1.getSequenceName().equals(v2.getSequenceName())) return false;
		if(v1.getType() != v2.getType()) return false;
		if(Math.abs(v1.getFirst() - v2.getFirst()) < MAX_CONSENSUS_LENGTH) return true; 
		if(Math.abs(v1.getLast() - v2.getLast()) < MAX_CONSENSUS_LENGTH) return true; 
		//if(Math.abs(v1.length() - v2.length()) < 20) return true; 
		return false;
	}
	public boolean testDownstreamSignatureCompatibility(GenomicVariant s1, GenomicVariant s2) {
		return s2.getFirst() - s1.getLast() < MAX_DOWNSTREAM_CONSENSUS_LENGTH;
	}
	public Map<String, List<List<Integer>>> findVariantClusters() {
		Map<String, List<List<Integer>>> clusters = new LinkedHashMap<>();
		List<String> keys = signatures.getSequenceNames().getNamesStringList();
		for(String k:keys) System.out.println("$ sequence: " + k);
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
				String type = GenomicVariantImpl.getVariantTypeName(sign.getType());
				if(GenomicVariant.TYPENAME_INDEL.equals(type)) indel.add(sign);
				else if(GenomicVariant.TYPENAME_LARGEDEL.equals(type)) del.add(sign);
				else if(GenomicVariant.TYPENAME_LARGEINS.equals(type)) ins.add(sign);
				else if(GenomicVariant.TYPENAME_DUPLICATION.equals(type)) dup.add(sign);
				else if(GenomicVariant.TYPENAME_INVERSION.equals(type)) inv.add(sign);
				else if(GenomicVariant.TYPENAME_CNV.equals(type)) tan.add(sign);
				else if(GenomicVariant.TYPENAME_UNDETERMINED.equals(type)) und.add(sign);
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
				System.out.println("#: " + typePartition.size());
				List<GenomicVariant> toCluster = new ArrayList<>();
				for(int i = 0; i < typePartition.size() - 1; i++) {
					GenomicVariant current = typePartition.get(i);
					GenomicVariant next = typePartition.get(i+1);
					toCluster.add(current);
					if(!testDownstreamSignatureCompatibility(current,next) || toCluster.size() >= 100) {
						List<Integer> idxs = new ArrayList<>();
						for(GenomicVariant sign:toCluster) idxs.add(signList.indexOf(sign));
						if(toCluster.size() < 3) {
							chrClusters.add(idxs);
							toCluster.clear();
						}
						else{
							boolean[][] adjMatrix = calculateAdjacencyMatrix(toCluster);
							System.out.println("#partition size: " + toCluster.size());
							System.out.println("Processing partition with types: " + typePartition.get(0).getType() 
									+ " from: " + idxs.get(0) + " to: " + idxs.get(idxs.size()-1));
							ArrayList<ArrayList<Integer>> maxCliques = MaximalCliquesFinder
									.callMaxCliqueFinder(idxs, adjMatrix);
							System.out.println("Cliques found");
							chrClusters.addAll(maxCliques);
							toCluster.clear();
						}
					}
				}
			}
			clusters.put(k, chrClusters);
		}
		System.out.println(" Finished clustering process");
		return clusters;
	}
	public GenomicVariantImpl processClusterToVariant(List<GenomicVariant> clusterSigns, String sequenceName){
		//SampleStatistics calcEnd = new SampleStatistics();
		//for(GenomicVariant v:clusterSigns) calcEnd.update(v.getLast());;
		GenomicVariant first = clusterSigns.get(0);
		GenomicVariant last = clusterSigns.get(clusterSigns.size()-1);
		int begin = first.getFirst();
		int end = last.getLast();
		byte type = first.getType();
		short variantScore = calculateVariantScore(clusterSigns);
		GenomicVariantImpl variant = new GenomicVariantImpl(sequenceName, begin, end, type);
		variant.setVariantQS(variantScore);
		if(GenomicVariantImpl.TYPE_LARGEINS == variant.getType()) variant.setLast(begin+1);
		return variant;
	}
	public void filterVariantsBySimilarity(GenomicRegionSortedCollection<GenomicVariant> sortedVariants) {
		List<GenomicVariant> variantsToKeep = new ArrayList<>();
		List<GenomicVariant> variantsList = sortedVariants.asList();
		List<GenomicVariant> localCandidates = new ArrayList<>();
		for(int i = 0; i < sortedVariants.size() - 1; i++) {
			GenomicVariant current = variantsList.get(i);
			GenomicVariant next = variantsList.get(i+1);
			if(testCandidateCompatibility(current, next)) {
				localCandidates.add(current);
				localCandidates.add(next);
			}else {
				if(!localCandidates.isEmpty()) {
					Collections.sort(localCandidates, Comparator.comparingInt(v -> (int) v.getVariantQS()));
					variantsToKeep.add(localCandidates.get(localCandidates.size() - 1));
					localCandidates.clear();
				}else {
					variantsToKeep.add(current);
				}
			}
		}
		sortedVariants.retainAll(variantsToKeep);
		sortedVariants.forceSort();
	}
	@Override
	public GenomicRegionSortedCollection<GenomicVariant> callVariants(){
		List<GenomicVariant> variants = new ArrayList<>();
		Map<String, List<List<Integer>>> clusters = findVariantClusters();
		List<String> keys = new ArrayList<>(clusters.keySet());
		QualifiedSequenceList sequences = new QualifiedSequenceList();
		for(String k:keys) sequences.add(new QualifiedSequence(k));
		GenomicRegionSortedCollection<GenomicVariant> sortedVariants = new GenomicRegionSortedCollection<>(sequences);
		for(String k:keys) {
			List<List<Integer>> chrClusters = clusters.get(k);
			for(List<Integer> cluster:chrClusters) {
				List<GenomicVariant> clusterSigns = new ArrayList<>();
				for(int idx:cluster) clusterSigns.add(signatures.getSequenceRegions(k).asList().get(idx));
				Collections.sort(clusterSigns, Comparator.comparingInt(s -> s.getFirst()));
				System.out.println("Cluster with size: " + clusterSigns.size());
				/**for (GenomicVariant sign:clusterSigns) {
					System.out.println("Signature: " + sign.getSequenceName() + " begin: " + sign.getFirst() + " end: " +
							sign.getLast() + " length: " + sign.length() + " type: " + 
							GenomicVariantImpl.getVariantTypeName(sign.getType()));
				}**/
				GenomicVariantImpl variant = processClusterToVariant(clusterSigns, k);
				if(variant.length() >= 30) variants.add(variant);
			}
		}
		sortedVariants.addAll(variants);
		sortedVariants.forceSort();
		filterVariantsBySimilarity(sortedVariants);
		return sortedVariants;
	}
}
