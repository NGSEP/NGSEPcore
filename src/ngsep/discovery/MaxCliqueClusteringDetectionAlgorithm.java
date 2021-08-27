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
	private double pdNormFactor = DEFAULT_PD_NORM_FACTOR;
	private double edgeTreshold = DEFAULT_EDGE_TRESHOLD;
	
	private ReferenceGenome refGenome;
	private Map<String, List<GenomicVariant>> signatures;
	private int partitionSignatureSize = 10;
	
	@Override
	public void setSignatures(Map<String, List<GenomicVariant>> signatures) {
		// TODO Auto-generated method stub
		this.signatures = signatures;
	}
	public void setReferenceGenome(ReferenceGenome ref) {
		// TODO Auto-generated method stub
		this.refGenome = ref;
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
		double normPos = 10.0*(1-Math.min(1, 
						(Math.sqrt(calcPos.getVariance()))/calcSpan.getMean()));
		double normSpan = 20.0*(1-Math.min(1, 
						(Math.sqrt(calcSpan.getVariance()))/calcSpan.getMean()));
		double normNum = candidates.size() >= 40 ? 40:40*((double) candidates.size()/40);
		score = normPos + normSpan + normNum;
		score = 100*(score/70);
		return (short) score;
	}
	public boolean testCandidateCompatibility(GenomicVariant v1, GenomicVariant v2) {
		if(!v1.getSequenceName().equals(v2.getSequenceName())) return false;
		if(!(v1.getType() == v2.getType())) return false;
		if(Math.abs(v1.getFirst() - v2.getFirst()) < 1000) return true; 
		if(Math.abs(v1.getLast() - v2.getLast()) < 1000) return true; 
		if(Math.abs(v1.length() - v2.length()) < 1000) return true; 
		return false;
	}
	/**
	public Map<String, List<List<Integer>>> findVariantClusters() {
		Map<String, List<List<Integer>>> clusters = new LinkedHashMap<>();
		List<String> keys = new ArrayList<String>(signatures.keySet());
		for(String k:keys){
			List<GenomicVariant> signList = signatures.get(k);
			List<List<Integer>> chrClusters = new ArrayList<>();
			List<GenomicVariant> indel = new ArrayList<>();
			List<GenomicVariant> del = new ArrayList<>();
			List<GenomicVariant> ins = new ArrayList<>();
			List<GenomicVariant> inv = new ArrayList<>();
			List<GenomicVariant> dup = new ArrayList<>();
			List<GenomicVariant> tan = new ArrayList<>();
			List<GenomicVariant> und = new ArrayList<>();
			List<List<GenomicVariant>> signsPerType = new ArrayList<>();
			for(int i = 0; i < signList.size(); i++) {
				GenomicVariant current = signList.get(i);
				String type = GenomicVariantImpl.getVariantTypeName(current.getType());
				if(GenomicVariant.TYPENAME_INDEL.equals(type)) indel.add(current);
				else if(GenomicVariant.TYPENAME_LARGEDEL.equals(type)) del.add(current);
				else if(GenomicVariant.TYPENAME_LARGEINS.equals(type)) ins.add(current);
				else if(GenomicVariant.TYPENAME_DUPLICATION.equals(type)) dup.add(current);
				else if(GenomicVariant.TYPENAME_INVERSION.equals(type)) inv.add(current);
				else if(GenomicVariant.TYPENAME_CNV.equals(type)) tan.add(current);
				else if(GenomicVariant.TYPENAME_UNDETERMINED.equals(type)) und.add(current);
				if(i+1==signList.size() || !testCandidateCompatibility(current, signList.get(i+1))) {
					signsPerType.add(indel);
					signsPerType.add(del);
					signsPerType.add(ins);
					signsPerType.add(dup);
					signsPerType.add(inv);
					signsPerType.add(tan);
					signsPerType.add(und);
					for(List<GenomicVariant> typePartition:signsPerType) {
						List<Integer> idxs = new ArrayList<>();
						if(typePartition.isEmpty()) continue;
						for(GenomicVariant sign:typePartition) idxs.add(signList.indexOf(sign));
						boolean[][] adjMatrix = calculateAdjacencyMatrix(typePartition);
						System.out.println("#partition size: " + typePartition.size());
						System.out.println("Processing partition with types: " + typePartition.get(0).getType() + " from: "
								+ idxs.get(0) + " to: " + idxs.get(idxs.size()-1));
						ArrayList<ArrayList<Integer>> maxCliques = MaximalCliquesFinder.callMaxCliqueFinder(idxs, adjMatrix);
						//List<List<Integer>> cliques = CliquesFinder.findCliques(adjMatrix);
						System.out.println("Cliques found");
						chrClusters.addAll(maxCliques);
						}
					signsPerType.clear();
					indel.clear();
					del.clear();
					ins.clear();
					dup.clear();
					inv.clear();
					tan.clear();
					und.clear();
				}
			}
			clusters.put(k, chrClusters);
		}
		return clusters;
	}
		**/
	public Map<String, List<List<Integer>>> findVariantClusters() {
		Map<String, List<List<Integer>>> clusters = new LinkedHashMap<>();
		List<String> keys = new ArrayList<String>(signatures.keySet());
		for(String k:keys){
			List<GenomicVariant> signList = signatures.get(k);
			List<List<Integer>> chrClusters = new ArrayList<>();
			List<GenomicVariant> indel = new ArrayList<>();
			List<GenomicVariant> del = new ArrayList<>();
			List<GenomicVariant> ins = new ArrayList<>();
			List<GenomicVariant> inv = new ArrayList<>();
			List<GenomicVariant> dup = new ArrayList<>();
			List<GenomicVariant> tan = new ArrayList<>();
			List<GenomicVariant> und = new ArrayList<>();
			List<List<GenomicVariant>> signsPerType = new ArrayList<>();
			int j = 0;
			for(int i = 0; i < signList.size(); i++) {
				GenomicVariant current = signList.get(i);
				String type = GenomicVariantImpl.getVariantTypeName(current.getType());
				if(GenomicVariant.TYPENAME_INDEL.equals(type)) indel.add(current);
				else if(GenomicVariant.TYPENAME_LARGEDEL.equals(type)) del.add(current);
				else if(GenomicVariant.TYPENAME_LARGEINS.equals(type)) ins.add(current);
				else if(GenomicVariant.TYPENAME_DUPLICATION.equals(type)) dup.add(current);
				else if(GenomicVariant.TYPENAME_INVERSION.equals(type)) inv.add(current);
				else if(GenomicVariant.TYPENAME_CNV.equals(type)) tan.add(current);
				else if(GenomicVariant.TYPENAME_UNDETERMINED.equals(type)) und.add(current);
			    if(j == partitionSignatureSize || i+1==signList.size()) {
					signsPerType.add(indel);
					signsPerType.add(del);
					signsPerType.add(ins);
					signsPerType.add(dup);
					signsPerType.add(inv);
					signsPerType.add(tan);
					signsPerType.add(und);
					for(List<GenomicVariant> typePartition:signsPerType) {
						List<Integer> idxs = new ArrayList<>();
						if(typePartition.isEmpty()) continue;
						for(GenomicVariant sign:typePartition) idxs.add(signList.indexOf(sign));
						boolean[][] adjMatrix = calculateAdjacencyMatrix(typePartition);
						System.out.println("#partition size: " + typePartition.size());
						System.out.println("Processing partition with types: " + typePartition.get(0).getType() + " from: "
								+ idxs.get(0) + " to: " + idxs.get(idxs.size()-1));
						ArrayList<ArrayList<Integer>> maxCliques = MaximalCliquesFinder.callMaxCliqueFinder(idxs, adjMatrix);
						//List<List<Integer>> cliques = CliquesFinder.findCliques(adjMatrix);
						System.out.println("Cliques found");
						chrClusters.addAll(maxCliques);
						}
					signsPerType.clear();
					indel.clear();
					del.clear();
					ins.clear();
					dup.clear();
					inv.clear();
					tan.clear();
					und.clear();
					j = 0;
				} else j++;
			}
			clusters.put(k, chrClusters);
		}
		return clusters;
	}
	public GenomicVariantImpl processClusterToVariant(List<GenomicVariant> clusterSigns, String sequenceName){
		Collections.sort(clusterSigns, Comparator.comparingInt(s -> s.getFirst()));
		GenomicVariant first = clusterSigns.get(0);
		GenomicVariant last = clusterSigns.get(clusterSigns.size()-1);
		int begin = first.getFirst();
		int end = last.getLast();
		byte type = first.getType();
		short variantScore = calculateVariantScore(clusterSigns);
		/**CharSequence refAllele = refGenome.getReference(k,
				begin, end);
		List<String> alleles = new ArrayList<>();
		alleles.add(refAllele.toString());
		alleles.add("ALT");
		GenomicVariantImpl variant = new GenomicVariantImpl(k, begin, end, alleles);
		variant.setType(type);
		**/
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
				for(int idx:cluster) clusterSigns.add(signatures.get(k).get(idx));
				GenomicVariantImpl variant = processClusterToVariant(clusterSigns, k);
				variants.add(variant);
			}
		}
		sortedVariants.addAll(variants);
		sortedVariants.forceSort();
		filterVariantsBySimilarity(sortedVariants);
		return sortedVariants;
	}
}
