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

public class MaxCliqueClusteringDetectionAlgorithm implements LongReadVariantDetectorAlgorithm{
	
	public static final double DEFAULT_PD_NORM_FACTOR = 900;
	public static final double DEFAULT_EDGE_TRESHOLD = 0.7;
	public static final double MAX_DOWNSTREAM_CONSENSUS_LENGTH = 50;
	
	private double pdNormFactor = DEFAULT_PD_NORM_FACTOR;
	private double edgeTreshold = DEFAULT_EDGE_TRESHOLD;
	private GenomicRegionSortedCollection<GenomicVariant> signatures;
	private ReferenceGenome refGenome;
	private int lengthToDefineSVEvent;
	
	public MaxCliqueClusteringDetectionAlgorithm(ReferenceGenome genome, 
			GenomicRegionSortedCollection<GenomicVariant>  signatures, int SVEventLength){
		this.refGenome = genome;
		this.signatures = signatures;
		this.lengthToDefineSVEvent = SVEventLength;
	}
	@Override
	public GenomicRegionSortedCollection<GenomicVariant> callVariants(){
		List<GenomicVariant> variants = new ArrayList<>();
		Map<String, List<List<Integer>>> clusters = findVariantClusters();
		List<String> keys = new ArrayList<>(clusters.keySet());
		QualifiedSequenceList sequences = new QualifiedSequenceList(refGenome.getSequencesList());
		GenomicRegionSortedCollection<GenomicVariant> sortedVariants = new GenomicRegionSortedCollection<>(sequences);
		for(String k:keys) {
			List<List<Integer>> chrClusters = clusters.get(k);
			List<GenomicVariant> chrSignList = signatures.getSequenceRegions(k).asList();
			for(List<Integer> cluster:chrClusters) {
				List<GenomicVariant> clusterSigns = new ArrayList<>();
				for(int idx:cluster) clusterSigns.add(chrSignList.get(idx));
				Collections.sort(clusterSigns, Comparator.comparingInt(s -> s.getFirst()));
				/**
				System.out.println("$Cluster with size: " + clusterSigns.size());
				for (GenomicRegion sign:clusterSigns) {
					System.out.println("$Signature: " + sign.getSequenceName() + " begin: " + sign.getFirst() + " end: " +
							sign.getLast() + " length: " + sign.length() + " type: " + 
							GenomicVariantImpl.getVariantTypeName(getSignatureType(sign)));
				}**/
				GenomicVariantImpl variant = processClusterToVariant(clusterSigns, k);
				if(variant.length() >= lengthToDefineSVEvent) variants.add(variant);
			}
		}
		sortedVariants.addAll(variants);
		sortedVariants.forceSort();
		filterVariantsBySimilarity(sortedVariants);
		//System.out.println("Variants filtered");
		return sortedVariants;
	}
	public void filterVariantsBySimilarity(GenomicRegionSortedCollection<GenomicVariant> sortedVariants) {
		List<GenomicVariant> variantsToKeep = new ArrayList<>();
		List<GenomicVariant> variantsList = sortedVariants.asList();
		List<GenomicVariant> localCandidates = new ArrayList<>();
		List<List<GenomicVariant>> variantsPerType = new ArrayList<>();
		GenomicRegionComparator cmpClassInstance = new GenomicRegionComparator(sortedVariants.getSequenceNames());
		List<GenomicVariant> indel = new ArrayList<>();
		List<GenomicVariant> del = new ArrayList<>();
		List<GenomicVariant> ins = new ArrayList<>();
		List<GenomicVariant> inv = new ArrayList<>();
		List<GenomicVariant> dup = new ArrayList<>();
		List<GenomicVariant> tan = new ArrayList<>();
		List<GenomicVariant> und = new ArrayList<>();
		for(int i = 0; i < sortedVariants.size(); i++) {
			GenomicVariant current = variantsList.get(i);
			byte type = current.getType();
			if(GenomicVariant.TYPE_INDEL == type) indel.add(current);
			else if(GenomicVariant.TYPE_LARGEDEL == type) del.add(current);
			else if(GenomicVariant.TYPE_LARGEINS == type) ins.add(current);
			else if(GenomicVariant.TYPE_DUPLICATION  == type) dup.add(current);
			else if(GenomicVariant.TYPE_INVERSION  == type) inv.add(current);
			else if(GenomicVariant.TYPE_CNV == type) tan.add(current);
			else if(GenomicVariant.TYPE_UNDETERMINED == type) und.add(current);
		}
		variantsPerType.add(ins);
		variantsPerType.add(del);
		variantsPerType.add(dup);
		variantsPerType.add(inv);
		variantsPerType.add(tan);
		variantsPerType.add(indel);
		variantsPerType.add(und);
		for(List<GenomicVariant> typedVariantList : variantsPerType) {
			if(typedVariantList.isEmpty()) continue;
			for(int i = 0; i < typedVariantList.size() - 1; i++) {
				GenomicVariant current = typedVariantList.get(i);
				GenomicVariant next = typedVariantList.get(i+1);
				int cmp;
				if(current.getType() == GenomicVariant.TYPE_LARGEINS) {
					cmp = modifiedInsComparison(current, next, cmpClassInstance);
				}
				else cmp = cmpClassInstance.compare(current, next);
				if(cmp == -1 || cmp == 0 || cmp == 1) {
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
		}
		sortedVariants.retainAll(variantsToKeep);
		sortedVariants.forceSort();
	}
	public GenomicVariantImpl processClusterToVariant(List<GenomicVariant> clusterSigns, String sequenceName){
		GenomicVariant firstVar = clusterSigns.get(0);
		SampleStatistics calcFirst = new SampleStatistics();
		SampleStatistics calcLast = new SampleStatistics();
		int endOfSpan = -1;
		detectAndRemoveOutliers(clusterSigns);
		for(GenomicVariant v:clusterSigns) {
			int currentEndOfSpan = v.getFirst() + v.length() - 1;
			calcFirst.update(v.getFirst());
			calcLast.update(currentEndOfSpan);
			//if(currentEndOfSpan > endOfSpan) endOfSpan = currentEndOfSpan;
		}
		int first = (int) calcFirst.getMean();
		int last = (int) calcLast.getMean();
		endOfSpan = last;
		//System.out.println("first: " + first);
		//System.out.println("last: " + last);
		byte type = firstVar.getType();
		short variantScore = calculateVariantScore(clusterSigns);
		List<String> alleles = new ArrayList<>();
		String ref = "";
		String alt = "";
		if(GenomicVariant.TYPE_LARGEINS == type) {
			last = first+1;
			ref = "" + refGenome.getReferenceBase(sequenceName, first);
			alt = "<" + GenomicVariant.TYPENAME_LARGEINS + ">";
		}
		else if(GenomicVariant.TYPE_LARGEDEL == type) {
			ref = "" + refGenome.getReferenceBase(sequenceName, first);
			alt = "<" + GenomicVariant.TYPENAME_LARGEDEL + ">";
		}
		else if(GenomicVariant.TYPE_INVERSION == type) {
			ref = "" + refGenome.getReferenceBase(sequenceName, first);
			alt = "<" + GenomicVariant.TYPENAME_INVERSION + ">";
		}
		alleles.add(ref);
		alleles.add(alt);
		GenomicVariantImpl variant = new GenomicVariantImpl(sequenceName, first, alleles);
		variant.setLength(endOfSpan - first + 1);
		variant.setLast(last);
		variant.setType(type);
		variant.setVariantQS(variantScore);
		return variant;
	}
	private void detectAndRemoveOutliers(List<GenomicVariant> clusterSigns) {
		// TODO Auto-generated method stub
		List<GenomicVariant> outliersToRemove = new ArrayList<>();
		int n = clusterSigns.size();
		double [] pdValues = new double[n];
		SampleStatistics clusterCalc = new SampleStatistics();
		double spd;
		SampleStatistics signCalc = new SampleStatistics();
		for(int i = 0; i < n; i++) {
			GenomicVariant si = clusterSigns.get(i);
			signCalc = new SampleStatistics();
			for(int j = 0; j < n; j++) {
				GenomicVariant sj = clusterSigns.get(j);
				if(si.equals(sj)) continue;
				else spd = calculateSD(si.length(), sj.length());
				signCalc.update(spd);
			}
			double rowAvg = signCalc.getMean();
			clusterCalc.update(rowAvg);
			pdValues[i] = rowAvg;
		}
		double averageSPD = clusterCalc.getMean();
		double stdSPD = Math.sqrt(clusterCalc.getVariance());
		double limit = 2*stdSPD;
		for(int i = 0; i < n; i++) {
			double distance = Math.abs(pdValues[i] - averageSPD);
			if(distance > limit) {
				outliersToRemove.add(clusterSigns.get(i));
			}
		}
		clusterSigns.removeAll(outliersToRemove);
	}
	public Map<String, List<List<Integer>>> findVariantClusters() {
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
					if(!testDownstreamSignatureCompatibility(current,next) || toCluster.size() >= 300) {
						List<Integer> idxs = new ArrayList<>();
						for(GenomicVariant sign:toCluster) idxs.add(signList.indexOf(sign));
						if(toCluster.size() < 4) {
							toCluster.clear();
							continue;
						}
						else{
							boolean[][] adjMatrix = calculateAdjacencyMatrix(toCluster);
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
	public double calculateSPD(GenomicVariant sign1, GenomicVariant sign2) {
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
		SPD = SD + (double) PD/pdNormFactor;
		return SPD;
	}
	public double calculateSD(int span1, int span2) {
		double SD = Math.abs((span1)-(span2));
		SD = (double) SD/Math.max((span1), (span2));
		return SD;
	}
	public int calculatePD(int first1, int first2, int last1, int last2) {
		int PD = Math.min(Math.abs(first1-first2), Math.abs(last1-last2));
		PD = Math.min(PD, Math.abs(((first1-last1)/2)) - ((first2-last2)/2));
		return PD;
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
	public short calculateVariantScore(List<GenomicVariant> candidates) {
		double score = 0;
		SampleStatistics calcPos = new SampleStatistics();
		SampleStatistics calcSpan = new SampleStatistics();
		for(GenomicRegion sign:candidates) {
			calcPos.update(sign.getFirst());
			calcSpan.update(sign.length());
		}
		int n = candidates.size();
		double numSign = Math.min(80, n);
		double spanMean = calcSpan.getMean();
		double normPos = (numSign/8)*(1-Math.min(1, 
						(Math.sqrt(calcPos.getVariance()))/spanMean));
		double normSpan = (numSign/8)*(1-Math.min(1, 
						(Math.sqrt(calcSpan.getVariance()))/spanMean));
		score = normSpan + normPos + numSign;
		return (short) score;
	}
	private int modifiedInsComparison(GenomicVariant current, GenomicVariant next,
			GenomicRegionComparator cmpClassInstance) {
		// TODO Auto-generated method stub
		int tempCurrEnd = current.getFirst() + current.length();
		int tempNextEnd = next.getFirst() + next.length();
		GenomicVariant tempCurrent = new GenomicVariantImpl(current.getSequenceName(), current.getFirst(),
				tempCurrEnd, GenomicVariant.TYPE_LARGEINS);
		GenomicVariant tempNext = new GenomicVariantImpl(next.getSequenceName(), next.getFirst(),
				tempNextEnd, GenomicVariant.TYPE_LARGEINS);
		return cmpClassInstance.compare(tempCurrent, tempNext);
	}
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
