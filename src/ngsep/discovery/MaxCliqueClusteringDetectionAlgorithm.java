package ngsep.discovery;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import JSci.maths.statistics.SampleStatistics;
import ngsep.genome.GenomicRegion;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.MaximalCliquesFinder;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;

public class MaxCliqueClusteringDetectionAlgorithm implements LongReadVariantCallerAlgorithm{
	
	public static final double DEFAULT_PD_NORM_FACTOR = 900;
	public static final double DEFAULT_EDGE_TRESHHOLD = 0.7;
	private double pdNormFactor = DEFAULT_PD_NORM_FACTOR;
	private double edgeTreshold = DEFAULT_EDGE_TRESHHOLD;
	
	private ReferenceGenome refGenome;
	private Map<String, List<Signature>> signatures;
	private int partitionSize = 200;
	
	@Override
	public void setPartitionSize(int partitions) {
		// TODO Auto-generated method stub
		this.partitionSize = partitions;
	}
	@Override
	public void setSignatures(Map<String, List<Signature>> signatures) {
		// TODO Auto-generated method stub
		this.signatures = signatures;
	}
	public void setReferenceGenome(ReferenceGenome ref) {
		// TODO Auto-generated method stub
		this.refGenome = ref;
	}
	public double calculateSPD(Signature sign1, Signature sign2) {
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

	public boolean[][] calculateAdjacencyMatrix(List<Signature> candidateSignatures) {
		int n = candidateSignatures.size();
		boolean[][] adjacencyMatrix = new boolean[n][n];
		for(int i = 0; i < n; i++) {
			Signature si = candidateSignatures.get(i);
			for(int j = 0; j < n; j++) {
				Signature sj = candidateSignatures.get(j);
				double spd = calculateSPD(si, sj);
				if(spd < edgeTreshold && !si.equals(sj)) {
					adjacencyMatrix[i][j] = true;
				}
			}
		}
		return adjacencyMatrix;
	}

	public Map<String, List<List<Integer>>> findVariantClusters() {
		Map<String, List<List<Integer>>> clusters = new LinkedHashMap<>();
		int begin = 0;
		int end = 0;
		List<String> keys = new ArrayList<String>(signatures.keySet());
		for(String k:keys){
			List<Signature> signList = signatures.get(k);
			List<List<Integer>> chrClusters = new ArrayList<>();
			int n = signList.size();
			int partSize = n/partitionSize;
			int modSize = n%partitionSize;
			begin = 0;
			end = partSize;
			for(int i = 0; i < partitionSize; i++) {
				List<Signature> partition = signList.subList(begin, end);
				List<Signature> indel = new ArrayList<>();
				List<Signature> del = new ArrayList<>();
				List<Signature> ins = new ArrayList<>();
				List<Signature> inv = new ArrayList<>();
				List<Signature> dup = new ArrayList<>();
				List<Signature> tan = new ArrayList<>();
				List<Signature> brk = new ArrayList<>();
				List<Signature> und = new ArrayList<>();
				List<List<Signature>> signsPerType = new ArrayList<>();
				for(Signature s:partition) {
					String type = s.getType();
					if(Signature.TYPENAME_INDEL.equals(type)) indel.add(s);
					else if(Signature.TYPENAME_DELETION.equals(type)) del.add(s);
					else if(Signature.TYPENAME_INSERTION.equals(type)) ins.add(s);
					else if(Signature.TYPENAME_DUPLICATION.equals(type)) dup.add(s);
					else if(Signature.TYPENAME_INVERSION.equals(type)) inv.add(s);
					else if(Signature.TYPENAME_TANDEM.equals(type)) tan.add(s);
					else if(Signature.TYPENAME_TRANSLOCATION_BREAKPOINT.equals(type)) brk.add(s);
					else if(Signature.TYPENAME_UNDETERMINED.equals(type)) und.add(s);
				}
				signsPerType.add(indel);
				signsPerType.add(del);
				signsPerType.add(ins);
				signsPerType.add(dup);
				signsPerType.add(inv);
				signsPerType.add(tan);
				signsPerType.add(brk);
				signsPerType.add(und);
				for(List<Signature> typePartition:signsPerType) {
					List<Integer> idxs = new ArrayList<>();
					if(typePartition.isEmpty()) continue;
					for(Signature sign:typePartition) idxs.add(signList.indexOf(sign));
					boolean[][] adjMatrix = calculateAdjacencyMatrix(typePartition);
					System.out.println("Processing partition with types: " + typePartition.get(0).getType() + " from: " + begin + " to: " + end);
					ArrayList<ArrayList<Integer>> maxCliques = MaximalCliquesFinder.callMaxCliqueFinder(idxs, adjMatrix);
					//List<List<Integer>> cliques = CliquesFinder.findCliques(adjMatrix);
					System.out.println("Cliques found");
					chrClusters.addAll(maxCliques);
					}
				begin = end + 1;
				end += partSize;
				if(i + 1 == partitionSize - 1) end += modSize;
			}
			clusters.put(k, chrClusters);
		}
		return clusters;
	}
	
	public short calculateVariantScore(List<Signature> candidates) {
		double score = 0;
		SampleStatistics calcSpan = new SampleStatistics();
		SampleStatistics calcPos = new SampleStatistics();
		for(Signature sign:candidates) {
			calcPos.update(sign.getFirst());
			calcSpan.update(sign.length());
		}
		double normPos = 10*(1-Math.min(1, 
						(Math.sqrt(calcPos.getVariance()))/calcSpan.getMean()));
		double normSpan = 20*(1-Math.min(1, 
						(Math.sqrt(calcSpan.getVariance()))/calcSpan.getMean()));
		double normNum = candidates.size() >= 40 ? 40:40*((double) candidates.size()/40);
		score = normPos + normSpan + normNum;
		score = 100*(score/70);
		return (short) score;
	}
	
	public boolean testCandidateCompatibility(GenomicVariant v1, GenomicVariant v2) {
		if(v1.getFirst() == v2.getFirst() || 
			v1.getLast() == v2.getLast() ||
			v1.length() == v2.length()) return true;
		return false;
	}
	
	@Override
	public List<GenomicVariant> callVariants(){
		List<GenomicVariant> variants = new ArrayList<>();
		Map<String, List<List<Integer>>> clusters = findVariantClusters();
		List<String> keys = new ArrayList<>(clusters.keySet());
		for(String k:keys) {
			List<List<Integer>> chrClusters = clusters.get(k);
			for(List<Integer> cluster:chrClusters) {
				List<Signature> clusterSigns = new ArrayList<>();
				for(int idx:cluster) clusterSigns.add(signatures.get(k).get(idx));
				Collections.sort(clusterSigns, Comparator.comparingInt(s -> s.getFirst()));
				Signature first = clusterSigns.get(0);
				Signature last = clusterSigns.get(clusterSigns.size()-1);
				int begin = first.getFirst();
				int end = last.getLast();
				byte type = GenomicVariantImpl.getVariantTypeId(first.getType());
				short variantScore = calculateVariantScore(clusterSigns);
				/**CharSequence refAllele = refGenome.getReference(k,
						begin, end);
				List<String> alleles = new ArrayList<>();
				alleles.add(refAllele.toString());
				alleles.add("ALT");
				GenomicVariantImpl variant = new GenomicVariantImpl(k, begin, end, alleles);
				variant.setType(type);
				**/
				GenomicVariantImpl variant = new GenomicVariantImpl(k, begin, end, type);
				variant.setVariantQS(variantScore);
				if(GenomicVariantImpl.TYPE_LARGEINS == variant.getType()) variant.setLast(begin+1);
				if(!variants.isEmpty()) {
					int n = variants.size()-1;
					GenomicVariant previous = variants.get(n);
					if(testCandidateCompatibility(variant, previous)){
						if(variantScore >= previous.getVariantQS()) {
							variants.set(n, variant);
						}
					}
					else variants.add(variant);
				}
				else variants.add(variant);
			}
		}
		return variants;
	}
	
}
