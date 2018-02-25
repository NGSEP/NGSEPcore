package ngsep.discovery;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.SequenceDistanceMeasure;

public class AlleleCallClustersBuilder {

	private Map<Integer, List<String>> alleleCallsByLength = new HashMap<>();
	private int totalAlleleCalls = 0;
	private SequenceDistanceMeasure distanceMeasure = new HammingSequenceDistanceMeasure();
	private double minRelativeProportion = 0.2; 
	
	public void addAlleleCall (String call) {
		int l = call.length();
		List<String> allelesLength = alleleCallsByLength.get(l);
		if(allelesLength==null) {
			allelesLength = new ArrayList<>();
			alleleCallsByLength.put(l, allelesLength);
		}
		allelesLength.add(call);
		totalAlleleCalls++;
	}
	public Map<String,List<String>> clusterAlleleCalls (String [] suggestedAlleles, boolean allowNewAlleles) {
		Map<String,List<String>> alleleClusters = new TreeMap<>();
		Map<Integer, List<String>> filteredClusters = filterLengthClusters();
		int s = filteredClusters.size();
		if(s==0) return alleleClusters;
		for(int l:filteredClusters.keySet()) {
			List<String> callsL = filteredClusters.get(l);	
			Set<String> suggestedAllelesSet = getAllelesByLength(suggestedAlleles,l);
			Map<String,List<String>> lengthClusters;
			if(allowNewAlleles) {
				if(callsL.size()<5*suggestedAllelesSet.size()) {
					//With low coverage and suggested alleles, only those are taken into account
					lengthClusters = clusterAlleleCallsPivotAlleles(callsL,suggestedAllelesSet);
				} else {
					//With enough calls, suggested alleles are actually used as a suggestion and the consensus is considered
					String consensus = HammingSequenceDistanceMeasure.makeHammingConsensus(callsL);
					suggestedAllelesSet.add(consensus);
					if(suggestedAllelesSet.size() == 1 && callsL.size()<5) {
						lengthClusters = clusterAlleleCallsPivotAlleles(callsL,suggestedAllelesSet);
					} else {
						lengthClusters = clusterAlleleCallsByHammingDistance(callsL, suggestedAllelesSet);
					}
					
				}
			} else if (suggestedAllelesSet.size()>0) {
				lengthClusters = clusterAlleleCallsPivotAlleles(callsL,suggestedAllelesSet);
			} else continue;
			alleleClusters.putAll(lengthClusters);
		}
		return alleleClusters;
	}
	private Set<String> getAllelesByLength(String[] suggestedAlleles, int l) {
		Set<String> answer = new TreeSet<>();
		for(int i=0;i<suggestedAlleles.length;i++) {
			if(suggestedAlleles[i].length()==l) answer.add(suggestedAlleles[i]); 
		}
		return answer;
	}
	/**
	 * Filter severe disbalances in read count for clusters based on length. Useful to remove stutter products
	 * @return Map<Integer, List<String>> Length clusters passing the minimum proportion of reads
	 */
	private Map<Integer, List<String>> filterLengthClusters() {
		Map<Integer, List<String>> answer = new HashMap<>();
		double minCount  = minRelativeProportion*totalAlleleCalls; 
		for(int l:alleleCallsByLength.keySet()) {
			List<String> allelesL = alleleCallsByLength.get(l);
			if(minCount<=allelesL.size()) answer.put(l, allelesL);
		}
		return answer;
	}
	/**
	 * 
	 * @param calls
	 * @param suggestedAlleles PRE: Suggested alleles is not empty
	 * @return
	 */
	private Map<String, List<String>> clusterAlleleCallsPivotAlleles(List<String> calls, Set<String> suggestedAlleles) {
		assert suggestedAlleles.size() > 0;
		Map<String, List<String>> clusters = new TreeMap<>();
		for(String call:calls) {
			//Find closest allele
			double minDistance = 1.1;
			String minAllele = null;
			for(String allele:suggestedAlleles) {
				double d = distanceMeasure.calculateNormalizedDistance(call, allele);
				if(minAllele==null || minDistance>d ) {
					minAllele = allele;
					minDistance = d;
				}
			}
			if(minDistance<=0.5) {
				//If not too far from the best allele, add to cluster
				List<String> cluster = clusters.get(minAllele);
				if(cluster==null) {
					cluster = new ArrayList<>();
					clusters.put(minAllele, cluster);
				}
				cluster.add(call);
			}
		}	
		return clusters;
	}
	private Map<String, List<String>> clusterAlleleCallsByHammingDistance(List<String> calls, Set<String> suggestedAlleles) {
		Map<String,List<String>> answer = new TreeMap<>();
		if(calls.size()==0) return answer;
		
		List<String> allSequences = new ArrayList<>();
		allSequences.addAll(suggestedAlleles);
		allSequences.addAll(calls);
		double [] relativeFrequencies = HammingSequenceDistanceMeasure.calculateMinorRelativeFrequencies(allSequences);
		
		List<Integer> variantSites = new ArrayList<>();
		for(int i=0;i<relativeFrequencies.length;i++) {
			if(relativeFrequencies[i]>=minRelativeProportion) {
				variantSites.add(i);
			}
		}
		if(variantSites.size()==0) {
			//No real variation
			answer.put(allSequences.get(0), calls);
			return answer;
		}
		return clusterAlleleCallsByHaplotypes (allSequences,variantSites);
	}
	
	private Map<String, List<String>> clusterAlleleCallsByHaplotypes(List<String> sequences, List<Integer> variantSites) {
		assert sequences.size()>0;
		int l = sequences.get(0).length();
		for(CharSequence seq:sequences) assert seq.length()==l;
		
		//Build haplotypes
		StringBuilder [] haplotypes = new StringBuilder[sequences.size()];
		for(int i=0;i<haplotypes.length;i++) {
			haplotypes[i] = new StringBuilder();
		}
		for(Integer j:variantSites) {
			for(int i=0;i<haplotypes.length;i++) {
				CharSequence seq = sequences.get(i);
				haplotypes[i].append(seq.charAt(j));
			}
		}
		//Cluster sequences having the same haplotype 
		//TODO: Allow mismatches in the haplotypes
		Map<String,List<String>> haplotypesMap = new TreeMap<>();
		for(int i=0;i<haplotypes.length;i++) {
			String hap = haplotypes[i].toString();
			List<String> sequencesHap = haplotypesMap.get(hap);
			if(sequencesHap==null) {
				sequencesHap = new ArrayList<>();
				haplotypesMap.put(hap, sequencesHap);
			}
			sequencesHap.add(sequences.get(i));
		}
		Map<String, List<String>> answer = new TreeMap<>();
		for(List<String> cluster:haplotypesMap.values()) {
			answer.put(cluster.get(0), cluster);
		}
		
		return answer;
	}
	
}
