package ngsep.discovery;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.math.CountsRankHelper;

public class AlleleClustersBuilder {

	private Map<Integer, List<String>> allelesByLength = new HashMap<>();
	
	public void addAllele (String allele) {
		int l = allele.length();
		List<String> allelesLength = allelesByLength.get(l);
		if(allelesLength==null) {
			allelesLength = new ArrayList<>();
			allelesByLength.put(l, allelesLength);
		}
		allelesLength.add(allele);
	}
	public Map<String,Integer> buildAlleleClusters (String referenceAllele) {
		Map<String,Integer> alleleClusters = new TreeMap<>();
		int s = allelesByLength.size();
		if(s==0) return alleleClusters;
		if(s==1) {
			//Make clusters by differences with consensus
			//TODO: Make it multiallelic
			List<String> allelesL = allelesByLength.values().iterator().next();
			String consensus = makeHammingConsensus(allelesL);
			if (consensus.length()!=referenceAllele.length() ) {
				return buildClustersByHammingDistance(allelesL, consensus);
			} else if(!consensus.equals(referenceAllele)) {
				return buildTwoClusters(allelesL, referenceAllele,consensus);
			} else {
				return buildClustersByHammingDistance(allelesL, referenceAllele);
			}
			
		}
		for(int l:allelesByLength.keySet()) {
			List<String> allelesL = allelesByLength.get(l);
			String consensus = makeHammingConsensus(allelesL);
			alleleClusters.put(consensus, allelesL.size());
		}
		return alleleClusters;
	}
	private String makeHammingConsensus(List<String> alleles) {
		int l = alleles.get(0).length();
		StringBuilder sb = new StringBuilder();
		for(int i=0;i<l;i++) {
			CountsRankHelper<Character> countsChar = new CountsRankHelper<>();
			for(String allele:alleles) {
				countsChar.add(allele.charAt(i));
			}
			if(countsChar.getNumDifferent()>0) sb.append(countsChar.selectBest(1).iterator().next());
			else sb.append("N");
		}
		return sb.toString();
	}
	private Map<String, Integer> buildTwoClusters(List<String> alleles, String referenceAllele, String consensus) {
		int refCount = 0;
		int consensusCount = 0;
		for(String allele: alleles) {
			int diffRef = calculateHammingDistance(referenceAllele,allele);
			int diffCon = calculateHammingDistance(consensus,allele);
			if(diffCon>=diffRef) refCount++;
			else consensusCount++;
		}
		Map<String, Integer> answer = new TreeMap<>();
		answer.put(referenceAllele, refCount);
		answer.put(consensus, consensusCount);
		return answer;
	}
	private Map<String, Integer> buildClustersByHammingDistance(List<String> alleles, String referenceAllele) {
		// TODO Implement better
		Map<String, Integer> answer = new TreeMap<>();
		answer.put(referenceAllele, alleles.size());
		return answer;
	}
	private int calculateHammingDistance(String s1, String s2) {
		int answer = 0;
		int l = s1.length();
		for(int i=0;i<l;i++) {
			if(s1.charAt(i)!=s2.charAt(i)) answer++;
		}
		return answer;
	}
	
	
}
