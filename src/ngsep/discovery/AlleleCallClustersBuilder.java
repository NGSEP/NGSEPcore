/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.discovery;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import ngsep.math.CountsRankHelper;
import ngsep.sequences.DNASequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.SequenceDistanceMeasure;

/**
 * Clustering of allele calls for indel and STR discovery and genotyping
 * @author Jorge Duitama
 *
 */
public class AlleleCallClustersBuilder {

	private String sequenceName;
	private int position;
	private Map<Integer, List<PileupAlleleCall>> alleleCallsByLength = new HashMap<>();
	private int totalAlleleCalls = 0;
	private SequenceDistanceMeasure distanceMeasure = new HammingSequenceDistanceMeasure();
	private double minRelativeProportion = 0.2;
	private int posPrint = -1;
	
	
	
	
	public AlleleCallClustersBuilder(String sequenceName, int position) {
		super();
		this.sequenceName = sequenceName;
		this.position = position;
	}
	
	
	/**
	 * @return the sequenceName
	 */
	public String getSequenceName() {
		return sequenceName;
	}

	/**
	 * @return the position
	 */
	public int getPosition() {
		return position;
	}


	public void addAlleleCall (PileupAlleleCall call) {
		String allele = call.getAlleleString();
		int l = allele.length();
		List<PileupAlleleCall> allelesLength = alleleCallsByLength.get(l);
		if(allelesLength==null) {
			allelesLength = new ArrayList<>();
			alleleCallsByLength.put(l, allelesLength);
		}
		allelesLength.add(call);
		totalAlleleCalls++;
	}
	public Map<String,List<PileupAlleleCall>> clusterAlleleCalls (String [] suggestedAlleles, boolean allowNewAlleles) {
		Map<String,List<PileupAlleleCall>> alleleClusters = new TreeMap<>();
		if(position == posPrint) System.out.println("Total length clusters: "+alleleCallsByLength.size());
		Map<Integer, List<PileupAlleleCall>> filteredClusters = filterLengthClusters();
		if(position == posPrint) System.out.println("Filtered clusters: "+filteredClusters.size());
		int s = filteredClusters.size();
		if(s==0) return alleleClusters;
		for(int l:filteredClusters.keySet()) {
			List<PileupAlleleCall> callsL = filteredClusters.get(l);	
			Set<String> suggestedAllelesSet = getAllelesByLength(suggestedAlleles,l);
			Map<String,List<PileupAlleleCall>> lengthClusters;
			if(allowNewAlleles) {
				if(position == posPrint) System.out.println("Suggested alleles set: "+suggestedAllelesSet+" calls: "+callsL.size()+" initial suggested: "+suggestedAlleles.length);
				if(callsL.size()<5*suggestedAllelesSet.size()) {
					//With low coverage and suggested alleles (probably the reference), only those are taken into account
					lengthClusters = clusterAlleleCallsPivotAlleles(callsL,suggestedAllelesSet);
				} else {
					//With enough calls, the consensus is considered as a suggested allele
					List<String> allelesL = new ArrayList<>(callsL.size());
					//List<String> qualities = new ArrayList<>(callsL.size());
					for(PileupAlleleCall call:callsL) {
						allelesL.add(call.getAlleleString());
						//qualities.add(call.getQualityScores());
					}
					String consensus = HammingSequenceDistanceMeasure.makeHammingConsensus(allelesL);
					if(position == posPrint) System.out.println("Consensus: "+consensus);
					suggestedAllelesSet.add(consensus);
					if(l<4 || suggestedAllelesSet.size()>1 || callsL.size()<10) {
						lengthClusters = clusterAlleleCallsPivotAlleles(callsL,suggestedAllelesSet);
					} else {
						//Only if large enough potential allele, enough read depth and no previously suggested alleles, try to break length cluster 
						lengthClusters = clusterAlleleCallsByVariantSites(callsL, consensus);
						if(position == posPrint) System.out.println("Clusters for length "+l+" "+lengthClusters.size());
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
	 * Filter severe imbalances in read count for clusters based on length. Useful to remove stutter products
	 * @return Map<Integer, List<PileupAlleleCall>> Length clusters passing the minimum proportion of reads
	 */
	private Map<Integer, List<PileupAlleleCall>> filterLengthClusters() {
		Map<Integer, List<PileupAlleleCall>> answer = new HashMap<>();
		//Only filter if more than 2 length clusters
		if(alleleCallsByLength.size()<3) return alleleCallsByLength;
		double minCount  = minRelativeProportion*totalAlleleCalls; 
		for(int l:alleleCallsByLength.keySet()) {
			List<PileupAlleleCall> allelesL = alleleCallsByLength.get(l);
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
	private Map<String, List<PileupAlleleCall>> clusterAlleleCallsPivotAlleles(List<PileupAlleleCall> calls, Set<String> suggestedAlleles) {
		assert suggestedAlleles.size() > 0;
		Map<String, List<PileupAlleleCall>> clusters = new TreeMap<>();
		for(PileupAlleleCall call:calls) {
			String alleleC = call.getAlleleString();
			//Find closest allele
			double minDistance = 1.1;
			String minAllele = null;
			for(String alleleS:suggestedAlleles) {
				double d = distanceMeasure.calculateNormalizedDistance(alleleC, alleleS);
				if(minAllele==null || minDistance>d ) {
					minAllele = alleleS;
					minDistance = d;
				}
			}
			if(minDistance<=0.5) {
				//If not too far from the best allele, add to cluster
				List<PileupAlleleCall> cluster = clusters.get(minAllele);
				if(cluster==null) {
					cluster = new ArrayList<>();
					clusters.put(minAllele, cluster);
				}
				cluster.add(call);
			}
		}	
		return clusters;
	}
	private Map<String, List<PileupAlleleCall>> clusterAlleleCallsByVariantSites(List<PileupAlleleCall> calls, String consensus) {
		Map<String,List<PileupAlleleCall>> answer = new TreeMap<>();
		if(calls.size()==0) return answer;
		
		double [] heterozygousPosteriors = calculateHetPosteriors(calls, consensus);
		
		List<Integer> variantSites = new ArrayList<>();
		for(int i=0;i<heterozygousPosteriors.length;i++) {
			// TODO: Use parameter
			if(heterozygousPosteriors[i]>=0.9) {
				variantSites.add(i);
				if(posPrint==position) System.out.println("Position: "+i+" heterozygous posterior: "+heterozygousPosteriors[i]);
			}
		}
		if(variantSites.size()==0) {
			//No real variation
			answer.put(consensus, calls);
			return answer;
		}
		return clusterAlleleCallsByHaplotypes (calls, variantSites);
	}
	
	private double[] calculateHetPosteriors(List<PileupAlleleCall> calls, String consensus) {
		int l = consensus.length();
		double[] answer = new double[l];
		for(int i=0;i<l;i++) {
			char c = consensus.charAt(i);
			int idxC = DNASequence.BASES_STRING.indexOf(c);
			//Check first if variable
			boolean variable = false;
			for(int j=0;j<calls.size() && !variable;j++) {
				String allele = calls.get(j).getAlleleString();
				variable = allele.charAt(i)!=c;
			}
			if(!variable) {
				answer[i] =0;
				continue;
			}
			CountsHelper helper = new CountsHelper();
			for(int j=0;j<calls.size();j++) {
				PileupAlleleCall call = calls.get(j);
				String allele = call.getAlleleString();
				String qualities = call.getQualityScores();
				//TODO: replace constant with parameter
				byte q = (byte)(Math.min(VariantPileupListener.DEF_MAX_BASE_QS, qualities.charAt(i)-33));
				helper.updateCounts(allele.substring(i, i+1), q, false);
			}
			double [][] posteriors = helper.getPosteriorProbabilities(VariantPileupListener.DEF_HETEROZYGOSITY_RATE_DIPLOID);
			answer[i] = 0;
			for(int k=0;k<DNASequence.BASES_ARRAY.length;k++) {
				double hetPost = posteriors[idxC][k]+posteriors[k][idxC]; 
				if(k!=idxC && hetPost>answer[i]) answer[i] = hetPost;
			}
		}
		return answer;
	}


	private Map<String, List<PileupAlleleCall>> clusterAlleleCallsByHaplotypes(List<PileupAlleleCall> calls, List<Integer> variantSites) {
		
		int n = calls.size();
		int m = variantSites.size();
		
		//Build haplotypes
		char [][] haplotypes = new char[n][m];
		//byte [][] qualities = new byte[n][m];
		for(int i=0;i<n;i++) {
			PileupAlleleCall call = calls.get(i);
			String seq = call.getAlleleString();
			//String scores = call.getQualityScores();
			for(int j=0;j<m;j++) {
				int k = variantSites.get(j);
				haplotypes[i][j]=seq.charAt(k);
				//TODO: replace constant with parameter
				//qualities[i][j] = (byte)(Math.min(VariantPileupListener.DEF_MAX_BASE_QS, scores.charAt(k)-33));
			}
		}
		
		
		//Cluster first sequences having the same haplotype on variant sites
		CountsRankHelper<String> hapCounts = new CountsRankHelper<>();
		for(int i=0;i<haplotypes.length;i++) {
			hapCounts.add(new String (haplotypes[i]));
		}
		int maxHaps = 2;
		if(m>3) maxHaps = m/2+1; 
		List<String> selectedHaplotypes = new ArrayList<>(hapCounts.selectBest(maxHaps).keySet());
		//Calculate hamming consensus for the chosen haplotypes
		List<List<String>> selectedSequences = new ArrayList<>();
		for(int i=0;i<selectedHaplotypes.size();i++) selectedSequences.add(new ArrayList<>());
		for(int i=0;i<haplotypes.length;i++) {
			PileupAlleleCall call = calls.get(i);
			String hap = new String (haplotypes[i]);
			int j = selectedHaplotypes.indexOf(hap);
			if(j>=0) selectedSequences.get(j).add(call.getAlleleString());
		}
		
		Set<String> consensus = new TreeSet<>();
		for(List<String> seqs:selectedSequences) {
			if(seqs.size()>0) consensus.add(HammingSequenceDistanceMeasure.makeHammingConsensus(seqs));
		}
		return clusterAlleleCallsPivotAlleles(calls, consensus);
	}
	
}
