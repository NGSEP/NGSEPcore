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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import ngsep.math.CountsRankHelper;
import ngsep.sequences.DNASequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;

/**
 * Clustering of allele calls for indel and STR discovery and genotyping
 * @author Jorge Duitama
 *
 */
public class AlleleCallClustersBuilder {

	private static final double DEF_MIN_RELATIVE_PROPORTION = 0.2;
	private static final double DEF_MIN_HET_POSTERIOR = 0.99;
	private String sequenceName;
	private int position;

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

	public String [] clusterAlleleCalls (PileupRecord pileup, List<PileupAlleleCall> calls, String reference, byte maxBaseQS) {
		Set<String> allelesSet = new TreeSet<>();
		Map<Integer, List<PileupAlleleCall>> alleleCallsByLength = new HashMap<>();
		for(PileupAlleleCall call:calls) {
			String allele = call.getAlleleString();
			int l = allele.length();
			List<PileupAlleleCall> allelesLength = alleleCallsByLength.get(l);
			if(allelesLength==null) {
				allelesLength = new ArrayList<>();
				alleleCallsByLength.put(l, allelesLength);
			}
			allelesLength.add(call);
		}
		
		if(position == posPrint) System.out.println("Loaded calls: "+calls.size()+" Total length clusters: "+alleleCallsByLength.size());
		Map<Integer, List<PileupAlleleCall>> filteredClusters = filterLengthClusters(alleleCallsByLength, calls.size());
		if(position == posPrint) System.out.println("Filtered clusters: "+filteredClusters.size());
		
		for(int l:filteredClusters.keySet()) {
			List<PileupAlleleCall> callsL = filteredClusters.get(l);	
			Set<String> suggestedAllelesSet = new HashSet<>();
			if(l==reference.length()) suggestedAllelesSet.add(reference);
			Set<String> lengthAlleles;
			if(position == posPrint) System.out.println("Suggested alleles set: "+suggestedAllelesSet+" calls: "+callsL.size());
			if(callsL.size()<5*suggestedAllelesSet.size()) {
				//With low coverage and suggested alleles (probably the reference), only those are taken into account
				lengthAlleles = suggestedAllelesSet;
				//lengthClusters = clusterAlleleCallsPivotAlleles(callsL,suggestedAllelesSet);
			} else {
				//With enough calls, the consensus is considered as a suggested allele
				List<String> allelesL = new ArrayList<>(callsL.size());
				//List<String> qualities = new ArrayList<>(callsL.size());
				for(PileupAlleleCall call:callsL) {
					allelesL.add(call.getAlleleString());
					//qualities.add(call.getQualityScores());
				}
				String consensus = HammingSequenceDistanceMeasure.makeHammingConsensus(allelesL);
				if(!DNASequence.isDNA(consensus)) {
					System.err.println ("Consensus allele for calls "+allelesL+" of length: "+l+" at "+pileup.getSequenceName()+":"+pileup.getPosition()+" is not DNA: "+consensus  );
				}
				if(position == posPrint) System.out.println("Consensus: "+consensus);
				suggestedAllelesSet.add(consensus);
				if(l<4 || suggestedAllelesSet.size()>1 || callsL.size()<10) {
					lengthAlleles = suggestedAllelesSet;
					//lengthClusters = clusterAlleleCallsPivotAlleles(callsL,suggestedAllelesSet);
				} else {
					//Only if large enough potential allele, enough read depth and no previously suggested alleles, try to break length cluster
					lengthAlleles = splitAllelesByVariantSites(callsL, consensus, maxBaseQS);
					if(position == posPrint) System.out.println("Alleles for length "+l+" "+lengthAlleles.size());
				}
				
			}
			allelesSet.addAll(lengthAlleles);
		}
		allelesSet.add(reference);
		if(allelesSet.size()>100) System.err.println("WARN: Number of alleles for site at "+pileup.getSequenceName()+":"+pileup.getPosition()+" is "+allelesSet.size()+" ref Allele: "+reference);
		String [] answer = new String [allelesSet.size()];
		answer[0] = reference;
		//if(pileup.getPosition()==posPrint) System.out.println("Reference allele for indel: "+referenceAllele);
		int i=1;
		for (String allele:allelesSet) {
			if(!allele.equals(reference)) {
				answer[i] = allele;
				//if(pileup.getPosition()==posPrint) System.out.println("Next alternative allele for indel: "+allele);
				i++;
			}
		}
		return answer;
	}

	/**
	 * Filter severe imbalances in read count for clusters based on length. Useful to remove stutter products
	 * @return Map<Integer, List<PileupAlleleCall>> Length clusters passing the minimum proportion of reads
	 */
	private Map<Integer, List<PileupAlleleCall>> filterLengthClusters(Map<Integer, List<PileupAlleleCall>> alleleCallsByLength, int totalAlleleCalls) {
		Map<Integer, List<PileupAlleleCall>> answer = new HashMap<>();
		//Only filter if more than 2 length clusters
		if(alleleCallsByLength.size()<3) return alleleCallsByLength;
		double minCount  = DEF_MIN_RELATIVE_PROPORTION*totalAlleleCalls; 
		for(int l:alleleCallsByLength.keySet()) {
			List<PileupAlleleCall> allelesL = alleleCallsByLength.get(l);
			if(minCount<=allelesL.size()) answer.put(l, allelesL);
		}
		return answer;
	}
	
	/**
	 * 
	 * @param callsL. Not empty
	 * @param consensus
	 * @return
	 */
	private Set<String> splitAllelesByVariantSites(List<PileupAlleleCall> calls, String consensus, byte maxBaseQS) {
		Set<String> answer = new TreeSet<>();
		
		double [] heterozygousPosteriors = calculateHetPosteriors(calls, consensus, maxBaseQS);
		
		List<Integer> variantSites = new ArrayList<>();
		for(int i=0;i<heterozygousPosteriors.length;i++) {
			// TODO: Use parameter
			if(heterozygousPosteriors[i]>=DEF_MIN_HET_POSTERIOR) {
				variantSites.add(i);
				if(posPrint==position) System.out.println("Position: "+i+" heterozygous posterior: "+heterozygousPosteriors[i]);
			}
		}
		if(variantSites.size()==0) {
			//No real variation
			answer.add(consensus);
			return answer;
		}
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
				//qualities[i][j] = (byte)(Math.min(maxBaseQS, scores.charAt(k)-33));
			}
		}
		
		
		//Cluster first sequences having the same haplotype on variant sites
		CountsRankHelper<String> hapCounts = new CountsRankHelper<>();
		for(int i=0;i<haplotypes.length;i++) {
			hapCounts.add(new String (haplotypes[i]));
		}
		int maxHaps = 2;
		if(m>3) maxHaps = Math.min(10, m/2+1); 
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
		for(List<String> seqs:selectedSequences) {
			if(seqs.size()>0) answer.add(HammingSequenceDistanceMeasure.makeHammingConsensus(seqs));
		}
		return answer;
	}
	
	private double[] calculateHetPosteriors(List<PileupAlleleCall> calls, String consensus, byte maxbaseQS) {
		int l = consensus.length();
		double[] answer = new double[l];
		for(int i=0;i<l;i++) {
			char c = consensus.charAt(i);
			int idxC = DNASequence.BASES_STRING.indexOf(c);
			if(idxC<0) {
				answer[i]=0;
				continue;
			}
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
				byte q = (byte)(Math.min(maxbaseQS, qualities.charAt(i)-33));
				helper.updateCounts(allele.substring(i, i+1), q, false);
			}
			double [][] posteriors = helper.getPosteriorProbabilities(CountsHelper.DEF_HETEROZYGOSITY_RATE_DIPLOID);
			answer[i] = 0;
			for(int k=0;k<DNASequence.BASES_ARRAY.length;k++) {
				double hetPost = posteriors[idxC][k]+posteriors[k][idxC]; 
				if(k!=idxC && hetPost>answer[i]) answer[i] = hetPost;
			}
		}
		return answer;
	}
	
}
