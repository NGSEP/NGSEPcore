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
package ngsep.genome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import JSci.maths.statistics.ChiSqrDistribution;
import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.alignments.PairwiseAlignment;
import ngsep.math.Distribution;
import ngsep.math.ShannonEntropyCalculator;
import ngsep.sequences.AbstractLimitedSequence;
import ngsep.sequences.AminoacidSequence;
import ngsep.sequences.KmersExtractor;

/**
 * 
 * @author Jorge Gomez
 *
 */
public class HomologRelationshipsFinder {
	
	public static final byte DEF_KMER_LENGTH = 6;
	public static final int DEF_MIN_PCT_KMERS = 11;
	public static final int DEF_MIN_WEIGHTED_COUNT = 5;
	public static final int DEF_MAX_LENGTH_RATIO = 2;
	
	private Logger log = Logger.getAnonymousLogger();
	
	private byte kmerLength = DEF_KMER_LENGTH;
	private int minPctKmers = DEF_MIN_PCT_KMERS;
	private int maxLengthRatio = DEF_MAX_LENGTH_RATIO;
	
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	public byte getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(byte kmerLength) {
		this.kmerLength = kmerLength;
	}
	public int getMinPctKmers() {
		return minPctKmers;
	}
	public void setMinPctKmers(int minPctKmers) {
		this.minPctKmers = minPctKmers;
	}
	public int getMaxLengthRatio() {
		return maxLengthRatio;
	}
	public void setMaxLengthRatio(int maxLengthRatio) {
		this.maxLengthRatio = maxLengthRatio;
	}
	
	public List<HomologyEdge> calculateParalogs (List<HomologyUnit> units) {
		Map<Long,Set<Integer>> unitsByKmer =indexKmersHomologyUnits(units);
		System.out.println("Indexed: "+unitsByKmer.size()+" kmers");
		return calculateHomologs(units, units, unitsByKmer);
	}
	
	public List<HomologyEdge> calculateOrthologs (List<HomologyUnit> queryUnits, List<HomologyUnit> subjectUnits) {
		Map<Long,Set<Integer>> subjectUnitsByKmer =indexKmersHomologyUnits(subjectUnits);
		System.out.println("Indexed: "+subjectUnitsByKmer.size()+" subject kmers");
		return calculateHomologs(queryUnits, subjectUnits, subjectUnitsByKmer);
	}
	
	public Map<Long,Set<Integer>> indexKmersHomologyUnits(List<HomologyUnit> units) {
		Map<Long,Set<Integer>> unitsByKmer = new HashMap<>();
		for(int i=0;i<units.size();i++) {
			HomologyUnit unit = units.get(i);
			Map<Long,Double> kmerCodes = unit.getKmerCodesWithEntropies(kmerLength, 1);
			//if("TcDm25H1_000257900".equals(unit.getId())) System.out.println("Unit: "+unit.getId()+" Kmer codes: "+kmerCodes.size());
			for(Long code:kmerCodes.keySet()) {
				Set<Integer> unitsKmer = unitsByKmer.computeIfAbsent(code, v->new HashSet<Integer>());
				unitsKmer.add(i);
			}
		}
		return unitsByKmer;
	}
	public List<HomologyEdge> calculateHomologs( List<HomologyUnit> queryUnits, List<HomologyUnit> subjectUnits, Map<Long,Set<Integer>> subjectUnitsByKmer) {
		Distribution kmerAbundances = new Distribution(1, 100, 1);
		//System.out.println("Kmer abundances distribution");
		for(Map.Entry<Long,Set<Integer>> entry:subjectUnitsByKmer.entrySet()) {
			int n = entry.getValue().size();
			kmerAbundances.processDatapoint(n);
			//if(n>50) System.out.println("Kmer: "+(new String(AbstractLimitedSequence.getSequence(entry.getKey(), kmerLength, new AminoacidSequence())))+" count: "+n);
		}
		//kmerAbundances.printDistributionInt(System.out);
		double average = kmerAbundances.getAverage();
		double mode = kmerAbundances.getLocalMode(1, 10);
		log.info("Kmer abundances distribution. Average: "+average+" Local mode: "+mode);
		ChiSqrDistribution chiSq = new ChiSqrDistribution((mode+average)/2);
		List<HomologyEdge> edges = new ArrayList<>();
		int totalQueryUnits = queryUnits.size();
		if(totalQueryUnits>1000) log.info("Calculating edges of catalog with size "+totalQueryUnits+" first unit "+queryUnits.get(0).getId());
		int processed = 0;
		for(HomologyUnit unit1:queryUnits) {
			boolean debug = subjectUnits.size()<3;
			//debug = debug || unit1.getId().equals("TcDm25H1_000257900");
			Map<Long,Double> kmerCodesWithEntropies = unit1.getKmerCodesWithEntropies(kmerLength, 1);
			int length1 = unit1.getUnitSequence().length();
			int n = kmerCodesWithEntropies.size();
			if(n==0) {
				System.out.println("Calculating edges of catalog with size "+totalQueryUnits+" Unit "+unit1.getUniqueKey()+" has zero kmers");
				continue;
			}
			Map<Integer,Double> rawCountsSubjectUnits = new HashMap<>();
			Map<Integer,Double> weightedCountsSubjectUnits = new HashMap<>();
			double totalWeight = 0;
			for(Map.Entry<Long,Double> entry:kmerCodesWithEntropies.entrySet()) {
				long code = entry.getKey();
				Set<Integer> subjectUnitIdxsKmer = subjectUnitsByKmer.get(code);
				//A kmer must add to the total weight even if it does not appear in the subject kmers
				double pValueCode = 1;
				if(subjectUnitIdxsKmer!=null) {
					//pValueCode = kmerAbundances.getEmpiricalPvalue(subjectUnitIdxsKmer.size());
					pValueCode = 1-chiSq.cumulative(subjectUnitIdxsKmer.size());
				}
				pValueCode = Math.max(pValueCode, 0.00001);
				double score = Math.max(1,-Math.log10(pValueCode));
				double weight = Math.max(0.5,entry.getValue())/score;
				totalWeight+=weight;
				if(subjectUnitIdxsKmer==null) continue;
				if(debug) {
					String kmer = new String(AbstractLimitedSequence.getSequence(code, kmerLength, new AminoacidSequence()));
					int limit = (queryUnits==subjectUnits)?1:0;
					if(subjectUnitIdxsKmer.size()>limit) System.out.println("next kmer: "+kmer+" entropy: "+ShannonEntropyCalculator.calculateEntropy(kmer)+" hits: "+subjectUnitIdxsKmer+" pvalue: "+pValueCode+" score: "+score+" weight: "+weight);
				}
				
				for(int i:subjectUnitIdxsKmer) {
					rawCountsSubjectUnits.compute(i, (k,v)->(v==null)?1:v+1);
					weightedCountsSubjectUnits.compute(i, (k,v)->(v==null)?weight:v+weight);
				}
			}
			
			for(Map.Entry<Integer, Double> entry:weightedCountsSubjectUnits.entrySet()) {
				int i = entry.getKey();
				double weightedCount = entry.getValue();
				double score = 100.0*weightedCount/totalWeight;
				double rawCount = rawCountsSubjectUnits.get(i);
				if(debug) {
					double rawPCT = 100.0*rawCount/n;
					System.out.println("next subject: "+subjectUnits.get(i).getId()+" weighted count: "+weightedCount+" total weight: "+totalWeight+" score: "+score+" raw count: "+rawCount+" raw score: "+rawPCT);
				}
				if(weightedCount<DEF_MIN_WEIGHTED_COUNT) continue;
				HomologyUnit unit2 = subjectUnits.get(i);
				if(unit1==unit2) continue;
				if(score <minPctKmers) continue;
				int length2 = unit2.getUnitSequence().length();
				if(length2>maxLengthRatio*length1) continue; 
				if((weightedCount<10 && score < 10) && !passSecondaryFilters(unit1,unit2,kmerLength,debug)) continue;
				HomologyEdge edge = new HomologyEdge(unit1, unit2, score);
				unit1.addHomologRelationship(edge);
				edges.add(edge);
				if(debug) System.out.println("Added relationship between "+unit1.getId()+" and "+unit2.getId());
			}
			processed++;
			if(totalQueryUnits>1000 && processed%1000==0) log.info("Processed "+processed+" of "+totalQueryUnits+" units. last unit "+unit1.getId());
		}
		if(totalQueryUnits>1000) log.info("Calculated edges of catalog with size "+totalQueryUnits+" first unit "+queryUnits.get(0).getId());
		return edges;
	}
	private static boolean passSecondaryFilters(HomologyUnit unit1, HomologyUnit unit2, int kmerLength, boolean debug) {
		PairwiseAlignerSimpleGap aligner = new PairwiseAlignerSimpleGap();
		aligner.setLocal(true);
		//String [] alignment = aligner.calculateAlignment(unit1.getUnitSequence(), unit2.getUnitSequence());
		String seq1 = unit1.getUnitSequence().toString();
		String seq2 = unit2.getUnitSequence().toString();
		String [] kmers1= KmersExtractor.extractKmers(seq1, kmerLength, 1, 0, seq1.length(), true, true, false);
		String [] kmers2 = KmersExtractor.extractKmers(seq2, kmerLength, 1, 0, seq2.length(), true, true, false);
		Map<String,Integer> lastOccKmerSites2 = new HashMap<>();
		for(int i=0;i<kmers2.length;i++) {
			lastOccKmerSites2.put(kmers2[i],i);
		}
		int first1=seq1.length();
		int last1=0;
		int first2=seq2.length();
		int last2=0;
		for(int i=0;i<kmers1.length;i++) {
			Integer j = lastOccKmerSites2.get(kmers1[i]);
			if(j!=null) {
				first1 = Math.min(first1, i);
				last1 = Math.max(last1, i+kmerLength);
				first2 = Math.min(first2, j);
				last2 = Math.max(last2, j+kmerLength);		
			}
		}
		if(first1>=last1 || first2>=last2) return false;
		PairwiseAlignment alignment = aligner.calculateAlignment(seq1.substring(first1,last1), seq2.substring(first2,last2));
		if(debug) System.out.println("Aln: "+alignment.getAlignedSequence1()+" "+alignment.getAlignedSequence2()+" score: "+alignment.getScore());
		
		return alignment.getScore()>=2*kmerLength;
	}
	
	//Old FMindex algorithm
	
	/*public List<HomologyEdge> calculateParalogs(AnnotatedReferenceGenome genome) {
	List<HomologyEdge> edges = new ArrayList<HomologyEdge>();
	List<HomologyUnit> units = genome.getHomologyUnits();
	for (HomologyUnit unit:units) {
		List<HomologyEdge> hits = findHomologs(unit, genome.getHomologyCatalog());
		edges.addAll(hits);
	}
	genome.selectUniqueOrthologyUnits();
	return edges;
	}
	
	public List<HomologyEdge> calculateParalogsOrganism(HomologyCatalog catalog) {
		List<HomologyEdge> edges = new ArrayList<HomologyEdge>();
		List<HomologyUnit> units = catalog.getHomologyUnits();
		for (HomologyUnit unit:units) {
			List<HomologyEdge> hits = findHomologs(unit, catalog);
			edges.addAll(hits);
		}
		return edges;
	}*/
	
	/**
	 * Finds orthologs of the orthology units in this genome in the given genome
	 * @param genome2 to search for orthologs
	 */
	/*public List<HomologyEdge> calculateOrthologs(HomologyCatalog catalog1, HomologyCatalog catalog2) {
		List<HomologyEdge> edges = new ArrayList<HomologyEdge>();
		List<HomologyUnit> units = catalog1.getHomologyUnits();
		for (HomologyUnit unit:units) {
			List<HomologyEdge> hits = findHomologs(unit, catalog2);
			edges.addAll(hits);
		}
		return edges;
	}*/
	
	/*private List<HomologyEdge> findHomologs(HomologyUnit unit, HomologyCatalog catalog) {
		List<HomologyEdge> edges = new ArrayList<HomologyEdge>();
		FMIndex indexCatalog = catalog.getIndexHomologyUnits();
		//Counts of k-mers mapping to each protein in the FM-index
		Map<String,Integer> kmerSupportMap = new TreeMap<>();
		int totalKmers = 0;
		List<String> kmers = unit.getKmers(kmerLength,kmerLength);
		//Step 1: Generate k-mers to query the FM-Index looking for homologous transcripts to calculate the kmer counts
		for(String kmer:kmers) {
			List <UngappedSearchHit> kmerHits = indexCatalog.exactSearch(kmer);
			for(UngappedSearchHit hit:kmerHits) {
				String name = hit.getSequenceName();
				kmerSupportMap.compute(name, (k,v)->v!=null?v+1:1);
			}
			totalKmers++;
		}

		//Step 2: Fill list traversing the counts and choosing transcripts for which at least x% of the k-mers support the match
		for(Map.Entry<String,Integer> entry : kmerSupportMap.entrySet()) {
			String homologId = entry.getKey();
			double transcriptKmers = entry.getValue();
			if(transcriptKmers<1+DEF_MIN_NUM_KMERS/2) continue;
			double percent = (transcriptKmers/totalKmers)*100;
			if(percent < minPctKmers) continue;
			HomologyUnit homolog = catalog.getHomologyUnit(homologId);
			if(homolog==null) continue;
			if(homolog==unit) continue;
			// TODO: calculate score
			double score = percent;
			HomologyEdge edge = new HomologyEdge(unit, homolog, score);
			unit.addHomologRelationship(edge);
			edges.add(edge);
		}
		return edges;
	}*/

}
