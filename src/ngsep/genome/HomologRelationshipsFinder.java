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
import ngsep.math.Distribution;
import ngsep.math.PhredScoreHelper;
import ngsep.math.ShannonEntropyCalculator;
import ngsep.sequences.AbstractLimitedSequence;
import ngsep.sequences.AminoacidSequence;
import ngsep.sequences.KmersExtractor;


public class HomologRelationshipsFinder {
	
	public static final byte DEF_KMER_LENGTH = 6;
	public static final int DEF_MIN_PCT_KMERS = 5;
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
		return calculateHomologs(units, units, unitsByKmer);
	}
	
	public List<HomologyEdge> calculateOrthologs (List<HomologyUnit> queryUnits, List<HomologyUnit> subjectUnits) {
		Map<Long,Set<Integer>> subjectUnitsByKmer =indexKmersHomologyUnits(subjectUnits);
		return calculateHomologs(queryUnits, subjectUnits, subjectUnitsByKmer);
	}
	
	public Map<Long,Set<Integer>> indexKmersHomologyUnits(List<HomologyUnit> units) {
		Map<Long,Set<Integer>> unitsByKmer = new HashMap<>();
		for(int i=0;i<units.size();i++) {
			HomologyUnit unit = units.get(i);
			Map<Long,Double> kmerCodes = unit.getKmerCodesWithEntropies(kmerLength, 1);
			for(Long code:kmerCodes.keySet()) {
				Set<Integer> unitsKmer = unitsByKmer.computeIfAbsent(code, v->new HashSet<Integer>());
				unitsKmer.add(i);
			}
		}
		return unitsByKmer;
	}
	public List<HomologyEdge> calculateHomologs( List<HomologyUnit> queryUnits, List<HomologyUnit> subjectUnits, Map<Long,Set<Integer>> subjectUnitsByKmer) {
		Distribution kmerAbundances = new Distribution(1, 100, 1);
		System.out.println("Kmer abundances distribution");
		for(Map.Entry<Long,Set<Integer>> entry:subjectUnitsByKmer.entrySet()) {
			int n = entry.getValue().size();
			kmerAbundances.processDatapoint(n);
			//if(n>50) System.out.println("Kmer: "+(new String(AbstractLimitedSequence.getSequence(entry.getKey(), kmerLength, new AminoacidSequence())))+" count: "+n);
		}
		kmerAbundances.printDistributionInt(System.out);
		double mode = kmerAbundances.getLocalMode(1, 10);
		System.out.println("Local mode: "+mode);
		ChiSqrDistribution chiSq = new ChiSqrDistribution(2*mode);
		List<HomologyEdge> edges = new ArrayList<>();
		int totalQueryUnits = queryUnits.size();
		if(totalQueryUnits>1000) System.out.println("Calculating edges of catalog with size "+totalQueryUnits+" first unit "+queryUnits.get(0).getId());
		int processed = 0;
		for(HomologyUnit unit1:queryUnits) {
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
				if(subjectUnitIdxsKmer==null) continue;
				//double pValueCode = kmerAbundances.getEmpiricalPvalue(subjectUnitIdxsKmer.size());
				double pValueCode = 1-chiSq.cumulative(subjectUnitIdxsKmer.size());
				double score = Math.max(1,-Math.log10(pValueCode));
				double weight = entry.getValue()/score+0.01;
				totalWeight+=weight;
				if(subjectUnits.size()<3) {
					String kmer = new String(AbstractLimitedSequence.getSequence(code, kmerLength, new AminoacidSequence()));
					if(subjectUnitIdxsKmer.size()>1) System.out.println("next kmer: "+kmer+" entropy: "+ShannonEntropyCalculator.calculateEntropy(kmer)+" hits: "+subjectUnitIdxsKmer+" weight: "+weight);
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
				if(subjectUnits.size()<3) {
					double rawPCT = 100.0*rawCount/n;
					System.out.println("next subject: "+subjectUnits.get(i).getId()+" weighted count: "+weightedCount+" total weight: "+totalWeight+" score: "+score+" raw count: "+rawCount+" raw score: "+rawPCT);
				}
				if(weightedCount<DEF_MIN_WEIGHTED_COUNT) continue;
				HomologyUnit unit2 = subjectUnits.get(i);
				if(unit1==unit2) continue;
				if(score <minPctKmers) continue;
				int length2 = unit2.getUnitSequence().length();
				if(length2>maxLengthRatio*length1) continue; 
				if((weightedCount<10 && score < 10) && !passSecondaryFilters(unit1,unit2,kmerLength)) continue;
				HomologyEdge edge = new HomologyEdge(unit1, unit2, score);
				unit1.addHomologRelationship(edge);
				edges.add(edge);
			}
			processed++;
			if(totalQueryUnits>1000 && processed%1000==0) System.out.println("Processed "+processed+" of "+totalQueryUnits+" units. last unit "+unit1.getId());
		}
		if(totalQueryUnits>1000) System.out.println("Calculated edges of catalog with size "+totalQueryUnits+" first unit "+queryUnits.get(0).getId());
		return edges;
	}
	private static boolean passSecondaryFilters(HomologyUnit unit1, HomologyUnit unit2, int kmerLength) {
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
		String [] alignment = aligner.calculateAlignment(seq1.substring(first1,last1), seq2.substring(first2,last2));
		
		//aligner.printAlignmentMatrix(aligner.getMatchScores(), unit1.getUnitSequence().toString(), unit2.getUnitSequence().toString());
		//System.out.println("Aln: "+alignment[0]+" "+alignment[1]+" score: "+aligner.getMaxScore());
		return aligner.getMaxScore()>=2*kmerLength;
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
