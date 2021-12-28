package ngsep.genome;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;


public class HomologRelationshipsFinder {
	
	public static final byte DEF_KMER_LENGTH = 6;
	public static final int DEF_MIN_PCT_KMERS = 5;
	public static final int DEF_MIN_NUM_KMERS = 3;
	
	private Logger log = Logger.getAnonymousLogger();
	
	private byte kmerLength = DEF_KMER_LENGTH;
	private int minPctKmers = DEF_MIN_PCT_KMERS;
	
	
	
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
	
	
	
	public List<HomologyEdge> calculateParalogs (List<HomologyUnit> units) {
		Map<Long,Set<Integer>> unitsByKmer =indexKmersHomologyUnits(units,kmerLength);
		return calculateHomologs(units, units, unitsByKmer, kmerLength, minPctKmers);
	}
	
	public List<HomologyEdge> calculateOrthologs (List<HomologyUnit> queryUnits, List<HomologyUnit> subjectUnits) {
		Map<Long,Set<Integer>> subjectUnitsByKmer =indexKmersHomologyUnits(subjectUnits, kmerLength);
		return calculateHomologs(queryUnits, subjectUnits, subjectUnitsByKmer, kmerLength, minPctKmers);
	}
	
	public static Map<Long,Set<Integer>> indexKmersHomologyUnits(List<HomologyUnit> units, int kmerLength) {
		Map<Long,Set<Integer>> unitsByKmer = new HashMap<>();
		for(int i=0;i<units.size();i++) {
			HomologyUnit unit = units.get(i);
			List<Long> kmerCodes = unit.getKmerCodes(kmerLength, 1);
			for(Long code:kmerCodes) {
				Set<Integer> unitsKmer = unitsByKmer.computeIfAbsent(code, v->new HashSet<Integer>());
				unitsKmer.add(i);
			}
		}
		return unitsByKmer;
	}
	public static List<HomologyEdge> calculateHomologs( List<HomologyUnit> queryUnits, List<HomologyUnit> subjectUnits, Map<Long,Set<Integer>> subjectUnitsByKmer, int kmerLength, int minPctKmers) {
		List<HomologyEdge> edges = new ArrayList<>();
		int totalQueryUnits = queryUnits.size();
		if(totalQueryUnits>1000) System.out.println("Calculating edges of catalog with size "+totalQueryUnits+" first unit "+queryUnits.get(0).getId());
		int processed = 0;
		for(HomologyUnit unit1:queryUnits) {
			Set<Long> kmerCodes = new HashSet<Long>(unit1.getKmerCodes(kmerLength, 1));
			int n = kmerCodes.size();
			if(n==0) {
				System.out.println("Calculating edges of catalog with size "+totalQueryUnits+" Unit "+unit1.getUniqueKey()+" has zero kmers");
				continue;
			}
			int [] countsSubjectUnits = new int [subjectUnits.size()];
			Arrays.fill(countsSubjectUnits, 0);
			for(long code:kmerCodes) {
				Set<Integer> subjectUnitIdxsKmer = subjectUnitsByKmer.get(code);
				if(subjectUnitIdxsKmer==null) continue;
				for(int i:subjectUnitIdxsKmer) (countsSubjectUnits[i])++;
			}
			
			for(int i=0;i<countsSubjectUnits.length;i++) {
				int count = countsSubjectUnits[i];
				if(count<DEF_MIN_NUM_KMERS) continue;
				HomologyUnit unit2 = subjectUnits.get(i);
				if(unit1==unit2) continue;
				
				double score = 100.0*count/n;
				if(score <minPctKmers) continue;
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
