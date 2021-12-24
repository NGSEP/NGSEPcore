package ngsep.genome;

import java.util.ArrayList;
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
		Map<String,HomologyUnit> unitsByKey = new HashMap<String, HomologyUnit>();
		
		for(HomologyUnit unit:units) unitsByKey.put(unit.getUniqueKey(), unit);
		Map<String, Set<String>> kmersByUnit = calculateKmers(units);
		Map<String,Set<String>> unitsByKmer =indexKmersHomologyUnits(kmersByUnit);
		return calculateHomologs(units, kmersByUnit, unitsByKey, unitsByKmer);
	}
	
	public List<HomologyEdge> calculateOrthologs (List<HomologyUnit> queryUnits, List<HomologyUnit> subjectUnits) {
		Map<String,HomologyUnit> subjectUnitsByKey = new HashMap<String, HomologyUnit>();
		for(HomologyUnit unit:subjectUnits) subjectUnitsByKey.put(unit.getUniqueKey(), unit);
		Map<String,Set<String>> kmersByQueryUnit = calculateKmers(queryUnits);
		Map<String,Set<String>> kmersBySubjectUnit = calculateKmers(subjectUnits);
		Map<String,Set<String>> subjectUnitsByKmer =indexKmersHomologyUnits(kmersBySubjectUnit);
		return calculateHomologs(queryUnits, kmersByQueryUnit, subjectUnitsByKey, subjectUnitsByKmer);
	}
	
	private Map<String, Set<String>> calculateKmers(List<HomologyUnit> units) {
		Map<String,Set<String>> kmersByUnit = new HashMap<String, Set<String>>();
		for(HomologyUnit unit:units) {
			List<String> kmers = unit.getKmers(kmerLength,1);
			Set<String> kmersSet = new HashSet<String>(kmers);
			kmersByUnit.put(unit.getUniqueKey(), kmersSet);
		}
		return kmersByUnit;
	}
	private Map<String, Set<String>> indexKmersHomologyUnits(Map<String, Set<String>> kmersByUnit) {
		Map<String, Set<String>> unitsByKmer = new HashMap<>();
		for(Map.Entry<String,Set<String>> entry:kmersByUnit.entrySet()) {
			String uniqueKey = entry.getKey();
			Set<String> kmersSet = entry.getValue();
			for(String kmer:kmersSet) {
				Set<String> unitsKmer = unitsByKmer.computeIfAbsent(kmer, v->new HashSet<String>());
				unitsKmer.add(uniqueKey);
			}
		}
		return unitsByKmer;
	}
	private List<HomologyEdge> calculateHomologs( List<HomologyUnit> queryUnits, Map<String, Set<String>> kmersByQueryUnit, Map<String, HomologyUnit> subjectUnitsByKey, Map<String, Set<String>> subjectUnitsByKmer) {
		List<HomologyEdge> edges = new ArrayList<>();
		int totalQueryUnits = queryUnits.size();
		if(totalQueryUnits>1000) log.info("Calculating edges of catalog with size "+totalQueryUnits+" first unit "+queryUnits.get(0).getId());
		int processed = 0;
		for(HomologyUnit unit1:queryUnits) {
			//TODO: Parameters
			Set<String> kmers1 = kmersByQueryUnit.get(unit1.getUniqueKey());
			int n = kmers1.size();
			if(n==0) {
				log.warning("Calculating edges of catalog with size "+totalQueryUnits+" Unit "+unit1.getUniqueKey()+" has zero kmers");
				continue;
			}
			Map<String,Integer> countsSubjectUnits = new HashMap<String, Integer>();
			for(String kmer:kmers1) {
				Set<String> subjectUnitsKmer = subjectUnitsByKmer.get(kmer);
				if(subjectUnitsKmer==null) continue;
				for(String subjectUnitKey:subjectUnitsKmer) {
					countsSubjectUnits.compute(subjectUnitKey, (k,v)->v!=null?v+1:1);
				}
			}
			
			for(Map.Entry<String,Integer> entry:countsSubjectUnits.entrySet()) {
				String subjectKey = entry.getKey();
				int count = entry.getValue();
				if(count<DEF_MIN_NUM_KMERS) continue;
				HomologyUnit unit2 = subjectUnitsByKey.get(subjectKey);
				if(unit1==unit2) continue;
				
				double score = 100.0*count/n;
				if(score <minPctKmers) continue;
				HomologyEdge edge = new HomologyEdge(unit1, unit2, score);
				unit1.addHomologRelationship(edge);
				edges.add(edge);
			}
			processed++;
			if(totalQueryUnits>1000 && processed%1000==0) log.info("Processed "+processed+" of "+totalQueryUnits+" units. last unit "+unit1.getId());
		}
		if(totalQueryUnits>1000) log.info("Calculated edges of catalog with size "+totalQueryUnits+" first unit "+queryUnits.get(0).getId());
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
