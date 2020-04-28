package ngsep.genome;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.sequences.FMIndex;
import ngsep.sequences.UngappedSearchHit;

public class HomologRelationshipsFinder {
	public static final byte DEF_KMER_LENGTH = 10;
	public static final int DEF_MIN_PCT_KMERS = 50;
	
	private byte kmerLength = DEF_KMER_LENGTH;
	private int minPctKmers = DEF_MIN_PCT_KMERS;
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
	
	public List<HomologyEdge> calculateParalogs(AnnotatedReferenceGenome genome) {
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
	}
	
	/**
	 * Finds orthologs of the orthology units in this genome in the given genome
	 * @param genome2 to search for orthologs
	 */
	public List<HomologyEdge> calculateOrthologs(HomologyCatalog catalog1, HomologyCatalog catalog2) {
		List<HomologyEdge> edges = new ArrayList<HomologyEdge>();
		List<HomologyUnit> units = catalog1.getHomologyUnits();
		for (HomologyUnit unit:units) {
			List<HomologyEdge> hits = findHomologs(unit, catalog2);
			edges.addAll(hits);
		}
		return edges;
	}
	
	private List<HomologyEdge> findHomologs(HomologyUnit unit, HomologyCatalog catalog) {
		List<HomologyEdge> edges = new ArrayList<HomologyEdge>();
		FMIndex indexCatalog = catalog.getIndexHomologyUnits();
		//Counts of k-mers mapping to each protein in the FM-index
		Map<String,Integer> kmerSupportMap = new TreeMap<>();
		int totalKmers = 0;
		String searchSequence = unit.getUnitSequence();
		//Step 1: Generate k-mers to query the FM-Index looking for homologous transcripts to calculate the kmer counts
		for(int i=0; i<searchSequence.length()-kmerLength+1; i+=kmerLength) {
			String kmer = searchSequence.substring(i, i+kmerLength);
			
			List <UngappedSearchHit> kmerHits = indexCatalog.exactSearch(kmer);
			for(UngappedSearchHit hit:kmerHits) {
				String name = hit.getSequenceName();
				if(kmerSupportMap.containsKey(name)) {
					int value = 1 + kmerSupportMap.get(name);
					kmerSupportMap.put(name, value);
				} else {
					kmerSupportMap.put(name, 1);
				}
			}
			totalKmers++;
		}

		//Step 2: Fill list traversing the counts and choosing transcripts for which at least x% of the k-mers support the match
		for(Map.Entry<String,Integer> entry : kmerSupportMap.entrySet()) {
			String homologId = entry.getKey();
			double transcriptKmers = entry.getValue();
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
	}

}
