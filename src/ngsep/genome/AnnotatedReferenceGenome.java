package ngsep.genome;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.ProteinTranslator;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;

public class AnnotatedReferenceGenome {
	
	private int id;
	private ProteinTranslator translator = new ProteinTranslator();
	private ReferenceGenome genome;
	private Transcriptome transcriptome;
	private QualifiedSequenceList sequencesMetadata;
	private List<OrthologyUnit> orthologyUnitsList;
	private GenomicRegionSortedCollection<OrthologyUnit> orthologyUnitsBySequence;
	private Map<String, OrthologyUnit> orthologyUnitsMap;
	private FMIndex indexOrthologyUnits;
	
	private GenomicRegionSortedCollection<OrthologyUnit> uniqueOrthologyUnitsBySequence;
	
	public AnnotatedReferenceGenome(int id, ReferenceGenome genome, Transcriptome transcriptome) {
		super();
		this.id = id;
		this.genome = genome;
		this.transcriptome = transcriptome;
		this.sequencesMetadata = genome.getSequencesMetadata();
		//Create orthology units based on transcripts
		//Query each orthology unit against its own genome and select the units that are unique
		extractOrthologyUnits();
		//log.info("Genome total units "+orthologyUnitsList.size()+" Building FM Index");
		buildFMIndex ();
	}
	
	/**
	 * Builds the orthology units for this annotated genome based on the largest transcript of each gene
	 */
	private void extractOrthologyUnits() {
		orthologyUnitsList = new ArrayList<>();
		orthologyUnitsMap = new HashMap<>();
		orthologyUnitsBySequence = new GenomicRegionSortedCollection<>(sequencesMetadata);
		List<Transcript> allTranscripts = transcriptome.getAllTranscripts();
		Gene lastGene = null; 
		List<Transcript> transcriptsGene = new ArrayList<>();
		for (Transcript tr:allTranscripts) {	
			Gene gene = tr.getGene();
			if(lastGene != gene) {
				if(lastGene!=null) {
					OrthologyUnit unit = buildOrthologyUnitGene(lastGene, transcriptsGene);
					if (unit!=null) {
						orthologyUnitsList.add(unit);
						orthologyUnitsMap.put(unit.getId(), unit);
						orthologyUnitsBySequence.add(unit);
					}
				}
				lastGene = gene;
				transcriptsGene.clear();
			}
			transcriptsGene.add(tr);
		}
		OrthologyUnit unit = buildOrthologyUnitGene(lastGene, transcriptsGene);
		if (unit!=null) {
			orthologyUnitsList.add(unit);
			orthologyUnitsMap.put(unit.getId(), unit);
			orthologyUnitsBySequence.add(unit);
		}
	}
	
	/**
	 * 
	 * @param gene
	 * @param transcriptsGene
	 * @return Orthology Unit related to the gene or null if the protein can not be assembled for any transcript
	 */
	private OrthologyUnit buildOrthologyUnitGene( Gene gene, List<Transcript> transcriptsGene) {
		String geneId = gene.getId();
		Transcript bestTranscript = null;
		int bestLength = 0;
		String bestProtein = null;
		
		for (Transcript tr:transcriptsGene) {
			String proteinSequence = tr.getProteinSequence(translator);
			if(proteinSequence==null) continue;
			int proteinLength = proteinSequence.length();
			if (proteinLength > bestLength) {
				bestLength = proteinLength;
				bestTranscript = tr;
				bestProtein = proteinSequence;
			}
		}
		
		if(bestTranscript==null) return null;
		int first = bestTranscript.getFirst();
		int last = bestTranscript.getLast();
		String sequenceName = bestTranscript.getSequenceName();
		OrthologyUnit unit = new OrthologyUnit(id, geneId, sequenceName, first, last);
		unit.setUnitSequence(bestProtein);
		return unit;
	}
	
	private void buildFMIndex() {
		indexOrthologyUnits = new FMIndex();
		QualifiedSequenceList unitSequences = new QualifiedSequenceList();
		for (OrthologyUnit ql:orthologyUnitsList) {
			String unitSequence = ql.getUnitSequence();
			String unitId = ql.getId();
			QualifiedSequence qualifiedSequence = new QualifiedSequence(unitId, unitSequence);
			unitSequences.add(qualifiedSequence);
		}
		indexOrthologyUnits.loadQualifiedSequenceList(unitSequences);
	}
	
	public void calculateParalogs(byte kmerSize, int minPctKmers) {
		for (OrthologyUnit unit:orthologyUnitsList) {
			//Orthology unit ids of similar proteins
			Set<String> hits = findOrthologyUnits(indexOrthologyUnits, unit.getUnitSequence(), kmerSize, minPctKmers);
			for(String paralogId:hits) {
				if(!paralogId.equals(unit.getId())) {
					OrthologyUnit paralog = orthologyUnitsMap.get(paralogId);
					unit.addParalog(paralog);
					paralog.addParalog(unit);
				}
			}
		}
		selectUniqueOrthologyUnits();
	}
	
	private void selectUniqueOrthologyUnits () {
		uniqueOrthologyUnitsBySequence = new GenomicRegionSortedCollection<>(sequencesMetadata);
		for(OrthologyUnit unit:orthologyUnitsList) {
			if(unit.isUnique()) uniqueOrthologyUnitsBySequence.add(unit);
		}
	}
	
	private Set<String> findOrthologyUnits(FMIndex index, String searchSequence, byte kmerSize, int minPctKmers) {	
		//Counts of k-mers mapping to each protein in the FM-index
		Map<String,Integer> kmerSupportMap = new TreeMap<>();
		int totalKmers = 0;
		//Step 1: Generate k-mers to query the FM-Index looking for homologous transcripts to calculate the kmer counts
		for(int i=0; i<searchSequence.length()-kmerSize+1; i+=kmerSize) {
			String kmer = searchSequence.substring(i, i+kmerSize);
			
			List <ReadAlignment> kmerHits = index.search(kmer);
			for(ReadAlignment alns:kmerHits) {
				String name = alns.getSequenceName();
				if(kmerSupportMap.containsKey(name))
				{
					int value = 1 + kmerSupportMap.get(name);
					kmerSupportMap.put(name, value);
				}
				else
				{
					kmerSupportMap.put(name, 1);
				}
			}
			totalKmers = totalKmers +1;
		}

		//Step 2: Fill list traversing the counts and choosing transcripts for which at least x% of the k-mers support the match
		Set<String> answer = new TreeSet<>();
		
		for(Map.Entry<String,Integer> entry : kmerSupportMap.entrySet()) {
			String name = entry.getKey();
			double transcriptKmers = entry.getValue();
			double percent = (transcriptKmers/totalKmers)*100;
			if(percent >= minPctKmers)
			{
				answer.add(name);
			}
		}
		return answer;
	}
	
	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}
	/**
	 * Finds orthologs of the orthology units in this genome in the given genome
	 * @param genome2 to search for orthologs
	 */
	public void calculateOrthologs (AnnotatedReferenceGenome genome2, byte kmerSize, int minPctKmers) {
		for (OrthologyUnit unit:orthologyUnitsList) {
			//Orthology unit ids of similar proteins
			Set<String> hits = findOrthologyUnits(genome2.indexOrthologyUnits, unit.getUnitSequence(), kmerSize, minPctKmers);
			for(String orthologId:hits) {
				OrthologyUnit ortholog = genome2.orthologyUnitsMap.get(orthologId);
				unit.addOrtholog(ortholog);
				ortholog.addOrtholog(unit);
			}
		}
	}
	/**
	 * @return List<OrthologyUnit> Orthology units for this genome
	 */
	public List<OrthologyUnit> getOrthologyUnits() {
		return Collections.unmodifiableList(orthologyUnitsList);
	}
	
	/**
	 * Returns the list of all orthology units for the given sequence name
	 * @param name Sequence name
	 * @return List<OrthologyUnit> units with the given sequence name
	 */
	public List<OrthologyUnit> getOrthologyUnits(String name) {
		return orthologyUnitsBySequence.getSequenceRegions(name).asList();
	}
	
	/**
	 * @return List<OrthologyUnit> List of unique orthology units in this genome
	 */
	public List<OrthologyUnit> getUniqueOrthologyUnits() {
		return uniqueOrthologyUnitsBySequence.asList();
	}
	
	/**
	 * Returns the list of unique orthology units for the given sequence name
	 * @param name Sequence name
	 * @return List<OrthologyUnit> unique units with the given sequence name
	 */
	public List<OrthologyUnit> getUniqueOrthologyUnits(String name) {
		return uniqueOrthologyUnitsBySequence.getSequenceRegions(name).asList();
	}

	/**
	 * @return QualifiedSequenceList Metadata of the sequences in this genome
	 */
	public QualifiedSequenceList getSequencesMetadata() {
		return sequencesMetadata;
	}
	/**
	 * @param sequenceName
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getSequenceByName(java.lang.String)
	 */
	public QualifiedSequence getSequenceByName(String sequenceName) {
		return genome.getSequenceByName(sequenceName);
	}
	/**
	 * @param index
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getSequenceByIndex(int)
	 */
	public QualifiedSequence getSequenceByIndex(int index) {
		return genome.getSequenceByIndex(index);
	}
	/**
	 * @param sequenceName
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getSequenceCharacters(java.lang.String)
	 */
	public CharSequence getSequenceCharacters(String sequenceName) {
		return genome.getSequenceCharacters(sequenceName);
	}
	/**
	 * @param index
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getSequenceCharacters(int)
	 */
	public CharSequence getSequenceCharacters(int index) {
		return genome.getSequenceCharacters(index);
	}
	/**
	 * @param r
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getReference(ngsep.genome.GenomicRegion)
	 */
	public CharSequence getReference(GenomicRegion r) {
		return genome.getReference(r);
	}
	/**
	 * @param sequenceName
	 * @param first
	 * @param last
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getReference(java.lang.String, int, int)
	 */
	public CharSequence getReference(String sequenceName, int first, int last) {
		return genome.getReference(sequenceName, first, last);
	}
	/**
	 * @param index
	 * @param first
	 * @param last
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getReference(int, int, int)
	 */
	public CharSequence getReference(int index, int first, int last) {
		return genome.getReference(index, first, last);
	}
	/**
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getTotalLength()
	 */
	public long getTotalLength() {
		return genome.getTotalLength();
	}
	/**
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getNumSequences()
	 */
	public int getNumSequences() {
		return genome.getNumSequences();
	}
	/**
	 * @return
	 * @see ngsep.genome.ReferenceGenome#getSequenceNamesStringList()
	 */
	public List<String> getSequenceNamesStringList() {
		return genome.getSequenceNamesStringList();
	}
}
