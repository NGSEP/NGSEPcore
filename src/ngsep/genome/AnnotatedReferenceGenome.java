package ngsep.genome;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

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
	private List<HomologyUnit> homologyUnitsList;
	private GenomicRegionSortedCollection<HomologyUnit> homologyUnitsBySequence;
	private HomologyCatalog catalog;
	
	private GenomicRegionSortedCollection<HomologyUnit> uniqueHomologyUnitsBySequence;
	
	public AnnotatedReferenceGenome(int id, ReferenceGenome genome, Transcriptome transcriptome) {
		super();
		this.id = id;
		this.genome = genome;
		this.transcriptome = transcriptome;
		this.sequencesMetadata = genome.getSequencesMetadata();
		//Create orthology units based on transcripts
		//Query each orthology unit against its own genome and select the units that are unique
		extractHomologyUnits();
		//log.info("Genome total units "+orthologyUnitsList.size()+" Building FM Index");
		catalog = new HomologyCatalog(homologyUnitsList); 
	}
	
	/**
	 * Builds the homology units for this annotated genome based on the largest transcript of each gene
	 */
	private void extractHomologyUnits() {
		homologyUnitsList = new ArrayList<>();
		homologyUnitsBySequence = new GenomicRegionSortedCollection<>(sequencesMetadata);
		List<Transcript> allTranscripts = transcriptome.getAllTranscripts();
		Gene lastGene = null; 
		List<Transcript> transcriptsGene = new ArrayList<>();
		for (Transcript tr:allTranscripts) {	
			Gene gene = tr.getGene();
			if(lastGene != gene) {
				if(lastGene!=null) {
					HomologyUnit unit = buildHomologyUnitGene(lastGene, transcriptsGene);
					if (unit!=null) {
						homologyUnitsList.add(unit);
						homologyUnitsBySequence.add(unit);
					}
				}
				lastGene = gene;
				transcriptsGene.clear();
			}
			transcriptsGene.add(tr);
		}
		HomologyUnit unit = buildHomologyUnitGene(lastGene, transcriptsGene);
		if (unit!=null) {
			homologyUnitsList.add(unit);
			homologyUnitsBySequence.add(unit);
		}
		
	}
	
	/**
	 * 
	 * @param gene
	 * @param transcriptsGene
	 * @return Orthology Unit related to the gene or null if the protein can not be assembled for any transcript
	 */
	private HomologyUnit buildHomologyUnitGene( Gene gene, List<Transcript> transcriptsGene) {
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
		HomologyUnit unit = new HomologyUnit(id, geneId, sequenceName, first, last);
		unit.setUnitSequence(bestProtein);
		return unit;
	}
	
	public void selectUniqueOrthologyUnits () {
		uniqueHomologyUnitsBySequence = new GenomicRegionSortedCollection<>(sequencesMetadata);
		for(HomologyUnit unit:homologyUnitsList) {
			if(unit.isUnique()) uniqueHomologyUnitsBySequence.add(unit);
		}
	}
	
	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}
	
	/**
	 * @return List<HomologyUnit> Homology units for this genome
	 */
	public List<HomologyUnit> getHomologyUnits() {
		return Collections.unmodifiableList(homologyUnitsList);
	}
	
	/**
	 * @return int total number of homology units
	 */
	public int getTotalHomologyUnits () {
		return homologyUnitsList.size();
	}
	
	/**
	 * Returns the list of all homology units for the given sequence name
	 * @param name Sequence name
	 * @return List<HomologyUnit> units with the given sequence name
	 */
	public List<HomologyUnit> getHomologyUnits(String name) {
		return homologyUnitsBySequence.getSequenceRegions(name).asList();
	}
	
	/**
	 * Returns an FM-index of the homology units
	 * @return FMIndex to search for homologs
	 */
	public FMIndex getIndexHomologyUnits() {
		return catalog.getIndexHomologyUnits();
	}
	
	/**
	 * Return the homology catalog for this reference genome
	 * @return OrganismHomologyCatalog to search for homology units
	 */
	public HomologyCatalog getHomologyCatalog() {
		return this.catalog;
	}
	
	public HomologyUnit getHomologyUnit(String unitId) {
		return catalog.getHomologyUnit(unitId);
	}

	/**
	 * @return List<HomologyUnit> List of unique homology units in this genome
	 */
	public List<HomologyUnit> getUniqueHomologyUnits() {
		return uniqueHomologyUnitsBySequence.asList();
	}
	
	/**
	 * @return int total number of unique homology units within the genome
	 */
	public int getNumberUniqueHomologyUnits () {
		return uniqueHomologyUnitsBySequence.size();
	}
	
	/**
	 * Returns the list of unique homology units for the given sequence name
	 * @param name Sequence name
	 * @return List<HomologyUnit> unique units with the given sequence name
	 */
	public List<HomologyUnit> getUniqueHomologyUnits(String name) {
		return uniqueHomologyUnitsBySequence.getSequenceRegions(name).asList();
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
