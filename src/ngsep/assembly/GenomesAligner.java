package ngsep.assembly;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegion;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.math.Distribution;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.ProteinTranslator;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;

public class GenomesAligner {

	private ReferenceGenome genome1;
	private ReferenceGenome genome2;
	private Transcriptome transcriptome1;
	private Transcriptome transcriptome2;
	private ProteinTranslator translator = new ProteinTranslator();
	
	
	
	public static void main(String[] args) throws Exception 
	{
		GenomesAligner instance = new GenomesAligner();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String fileGenome1 = args[i++];
		String fileTranscriptome1 = args[i++];
		String fileGenome2 = args[i++];
		String fileTranscriptome2 = args[i++];
		String outPrefix = args[i++]; 
		instance.loadGenomes(fileGenome1,fileTranscriptome1, fileGenome2, fileTranscriptome2);
		instance.findUniqueCopyOrthologs();
		instance.alignOrthologUnits();
		try (PrintStream outAlignmnent = new PrintStream(outPrefix+"_lcs.txt");) {
			instance.printAlignmentOrthologUnits (outAlignmnent);
		}
	}

	
	
	
	public void loadGenomes(String fileGenome1, String fileTranscriptome1, String fileGenome2, String fileTranscriptome2) throws IOException {
		genome1 = new ReferenceGenome(fileGenome1);
		genome2 = new ReferenceGenome(fileGenome2);
		GFF3TranscriptomeHandler transcriptomeHandler = new GFF3TranscriptomeHandler();
		transcriptome1 = transcriptomeHandler.loadMap(fileTranscriptome1);
		transcriptome1.fillSequenceTranscripts(genome1);
		transcriptome2 = transcriptomeHandler.loadMap(fileTranscriptome2);
		transcriptome2.fillSequenceTranscripts(genome2);
		
	}
	
	
	public void findUniqueCopyOrthologs() {
		//Create orthology units based on transcripts
		
		List<OrthologyUnit> unitsT1 = extractOrthologyUnits(transcriptome1);
		List<OrthologyUnit> unitsT2 = extractOrthologyUnits(transcriptome2);
		QualifiedSequenceList proteinSequencesT1 = extractProteinSequences(unitsT1);
		QualifiedSequenceList proteinSequencesT2 = extractProteinSequences(unitsT2);	

		
		//Build two FMIndex objects for proteins in the two genomes
		
		//FMIndex for qualified sequence list for genome 1
		FMIndex FMIndex1 = new FMIndex();
		FMIndex1.loadQualifiedSequenceList(proteinSequencesT1);
		
		
		// FMIndex for qualified sequence list for genome 2
		FMIndex FMIndex2 = new FMIndex();
		FMIndex2.loadQualifiedSequenceList(proteinSequencesT2);
		
		
		

		
		//Query each orthology unit against its own genome and select the units that are unique
		QualifiedSequenceList uniquesUnitsT1 = searchUniqueUnits(FMIndex1, proteinSequencesT1);
		System.out.println("Total proteins 1: "+proteinSequencesT1.size()+" Unique units T1: "+uniquesUnitsT1.size());
		
		QualifiedSequenceList uniquesUnitsT2 = searchUniqueUnits(FMIndex2, proteinSequencesT2);
		System.out.println("Total proteins 2: "+proteinSequencesT2.size()+" Unique units T2: "+uniquesUnitsT2.size());
		
		
		
		//Query selected units against the other genome to assign mates
		QualifiedSequenceList mates1 = searchUniqueUnits(FMIndex1, uniquesUnitsT2);
		System.out.println("Mates T1: "+mates1.size());
		QualifiedSequenceList mates2 = searchUniqueUnits(FMIndex2, uniquesUnitsT1);
		System.out.println("Mates T2: "+mates2.size());
	}
	

	
	public QualifiedSequenceList searchUniqueUnits (FMIndex fmindex, QualifiedSequenceList proteinSequences)
	{
	
		QualifiedSequenceList uniques = new QualifiedSequenceList();
		
		for (int i=0; i<proteinSequences.size(); i++)
		{
			QualifiedSequence next = proteinSequences.get(i);
			String nextSequence = next.getCharacters().toString();
			
			
			List<ReadAlignment> searchResult = fmindex.search(nextSequence);
			if(searchResult.size() == 1)
			{
				uniques.add(next);
			}
		}
		
		return uniques;
	}

	

	/**
	 * @param transcriptome
	 * @return
	 */
	private List<OrthologyUnit> extractOrthologyUnits(Transcriptome transcriptome) {
		List<OrthologyUnit> orthologyUnits = new ArrayList<>();
		List<Transcript> allTranscripts = transcriptome.getAllTranscripts();
		Gene lastGene = null; 
		List<Transcript> transcriptsGene = new ArrayList<>();
		
		
		for (Transcript tr:allTranscripts)
		{	
			Gene gene = tr.getGene();
			if(lastGene != gene) 
			{
				if(lastGene!=null) 
				{
					OrthologyUnit unit = buildOrthologyUnitGene(lastGene, transcriptsGene);
					if (unit!=null) orthologyUnits.add(unit);
				}
				lastGene = gene;
				transcriptsGene.clear();
			}
			transcriptsGene.add(tr);
		}
		OrthologyUnit unit = buildOrthologyUnitGene(lastGene, transcriptsGene);
		if (unit!=null) orthologyUnits.add(unit);
		return orthologyUnits;
	}
	
	private QualifiedSequenceList extractProteinSequences(List<OrthologyUnit> units) {
		QualifiedSequenceList proteinSequences = new QualifiedSequenceList();
		
		for (OrthologyUnit ql:units)
		{
			
			String proteinSequence = ql.getProteinSequence();
			String proteinId = ql.getId();
			QualifiedSequence qualifiedSequence = new QualifiedSequence(proteinId, proteinSequence);
			
			proteinSequences.add(qualifiedSequence);
		}
		
		return proteinSequences;
	}
	/**
	 * 
	 * @param gene
	 * @param transcriptsGene
	 * @return Orthology Unit related to the gene or null if the protein can not be assembled
	 */
	private OrthologyUnit buildOrthologyUnitGene( Gene gene, List<Transcript> transcriptsGene) {
		String geneId = gene.getId();
		Transcript bestTranscript = chooseBestTranscript(transcriptsGene);
		String proteinSequence = bestTranscript.getProteinSequence(translator);
		if(proteinSequence==null)
		{
			return null;
		}
		int first = bestTranscript.getFirst();
		int last = bestTranscript.getLast();
		String sequenceName = bestTranscript.getSequenceName();
		OrthologyUnit unit = new OrthologyUnit(geneId, sequenceName, first, last);
		unit.setProteinSequence(proteinSequence);
		
		
		return unit;
	}



	private Transcript chooseBestTranscript(List<Transcript> transcriptsGene) {
	
		
		Transcript bestTranscript = null;
		int longestTranscript = 0;
		
		for (Transcript tr:transcriptsGene)
		{
			int trLength = tr.length();
			if (trLength > longestTranscript)
			{
				longestTranscript = trLength;
				bestTranscript = tr;
			}
		}
		
		return bestTranscript;
	}
//probar esto y todos los metodos que haga 
	public void alignOrthologUnits() {
		// TODO: Implement LCS algorithm to align ortholog units
		
	}



	private void printAlignmentOrthologUnits(PrintStream outAlignmnent) {
		// TODO Auto-generated method stub
		
	}

}
class OrthologyUnit implements GenomicRegion {
	private String id;
	private String sequenceName;
	private int first;
	private int last;
	private boolean negativeStrand = false;
	private String proteinSequence;
	private String mateId;
	
	public OrthologyUnit(String id, String sequenceName, int first, int last) {
		super();
		this.id = id;
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
	}
	

	public String getId() {
		return id;
	}

	@Override
	public String getSequenceName() {
		return sequenceName;
	}
	@Override
	public int getFirst() {
		return first;
	}
	@Override
	public int getLast() {
		return last;
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return !negativeStrand;
	}
	@Override
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	/**
	 * @param negativeStrand the negativeStrand to set
	 */
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
	/**
	 * @return the proteinSequence
	 */
	public String getProteinSequence() {
		return proteinSequence;
	}
	/**
	 * @param proteinSequence the proteinSequence to set
	 */
	public void setProteinSequence(String proteinSequence) {
		this.proteinSequence = proteinSequence;
	}
	/**
	 * @return the mateId
	 */
	public String getMateId() {
		return mateId;
	}
	/**
	 * @param mateId the mateId to set
	 */
	public void setMateId(String mateId) {
		this.mateId = mateId;
	}
}
