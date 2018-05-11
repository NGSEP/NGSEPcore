package ngsep.assembly;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
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
	private List<OrthologyUnit> alignedUnits = new ArrayList<>();
	
	
	
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
		instance.alignGenomes();
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
	
	
	public void alignGenomes() {
		//Create orthology units based on transcripts
		
		List<OrthologyUnit> unitsT1 = extractOrthologyUnits(transcriptome1);
		List<OrthologyUnit> unitsT2 = extractOrthologyUnits(transcriptome2);
		
		//Query each orthology unit against its own genome and select the units that are unique
		List<OrthologyUnit>  uniquesUnitsT1 = selectUniqueUnits(unitsT1);
		System.out.println("Total units 1: "+unitsT1.size()+" Unique units T1: "+uniquesUnitsT1.size());
		List<OrthologyUnit> uniquesUnitsT2 = selectUniqueUnits(unitsT2);
		System.out.println("Total proteins 2: "+unitsT2.size()+" Unique units T2: "+uniquesUnitsT2.size());
		
		QualifiedSequenceList genomeMetadata1 = genome1.getSequencesMetadata();
		GenomicRegionSortedCollection<OrthologyUnit> unitsC1 = new GenomicRegionSortedCollection<>(genomeMetadata1);
		unitsC1.addAll(uniquesUnitsT1);
		QualifiedSequenceList genomeMetadata2 = genome2.getSequencesMetadata();
		GenomicRegionSortedCollection<OrthologyUnit> unitsCollectionG2 = new GenomicRegionSortedCollection<>(genomeMetadata2);
		unitsCollectionG2.addAll(uniquesUnitsT2);
		alignedUnits = new ArrayList<>();
		
		for(QualifiedSequence chrG1:genomeMetadata1) {
			List<OrthologyUnit> unitsChrG1 = unitsC1.getSequenceRegions(chrG1.getName()).asList();
			fillMateIds(unitsChrG1,unitsCollectionG2);
			String chrNameG2 = findBestChromosome(unitsChrG1,unitsCollectionG2);
			if(chrNameG2!=null) {
				QualifiedSequence chrG2 = genomeMetadata2.get(chrNameG2);
				List<OrthologyUnit> unitsChrG2 = unitsCollectionG2.getSequenceRegions(chrG2.getName()).asList();
				List<OrthologyUnit> selectedUnits = alignOrthologyUnits(unitsChrG1,unitsChrG2);
				System.out.println("Chr g1: "+chrG1.getName()+" chrG2: "+chrG2.getName()+" units G1 "+unitsChrG1.size()+" units G2: "+unitsChrG2.size()+" LCS size: "+selectedUnits.size());
				alignedUnits.addAll(selectedUnits);
			} else {
				System.out.println("Mate chromosome not found for "+chrG1.getName()+" units G1 "+unitsChrG1.size());
			}
		}
	}
	
	/**
	 * FInds the unique mate in the second genome of the units in the first genome
	 * @param unitsG1
	 * @param unitsG2
	 */
	private void fillMateIds(List<OrthologyUnit> unitsG1, GenomicRegionSortedCollection<OrthologyUnit> unitsG2) {
		Map<String,OrthologyUnit> unitsC2Map = new HashMap<>();
		for(OrthologyUnit unit:unitsG2) unitsC2Map.put(unit.getId(), unit);
		
		//Build FMindex g2
		
		QualifiedSequenceList proteinSequencesG2 = extractProteinSequences(unitsG2.asList());
		FMIndex indexG2 = new FMIndex();
		indexG2.loadQualifiedSequenceList(proteinSequencesG2);
		//Search proteins from G1 using method findSimilarProteins

		for(int i=0; i<unitsG1.size(); i++)
		{
			OrthologyUnit unit = unitsG1.get(i);
			String protein = unit.getProteinSequence();
			Set<String> hits = findSimilarProteins(indexG2, protein);
		//If the hit is unique fill the mate information using setMateId from OrthologyUnit
			if(hits.size() == 1)
			{
				Iterator<String> iterator = hits.iterator();
				String mateId = iterator.next();
				OrthologyUnit mate = unitsC2Map.get(mateId);
				unit.setMate(mate);
			}
		}
	}
		



	/**
	 * Selects the chromosome in unitsC2 having the largest number of mates with the given units from genome 1
	 * @param unitsG1
	 * @param unitsC2
	 * @return String name of the chromosome in the second genome
	 */
	private String findBestChromosome(List<OrthologyUnit> unitsG1, GenomicRegionSortedCollection<OrthologyUnit> unitsG2) {
		
		
		Map<String,Integer> chrG2Counts = new HashMap<>();
		

		//Go over the orthology units of G1. Locate chromosome of each mate and update counts map
		
		for(int i=0; i<unitsG1.size(); i++)
		{
			OrthologyUnit unitG2 = unitsG1.get(i).getMate();
			if(unitG2 == null) continue;
			//Recupero la unidad de ortología del mapa unitsC2Map
			//Recupero el cromosoma
			String sequenceName = unitG2.getSequenceName();			
			if((chrG2Counts.containsKey(sequenceName)))
			{
				int value = chrG2Counts.get(sequenceName)+1;
				chrG2Counts.put(sequenceName, value);
			}
			else
			{
				chrG2Counts.put(sequenceName, 1);
			}
				
		}
		//Find the chromosome with the largest count
		int bestChromosomeCounts = 0;
		String bestChromosome = null;
		
		for(Map.Entry<String,Integer> entry : chrG2Counts.entrySet()) {
			String name = entry.getKey();
			int chromosomeProteins = entry.getValue();
			
			if(chromosomeProteins > bestChromosomeCounts)
			{
				bestChromosomeCounts = chromosomeProteins;
				bestChromosome = name;
			}
		}

		return bestChromosome;
	}




	/**
	 * Aligns the orthology units from two homologous chromosomes using LCS
	 * @param unitsChrG1 This list has only one chromosome and is sorted by position
	 * @param unitsChrG2 This list has only one chromosome and is sorted by position
	 * @return List<OrthologyUnit> List of selected units of the first list making the LCS relative to the second list
	 */
	private List<OrthologyUnit> alignOrthologyUnits(List<OrthologyUnit> unitsChrG1, List<OrthologyUnit> unitsChrG2) {
		List<OrthologyUnit> answer = new ArrayList<>();
		
		List<OrthologyUnit> unitsG1List = new ArrayList<>();
		List<OrthologyUnit> unitsG2List = new ArrayList<>();
		//Assigns mates of units in G2
		
		// Select orthology units in G1 having mate in g2 and assign the mate of units in G2.
	
		//At the same time create two new lists of the same size with the  units in g1 having its mate in g2
		
		for(int i=0; i<unitsChrG1.size(); i++)
		{
			OrthologyUnit unitG1 = unitsChrG1.get(i);
			OrthologyUnit unitG2 = unitG1.getMate();
			if(unitG2==null) continue;
			unitsG1List.add(unitG1);			
			unitsG2List.add(unitG2);
			unitG2.setMate(unitG1);
		}
		Collections.sort(unitsG2List, GenomicRegionPositionComparator.getInstance());		
		
		// Create int array with the positions in g2 for the units in g1. Input for LCS
		Map<String,Integer> sortedUnitsG2Pos = new HashMap<>();
		for(int i=0; i<unitsG2List.size(); i++) {
			OrthologyUnit unit=unitsG2List.get(i);
			sortedUnitsG2Pos.put(unit.getId(), i);
		}
		int []positions = new int[unitsG1List.size()];
		for(int i=0; i<unitsG1List.size(); i++)
		{
			OrthologyUnit unitG1 = unitsG1List.get(i);
			int j = sortedUnitsG2Pos.get(unitG1.getMate().getId());
			positions[i] = j;
		}
		// Run LCS
		
		Set<Integer> lcs = findLCS(positions);
		// Select the orthology units in G1 located at the indexes given by the output of LCS
		for(int i:lcs)
		{
			OrthologyUnit lcsResult = unitsG1List.get(i);
			answer.add(lcsResult);
		}
		
		
		return answer;
	}




	public List<OrthologyUnit> selectUniqueUnits (List<OrthologyUnit> units)
	{
	
		List<OrthologyUnit> uniques = new ArrayList<>();
		QualifiedSequenceList proteinSequences = extractProteinSequences(units);
		FMIndex index = new FMIndex();
		index.loadQualifiedSequenceList(proteinSequences);
		
		for (int i=0; i<proteinSequences.size(); i++)
		{
			QualifiedSequence next = proteinSequences.get(i);
			String nextSequence = next.getCharacters().toString();
			
			
			Set<String> hits = findSimilarProteins(index, nextSequence);
			if(hits.size() == 1)
			{
				uniques.add(units.get(i));
			}
		}
		
		return uniques;
	}
	
	private static Set<String> findSimilarProteins(FMIndex index, String protein1) {
		
		//el mapa es transcript y cuantos kmers lo soporta
		Map<String,Integer> kmerSupportMap = new TreeMap<>();
		int totalKmers = 0;
		int kmerSize = 10;
		//PASO 1: Generar n kmers y por cada uno consultar en el fm index los transcripts que tienen el kmer para llenar el mapa
		
		for(int i=0; i<protein1.length()-kmerSize+1; i++) {
			String kmer = protein1.substring(i, i+kmerSize);
			
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

		
		//PASO 2: Llenar la lista de respuesta recorriendo el mapa y escogiendo los ids de transcripts que
		//para los que al menos el 80% de totalKmers dice que el transcript es paralogo de protein1
		Set<String> answer = new TreeSet<>();
		
		for(Map.Entry<String,Integer> entry : kmerSupportMap.entrySet()) {
			String name = entry.getKey();
			double transcriptKmers = entry.getValue();
			double percent = (transcriptKmers/totalKmers)*100;
			if(percent >= 80)
			{
				answer.add(name);
			}
			
		}
		
		return answer;
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
	
	public SortedSet<Integer> findLCS (int [] indexesMap) {
		SortedSet<Integer> answer = new TreeSet<>();
		int n = indexesMap.length;
		int [] [] m = new int [n][n+1];
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<n+1; j++)
			{
				
				if(i==0 && indexesMap[i]==0)
				{
					m[i][j] = 1;
				}
				else if(i==0 && indexesMap[i]>0)
				{
					m[i][j] = 0;
				}
				else if(i>0 && j<=indexesMap[i])
				{
					m[i][j] = m[i-1][j];
				}
				else if(i>0 && j>indexesMap[i])
				{
					m[i][j] = Math.max((m[i-1][j]),(m[i-1][indexesMap[i]]+1));
				}
			}
		}

		
		//Hacer debugging a esta parte para ver dónde se estan perdiendo los datos 
		int i = m.length-1;
		int j = m[0].length-1;
		System.out.println("LCS matrix size: "+i+"-"+j);
		while(i > 0 && j > 0)
		{
			int up = m[i-1][j];
			int diag = -1;
			if(j>indexesMap[i]) {
				diag = m[i-1][indexesMap[i]]+1;
			}
			if(diag >= up )
			{
				answer.add(i);
				i--;
				j=indexesMap[i];
			} else {
				i--;
			}
		}
		if(i==0 && indexesMap[i]==0)
		{
			answer.add(i);
		}
		return answer;
	}
	
	
	private void printAlignmentOrthologUnits(PrintStream outAlignmnent) {
		for(OrthologyUnit unit:alignedUnits) {
			outAlignmnent.print(unit.getId()+"\t"+unit.getSequenceName()+"\t"+unit.getFirst()+"\t"+unit.getLast());
			OrthologyUnit mate = unit.getMate();
			if(mate != null) outAlignmnent.println("\t"+mate.getId()+"\t"+mate.getSequenceName()+"\t"+mate.getFirst()+"\t"+mate.getLast());
			else outAlignmnent.println("\t-\t-\t-\t-");
		}
	}

}
class OrthologyUnit implements GenomicRegion {
	private String id;
	private String sequenceName;
	private int first;
	private int last;
	private boolean negativeStrand = false;
	private String proteinSequence;
	private OrthologyUnit mate;
	
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


	public OrthologyUnit getMate() {
		return mate;
	}


	public void setMate(OrthologyUnit mate) {
		this.mate = mate;
	}

}
