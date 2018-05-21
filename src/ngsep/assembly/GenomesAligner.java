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
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.ProteinTranslator;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;

/**
 * @author Daniel Tello
 * @author Jorge Duitama
 */
public class GenomesAligner {

	private ReferenceGenome genome1;
	private ReferenceGenome genome2;
	private Transcriptome transcriptome1;
	private Transcriptome transcriptome2;
	private ProteinTranslator translator = new ProteinTranslator();
	private List<OrthologyUnit> alignedUnits = new ArrayList<>();
	private Logger log = Logger.getLogger(GenomesAligner.class.getName());
	private ProgressNotifier progressNotifier=null;
	
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
	
	/**
	 * @return the log
	 */
	public Logger getLog() {
		return log;
	}
	
	/**
	 * @param log the log to set
	 */
	public void setLog(Logger log) {
		this.log = log;
	}

	/**
	 * @return the progressNotifier
	 */
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}

	/**
	 * @param progressNotifier the progressNotifier to set
	 */
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}

	public void loadGenomes(String fileGenome1, String fileTranscriptome1, String fileGenome2, String fileTranscriptome2) throws IOException {
		genome1 = new ReferenceGenome(fileGenome1);
		log.info("Loaded genome "+fileGenome1);
		genome2 = new ReferenceGenome(fileGenome2);
		log.info("Loaded genome "+fileGenome2);
		GFF3TranscriptomeHandler transcriptomeHandler = new GFF3TranscriptomeHandler();
		transcriptome1 = transcriptomeHandler.loadMap(fileTranscriptome1);
		transcriptome1.fillSequenceTranscripts(genome1);
		log.info("Loaded transcriptome "+fileTranscriptome1);
		transcriptome2 = transcriptomeHandler.loadMap(fileTranscriptome2);
		transcriptome2.fillSequenceTranscripts(genome2);
		log.info("Loaded transcriptome "+fileTranscriptome2);
	}
	
	
	public void alignGenomes() {
		//Create orthology units based on transcripts
		
		List<OrthologyUnit> unitsT1 = extractOrthologyUnits(transcriptome1);
		
		List<OrthologyUnit> unitsT2 = extractOrthologyUnits(transcriptome2);
		
		
		//Query each orthology unit against its own genome and select the units that are unique
		List<OrthologyUnit>  uniquesUnitsT1 = selectUniqueUnits(unitsT1);
		log.info("First genome has "+unitsT1.size()+" total orthology units. Unique units: "+uniquesUnitsT1.size());
		List<OrthologyUnit> uniquesUnitsT2 = selectUniqueUnits(unitsT2);
		log.info("Second genome has "+unitsT2.size()+" total orthology units. Unique units: "+uniquesUnitsT2.size());
		
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
				log.info("Sequence "+chrG1.getName()+" in first genome aligned to sequence "+chrG2.getName()+" in the second genome. Orthology units sequence genome 1 "+unitsChrG1.size()+". Orthology units sequence genome 2: "+unitsChrG2.size()+" LCS size: "+selectedUnits.size());
				alignedUnits.addAll(selectedUnits);
			} else {
				log.info("Mate sequence not found for "+chrG1.getName()+" Sequence orthology units: "+unitsChrG1.size());
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
			// If the hit is unique fill the mate information using setMateId from OrthologyUnit
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
		
		
		Set<Integer> lcsForward = findLCS(unitsG1List, unitsG2List);
		
		Collections.reverse(unitsG2List);
		
		Set<Integer> lcsReverse = findLCS(unitsG1List, unitsG2List);
		
		Set<Integer> lcs = lcsForward;
		if(lcsReverse.size()>lcsForward.size()) lcs = lcsReverse;
		//System.out.println("Positions for LCS: "+positions.length+" LCS: "+lcs.size());
		// Select the orthology units in G1 located at the indexes given by the output of LCS
		for(int i:lcs)
		{
			OrthologyUnit lcsResult = unitsG1List.get(i);
			answer.add(lcsResult);
		}
		
		
		return answer;
	}

	private Set<Integer> findLCS(List<OrthologyUnit> unitsG1List, List<OrthologyUnit> unitsG2List) {
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
			//System.out.println("Positions [ "+i+"]:"+j);
		}
		// Run LCS
		
		Set<Integer> lcs = findLCS(positions);
		return lcs;
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
			if(i%100==0)log.info("Processed "+i+" proteins. Unique: "+uniques.size());
		}
		
		return uniques;
	}
	
	private static Set<String> findSimilarProteins(FMIndex index, String protein1) {
		
		//Counts of k-mers mapping to each protein in the FM-index
		Map<String,Integer> kmerSupportMap = new TreeMap<>();
		int totalKmers = 0;
		int kmerSize = 10;
		//Step 1: Generate k-mers to query the FM-Index looking for homologous transcripts to calculate the kmer counts
		for(int i=0; i<protein1.length()-kmerSize+1; i+=kmerSize) {
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

		//Step 2: Fill list traversing the counts and choosing transcripts for which at least 80% of the k-mers support the match
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
				if(i==0)
				{
					if(j==0 || indexesMap[i]>0) {
						m[i][j] = 0;
					}
					else {
						m[i][j] = 1;
					}
					
				}
				else if(j<=indexesMap[i])
				{
					m[i][j] = m[i-1][j];
				}
				else
				{
					m[i][j] = Math.max((m[i-1][j]),(m[i-1][indexesMap[i]]+1));
				}
				//System.out.print(" "+m[i][j]);
			}
			//System.out.println();
		}

		int i = m.length-1;
		int j = m[0].length-1;
		//System.out.println("LCS matrix size: "+i+"-"+j);
		while(i > 0 && j > 0)
		{
			//System.out.print("Position: "+i+"-"+j);
			int up = m[i-1][j];
			int diag = -1;
			if(j>indexesMap[i]) {
				diag = m[i-1][indexesMap[i]]+1;
			}
			if(diag >= up )
			{
				answer.add(i);
				j=indexesMap[i];
				i--;
			} else {
				i--;
			}
			//System.out.println(" Diag score: "+diag+" up score: "+up+" size answer: "+answer.size()+" next Position: "+i+"-"+j);
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
