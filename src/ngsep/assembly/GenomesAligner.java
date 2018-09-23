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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
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

	public static final String DEF_OUT_PREFIX = "genomesAlignment";
	
	private Logger log = Logger.getLogger(GenomesAligner.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private List<AnnotatedReferenceGenome> genomes = new ArrayList<>();
	private String outPrefix = DEF_OUT_PREFIX;
	
	private List<OrthologyUnit> alignedUnitsFirstGenome = new ArrayList<>();
	
	public static void main(String[] args) throws Exception 
	{
		GenomesAligner instance = new GenomesAligner();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		while(i<args.length-1) {
			String fileGenome = args[i++];
			String fileTranscriptome = args[i++];
			instance.loadGenome(fileGenome, fileTranscriptome);
		}
		instance.alignGenomes();
		instance.printAlignmentResults ();
		instance.log.info("Process finished");
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
	
	

	/**
	 * @return the outPrefix
	 */
	public String getOutPrefix() {
		return outPrefix;
	}

	/**
	 * @param outPrefix the outPrefix to set
	 */
	public void setOutPrefix(String outPrefix) {
		this.outPrefix = outPrefix;
	}

	public void loadGenome(String fileGenome, String fileTranscriptome) throws IOException {
		ReferenceGenome genome = new ReferenceGenome(fileGenome);
		log.info("Loaded genome "+fileGenome);
		GFF3TranscriptomeHandler transcriptomeHandler = new GFF3TranscriptomeHandler(genome.getSequencesMetadata());
		Transcriptome transcriptome = transcriptomeHandler.loadMap(fileTranscriptome);
		transcriptome.fillSequenceTranscripts(genome);
		log.info("Loaded transcriptome "+fileTranscriptome+ " number of transcripts: "+transcriptome.getAllTranscripts().size());
		AnnotatedReferenceGenome annGenome = new AnnotatedReferenceGenome(genomes.size(), genome, transcriptome);
		log.info("Genome: "+annGenome.getId()+" has "+annGenome.getOrthologyUnits().size()+" total orthology units. Unique: "+annGenome.getUniqueOrthologyUnits().size());
		genomes.add(annGenome);
	}
	
	
	public void alignGenomes() {
		if(genomes.size()==0) {
			log.severe("At least one genome is required");
			return;
		}
		
		for(int i=0;i<genomes.size();i++) {
			for (int j=0;j<genomes.size();j++) {
				if(i!=j) genomes.get(i).findOrthologs(genomes.get(j));
			}
		}
		// By now this is still done for two genomes
		AnnotatedReferenceGenome genome1 = genomes.get(0);
		AnnotatedReferenceGenome genome2 = genomes.get(1);
		alignedUnitsFirstGenome = new ArrayList<>();
		QualifiedSequenceList sequencesG1 = genome1.getSequencesMetadata();
		for(QualifiedSequence chrG1:sequencesG1) {
			List<OrthologyUnit> unitsChrG1 = genome1.getUniqueOrthologyUnits(chrG1.getName());
			log.info("Unique units G1 for "+chrG1.getName()+": "+unitsChrG1.size());
			String chrNameG2 = findBestChromosome(unitsChrG1, genome2.getId());
			if(chrNameG2!=null) {
				List<OrthologyUnit> selectedUnits = alignOrthologyUnits(genome1.getId(),unitsChrG1,genome2.getId(),chrNameG2);
				alignedUnitsFirstGenome.addAll(selectedUnits);
				List<OrthologyUnit> unitsChrG2 = genome2.getUniqueOrthologyUnits(chrNameG2);
				log.info("Sequence "+chrG1.getName()+" in first genome aligned to sequence "+chrNameG2+" in the second genome. Orthology units sequence genome 1 "+unitsChrG1.size()+". Orthology units sequence genome 2: "+unitsChrG2.size()+" LCS size: "+selectedUnits.size());
				
			} else {
				log.info("Mate sequence not found for "+chrG1.getName()+" Sequence orthology units: "+unitsChrG1.size());
			}
		}
	}
	
	/**
	 * Selects the chromosome having the largest number of mates with the given units
	 * @param units to select the chromosome with the best fit
	 * @param genomeId Id of the genome to query orthologs
	 * @return String name of the most frequent chromosome in the mates of the given units
	 */
	private String findBestChromosome(List<OrthologyUnit> units, int genomeId) {
		Map<String,Integer> chrMateCounts = new HashMap<>();

		//Go over the orthology units. Locate chromosome of each unique ortholog and update counts map
		
		for(OrthologyUnit unit:units)
		{
			OrthologyUnit mate = unit.getUniqueOrtholog(genomeId);
			if(mate == null) continue;
			String sequenceName = mate.getSequenceName();			
			if((chrMateCounts.containsKey(sequenceName)))
			{
				int value = chrMateCounts.get(sequenceName)+1;
				chrMateCounts.put(sequenceName, value);
			}
			else
			{
				chrMateCounts.put(sequenceName, 1);
			}
				
		}
		//Find the chromosome with the largest count
		int bestChromosomeCounts = 0;
		String bestChromosome = null;
		
		for(Map.Entry<String,Integer> entry : chrMateCounts.entrySet()) {
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
	 * @param genome1Id Id of the first genome
	 * @param unitsChrG1 This list has only one chromosome and is sorted by position
	 * @param genome2Id Id of the second genome
	 * @param unitsChrG2 This list has only one chromosome and is sorted by position
	 * @return List<OrthologyUnit> List of selected units of the first list making the LCS relative to the second list
	 */
	private List<OrthologyUnit> alignOrthologyUnits(int genome1Id, List<OrthologyUnit> unitsChrG1, int genome2Id, String seqName2) {
		List<OrthologyUnit> answer = new ArrayList<>();
		
		List<OrthologyUnit> unitsG1List = new ArrayList<>();
		List<OrthologyUnit> unitsG2List = new ArrayList<>();
		//Assigns mates of units in G2
		
		// Select orthology units in G1 having mate in g2 and assign the mate of units in G2.
	
		//At the same time create two new lists of the same size with the  units in g1 having its mate in g2
		
		for(int i=0; i<unitsChrG1.size(); i++)
		{
			OrthologyUnit unitG1 = unitsChrG1.get(i);
			OrthologyUnit unitG2 = unitG1.getUniqueOrtholog(genome2Id);
			if(unitG2==null) continue;
			if(!unitG2.getSequenceName().equals(seqName2)) continue;
			unitsG1List.add(unitG1);			
			unitsG2List.add(unitG2);
		}
		Collections.sort(unitsG2List, GenomicRegionPositionComparator.getInstance());		
		
		
		Set<Integer> lcsForward = findLCS(unitsG1List, unitsG2List, genome2Id);
		
		Collections.reverse(unitsG2List);
		
		Set<Integer> lcsReverse = findLCS(unitsG1List, unitsG2List, genome2Id);
		
		Set<Integer> lcs = lcsForward;
		if(lcsReverse.size()>lcsForward.size()) lcs = lcsReverse;
		//System.out.println("Positions for LCS: "+positions.length+" LCS: "+lcs.size());
		// Select the orthology units in G1 located at the indexes given by the output of LCS
		for(int i:lcs)
		{
			OrthologyUnit lcsResult = unitsG1List.get(i);
			lcsResult.setInLCS(genome2Id);
			OrthologyUnit unitG2 = lcsResult.getUniqueOrtholog(genome2Id);
			unitG2.setInLCS(genome1Id);
			answer.add(lcsResult);
		}
		
		
		return answer;
	}

	private Set<Integer> findLCS(List<OrthologyUnit> unitsG1List, List<OrthologyUnit> unitsG2List, int genome2Id) {
		//Reverse map with unit ids as keys and positions as values
		Map<String,Integer> sortedUnitsG2Pos = new HashMap<>();
		for(int i=0; i<unitsG2List.size(); i++) {
			OrthologyUnit unit=unitsG2List.get(i);
			sortedUnitsG2Pos.put(unit.getId(), i);
		}
		// Create int array with the positions in g2 for the units in g1. Input for LCS
		int []positions = new int[unitsG1List.size()];
		for(int i=0; i<unitsG1List.size(); i++)
		{
			OrthologyUnit unitG1 = unitsG1List.get(i);
			int j = sortedUnitsG2Pos.get(unitG1.getUniqueOrtholog(genome2Id).getId());
			positions[i] = j;
			//System.out.println("Positions [ "+i+"]:"+j);
		}
		// Run LCS
		
		Set<Integer> lcs = findLCS(positions);
		return lcs;
	}

	/**
	 * Calculates the longest common subsequence (LCS) of sorted entries in the given indexes array using dynamic programming 
	 * @param indexesMap Indexes to find the LCS
	 * @return SortedSet<Integer> Positions making the LCS
	 */
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
	
	public void printAlignmentResults() throws IOException {
		AnnotatedReferenceGenome genome1 = genomes.get(0);
		AnnotatedReferenceGenome genome2 = genomes.get(1);
		
		// Print Unique units genome 1 outUniqueG1
		try (PrintStream outUniqueG1 = new PrintStream(outPrefix+"_uniqueG1.tsv");) {
			outUniqueG1.println("geneIdG1\tchromosomeG1\tgeneStartG1\tgeneEndG1\tgeneIdG2\tchromosomeG2\tgeneStartG2\tgeneEndG2\ttype");
			for(OrthologyUnit unit:genome1.getUniqueOrthologyUnits()) {
				printOrthologyUnit(unit, genome2.getId(), outUniqueG1);
			}
		}
		try (PrintStream outUniqueG2 = new PrintStream(outPrefix+"_uniqueG2.tsv");) {
			for(OrthologyUnit unit:genome2.getUniqueOrthologyUnits()) {
				printOrthologyUnit(unit, genome1.getId(), outUniqueG2);
			}
		}
		
		// Print paralogs
		
		try (PrintStream outParalogsG1 = new PrintStream(outPrefix+"_paralogsG1.tsv");) {
			outParalogsG1.println("geneId\tchromosome\tgeneStart\tgeneEnd\tparalogId\tparalogChr\tparalogStart\tparalogEnd");
			for(OrthologyUnit unit:genome1.getOrthologyUnits()) {
				for(OrthologyUnit paralog: unit.getParalogs()) {
					outParalogsG1.print(unit.getId()+"\t"+unit.getSequenceName()+"\t"+unit.getFirst()+"\t"+unit.getLast());
					outParalogsG1.println("\t"+paralog.getId()+"\t"+paralog.getSequenceName()+"\t"+paralog.getFirst()+"\t"+paralog.getLast());
				}
			}
		}
		try (PrintStream outParalogsG2 = new PrintStream(outPrefix+"_paralogsG2.tsv");) {
			outParalogsG2.println("geneId\tchromosome\tgeneStart\tgeneEnd\tparalogId\tparalogChr\tparalogStart\tparalogEnd");
			for(OrthologyUnit unit:genome2.getOrthologyUnits()) {
				for(OrthologyUnit paralog: unit.getParalogs()) {
					outParalogsG2.print(unit.getId()+"\t"+unit.getSequenceName()+"\t"+unit.getFirst()+"\t"+unit.getLast());
					outParalogsG2.println("\t"+paralog.getId()+"\t"+paralog.getSequenceName()+"\t"+paralog.getFirst()+"\t"+paralog.getLast());
				}
			}
		}
		//Print LCS
		//try (PrintStream outAlignmnent = new PrintStream(outPrefix+"_lcs.tsv");) {
			//outAlignmnent.println("geneIdG1\tchromosomeG1\tgeneStartG1\tgeneEndG1\tgeneIdG2\tchromosomeG2\tgeneStartG2\tgeneEndG2");
			//for(OrthologyUnit unit:alignedUnitsFirstGenome) {
				//printOrthologyUnit(unit, genome2.getId(), outAlignmnent);
		//	}
	//	}
		//Print metadata genome 1
		try (PrintStream outGenome1 = new PrintStream(outPrefix+"_genome1.tsv");) {
			printGenomeMetadata(outGenome1, genome1.getSequencesMetadata());
		}
		
		//Print metadata genome 2
		try (PrintStream outGenome2 = new PrintStream(outPrefix+"_genome2.tsv");) {
			printGenomeMetadata(outGenome2, genome2.getSequencesMetadata());
		}
		
		//Print D3 linear visualization
		try (PrintStream outD3Linear = new PrintStream(outPrefix+"_linearView.html");) {
			printD3Visualization(outD3Linear,"GenomesAlignerLinearVisualizer.js");
		}
		
	}

	

	private void printOrthologyUnit(OrthologyUnit unit, int genomeId, PrintStream out) {
		List<OrthologyUnit> orthologs = unit.getOrthologs(genomeId);
		char type = 'U';
		if(orthologs.size()>1) type = 'M';
		else if (unit.isInLCS(genomeId)) type = 'L';
		for(OrthologyUnit ortholog:orthologs) {
			out.print(unit.getId()+"\t"+unit.getSequenceName()+"\t"+unit.getFirst()+"\t"+unit.getLast());
			out.println("\t"+ortholog.getId()+"\t"+ortholog.getSequenceName()+"\t"+ortholog.getFirst()+"\t"+ortholog.getLast()+"\t"+type);
		}
		
	}
	private void printGenomeMetadata(PrintStream out, QualifiedSequenceList sequencesMetadata) {
		out.println("Name\tLength");
		for(QualifiedSequence seq:sequencesMetadata) {
			out.println(""+seq.getName()+"\t"+seq.getLength());
		}
	}
	private void printD3Visualization(PrintStream outD3Linear, String jsFile) throws IOException {
		outD3Linear.println("<!DOCTYPE html>");
		outD3Linear.println("<meta charset=\"utf-8\">");
		outD3Linear.println("<head>");
		outD3Linear.println("<FORM>\n" + 
				"<h1> Genomes Aligner v1.0 &emsp;&emsp;\n" + 
				"<INPUT TYPE=\"button\" onClick=\"history.go(0)\" VALUE=\"Start Again!\">\n" + 
				"</h1>\n" + 
				"</FORM>");
		outD3Linear.println("</titled>");
		outD3Linear.println(htmlStyleCode());
		outD3Linear.println("<body>");
		//Adds the buttons for lcs, multiple and uniques
		outD3Linear.println("<div id=\"option\">\n" + 
				"    <input id=\"LCSbutton\" \n" + 
				"           type=\"button\" \n" + 
				"           value=\"LCS\"/>\n" + 
				"    <input id=\"Multiplebutton\" \n" + 
				"           type=\"button\" \n" + 
				"           value=\"Multiple\"/>\n" + 
				"    <input id=\"Uniquesbutton\" \n" + 
				"           type=\"button\" \n" + 
				"           value=\"Uniques\"/>\n" + 
				"</div>");
		outD3Linear.println("<script src=\"http://d3js.org/d3.v3.min.js\"></script>");
		outD3Linear.println("<script>");
		//Print D3 script
		Class<? extends GenomesAligner> c = this.getClass();
		String resource = "/ngsep/assembly/"+jsFile;
		System.out.println("Loading resource: "+resource);
		try (InputStream is = c.getResourceAsStream(resource);
			 BufferedReader in = new BufferedReader(new InputStreamReader(is))) {
			String line=in.readLine();
			while(line!=null) {
				if(line.contains("InputFileLCS.tsv")) {
					outD3Linear.println("d3.tsv(\""+outPrefix+"_uniqueG1.tsv\", function(error, lcs)");
				} else if(line.contains("InputFileGenome1.tsv")) {
					outD3Linear.println("  d3.tsv(\""+outPrefix+"_genome1.tsv\", function(error, crmsA)");
				} else if(line.contains("InputFileGenome2.tsv")) {
					outD3Linear.println("    d3.tsv(\""+outPrefix+"_genome2.tsv\", function(error, crmsB)");
				} else {
					outD3Linear.println(line);
				}
				line=in.readLine();
			}
		}
		outD3Linear.println("</script>");
		outD3Linear.println("</body>");
	}
	
	private String htmlStyleCode() {
		String code = "<style>\n" + 
				"		svg {\n" + 
				"		  font: 18px sans-serif;\n" + 
				"		}\n" + 
				"		.background path {\n" + 
				"		  fill: none;\n" + 
				"		  stroke: #ddd;\n" + 
				"		  shape-rendering: crispEdges;\n" + 
				"		}\n" + 
				"\n" + 
				"\n" + 
				"\n" + 
				"		.brush .extent {\n" + 
				"		  fill-opacity: .3;\n" + 
				"		  stroke: #fff;\n" + 
				"		  shape-rendering: crispEdges;\n" + 
				"		}\n" + 
				"\n" + 
				"		.axis line,\n" + 
				"		.axis path {\n" + 
				"		  fill: none;\n" + 
				"		  stroke: #000;\n" + 
				"		  shape-rendering: crispEdges;\n" + 
				"		}\n" + 
				"\n" + 
				"		.axis text {\n" + 
				"		  text-shadow: 0 1px 0 #fff, 1px 0 0 #fff, 0 -1px 0 #fff, -1px 0 0 #fff;\n" + 
				"		  cursor: move;\n" + 
				"		}\n" + 
				"\n" + 
				"		</style>";
		
		return code;
	}
}


class AnnotatedReferenceGenome {
	private static final int DEF_PCT_KMERS = 50;
	
	private int id;
	private ProteinTranslator translator = new ProteinTranslator();
	private ReferenceGenome genome;
	private Transcriptome transcriptome;
	private QualifiedSequenceList sequencesMetadata;
	private List<OrthologyUnit> orthologyUnitsList;
	private GenomicRegionSortedCollection<OrthologyUnit> orthologyUnitsBySequence;
	private GenomicRegionSortedCollection<OrthologyUnit> uniqueOrthologyUnitsBySequence;
	private Map<String, OrthologyUnit> orthologyUnitsMap;
	private FMIndex indexOrthologyUnits;
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
		findParalogs();
		selectUniqueOrthologyUnits();
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
		unit.setProteinSequence(bestProtein);
		return unit;
	}
	
	private void buildFMIndex() {
		indexOrthologyUnits = new FMIndex();
		QualifiedSequenceList proteinSequences = new QualifiedSequenceList();
		for (OrthologyUnit ql:orthologyUnitsList) {
			String proteinSequence = ql.getProteinSequence();
			String proteinId = ql.getId();
			QualifiedSequence qualifiedSequence = new QualifiedSequence(proteinId, proteinSequence);
			proteinSequences.add(qualifiedSequence);
		}
		indexOrthologyUnits.loadQualifiedSequenceList(proteinSequences);
	}
	
	private void findParalogs() {
		for (OrthologyUnit unit:orthologyUnitsList) {
			//Orthology unit ids of similar proteins
			Set<String> hits = findSimilarProteins(indexOrthologyUnits, unit.getProteinSequence());
			for(String paralogId:hits) {
				if(!paralogId.equals(unit.getId())) {
					OrthologyUnit paralog = orthologyUnitsMap.get(paralogId);
					unit.addParalog(paralog);
				}
			}
		}
	}
	
	private void selectUniqueOrthologyUnits () {
		uniqueOrthologyUnitsBySequence = new GenomicRegionSortedCollection<>(sequencesMetadata);
		for(OrthologyUnit unit:orthologyUnitsList) {
			if(unit.isUnique()) uniqueOrthologyUnitsBySequence.add(unit);
		}
	}
	
	private Set<String> findSimilarProteins(FMIndex index, String protein1) {	
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

		//Step 2: Fill list traversing the counts and choosing transcripts for which at least x% of the k-mers support the match
		Set<String> answer = new TreeSet<>();
		
		for(Map.Entry<String,Integer> entry : kmerSupportMap.entrySet()) {
			String name = entry.getKey();
			double transcriptKmers = entry.getValue();
			double percent = (transcriptKmers/totalKmers)*100;
			if(percent >= DEF_PCT_KMERS)
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
	public void findOrthologs (AnnotatedReferenceGenome genome2) {
		for (OrthologyUnit unit:orthologyUnitsList) {
			//Orthology unit ids of similar proteins
			Set<String> hits = findSimilarProteins(genome2.indexOrthologyUnits, unit.getProteinSequence());
			for(String orthologId:hits) {
				OrthologyUnit ortholog = genome2.orthologyUnitsMap.get(orthologId);
				unit.addOrtholog(ortholog);
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
	 * @return List<OrthologyUnit> List of unique orthology units in this genome
	 */
	public List<OrthologyUnit> getUniqueOrthologyUnits() {
		return uniqueOrthologyUnitsBySequence.asList();
	}
	
	/**
	 * Returns the list of orthology units for the given sequence name
	 * @param name Sequence name
	 * @return List<OrthologyUnit> units with the given sequence name
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
