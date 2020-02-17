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
package ngsep.genome;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.math.Distribution;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;

/**
 * @author Daniel Tello
 * @author Jorge Duitama
 */
public class GenomesAligner {

	// Constants for default values
	public static final String DEF_OUT_PREFIX = "genomesAlignment";
	public static final byte DEF_KMER_LENGTH = 10;
	public static final int DEF_MIN_PCT_KMERS = 50;
	public static final int DEF_MAX_HOMOLOGS_UNIT = 3;

	// Logging and progress
	private Logger log = Logger.getLogger(GenomesAligner.class.getName());
	private ProgressNotifier progressNotifier=null;

	// Parameters
	private List<AnnotatedReferenceGenome> genomes = new ArrayList<>();
	private String outputPrefix = DEF_OUT_PREFIX;
	private String jsFilename = outputPrefix + "_vizVariables.js";
	private byte kmerLength = DEF_KMER_LENGTH;
	private int minPctKmers = DEF_MIN_PCT_KMERS;
	private int maxHomologsUnit = DEF_MAX_HOMOLOGS_UNIT;

	private List<List<OrthologyUnit>> orthologyUnitClusters=new ArrayList<>();

	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}

	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}

	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
		jsFilename = outputPrefix + "_vizVariables.js";
	}

	public byte getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(byte kmerLength) {
		this.kmerLength = kmerLength;
	}
	public void setKmerLength(String value) {
		setKmerLength((byte)OptionValuesDecoder.decode(value, Byte.class));
	}

	public int getMinPctKmers() {
		return minPctKmers;
	}
	public void setMinPctKmers(int minPctKmers) {
		this.minPctKmers = minPctKmers;
	}
	public void setMinPctKmers(String value) {
		setMinPctKmers((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getMaxHomologsUnit() {
		return maxHomologsUnit;
	}
	public void setMaxHomologsUnit(int maxHomologsUnit) {
		this.maxHomologsUnit = maxHomologsUnit;
	}
	public void setMaxHomologsUnit(String value) {
		setMaxHomologsUnit((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public static void main(String[] args) throws Exception 
	{
		GenomesAligner instance = new GenomesAligner();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		while(i<args.length-1) {
			String fileGenome = args[i++];
			String fileTranscriptome = args[i++];
			instance.loadGenome(fileGenome, fileTranscriptome);
		}
		instance.run();
	}
	
	public void run () throws IOException {
		logParameters ();
		if(genomes.size()==0) throw new IOException("At least one genome and its annotation should be provided");
		if(outputPrefix==null) throw new IOException("A prefix for output files is required");
		alignGenomes();
		printAlignmentResults ();
		log.info("Process finished");
	}
	

	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Loaded: "+ genomes.size()+" annotated genomes");
		out.println("Output prefix:"+ outputPrefix);
		out.println("K-mer length: "+ getKmerLength());
		out.println("Minimum percentage of k-mers to call orthologs: "+ getMinPctKmers());
		log.info(os.toString());
	}
	public void loadGenome(String fileGenome, String fileTranscriptome) throws IOException {
		ReferenceGenome genome = new ReferenceGenome(fileGenome);
		log.info("Loaded genome "+fileGenome);
		GFF3TranscriptomeHandler transcriptomeHandler = new GFF3TranscriptomeHandler(genome.getSequencesMetadata());
		Transcriptome transcriptome = transcriptomeHandler.loadMap(fileTranscriptome);
		transcriptome.fillSequenceTranscripts(genome, log);
		log.info("Loaded transcriptome "+fileTranscriptome+ " number of transcripts: "+transcriptome.getAllTranscripts().size());
		AnnotatedReferenceGenome annGenome = new AnnotatedReferenceGenome(genomes.size()+1, genome, transcriptome);
		log.info("Genome: "+annGenome.getId()+" has "+annGenome.getOrthologyUnits().size()+" total orthology units. Calculating Paralogs");
		annGenome.calculateParalogs(kmerLength, minPctKmers);
		log.info("Paralogs found for Genome: "+annGenome.getId()+" Unique orthology units: "+annGenome.getUniqueOrthologyUnits().size());
		genomes.add(annGenome);
	}


	public void alignGenomes() {

		for(int i=0;i<genomes.size();i++) {
			for (int j=0;j<genomes.size();j++) {
				if(i!=j) genomes.get(i).calculateOrthologs(genomes.get(j), kmerLength, minPctKmers);
			}
		}
		calculateOrthologClusters();
		if(genomes.size()<2) return;
		// By now this is still done for two genomes
		AnnotatedReferenceGenome genome1 = genomes.get(0);
		AnnotatedReferenceGenome genome2 = genomes.get(1);
		QualifiedSequenceList sequencesG1 = genome1.getSequencesMetadata();
		for(QualifiedSequence chrG1:sequencesG1) {
			List<OrthologyUnit> unitsChrG1 = genome1.getUniqueOrthologyUnits(chrG1.getName());
			log.info("Unique units G1 for "+chrG1.getName()+": "+unitsChrG1.size());
			String chrNameG2 = findBestChromosome(unitsChrG1, genome2.getId());
			if(chrNameG2!=null) {
				List<OrthologyUnit> selectedUnits = alignOrthologyUnits(genome1.getId(),unitsChrG1,genome2.getId(),chrNameG2);

				List<OrthologyUnit> unitsChrG2 = genome2.getUniqueOrthologyUnits(chrNameG2);
				log.info("Sequence "+chrG1.getName()+" in first genome aligned to sequence "+chrNameG2+" in the second genome. Orthology units sequence genome 1 "+unitsChrG1.size()+". Orthology units sequence genome 2: "+unitsChrG2.size()+" LCS size: "+selectedUnits.size());
				completeLCS(genome1.getId(),genome1.getOrthologyUnits(chrG1.getName()),genome2.getId());
			} else {
				log.info("Mate sequence not found for "+chrG1.getName()+" Sequence orthology units: "+unitsChrG1.size());
			}
		}
	}



	private void calculateOrthologClusters() {
		log.info("Clustering orthologs and paralogs");
		orthologyUnitClusters=new ArrayList<>();
		List<OrthologyUnit> unitsWithOrthologs = new ArrayList<>();
		Map<String, Integer> unitsPositionMap = new HashMap<>();
		for(AnnotatedReferenceGenome genome:genomes) {
			for(OrthologyUnit unit:genome.getOrthologyUnits()) {
				if(unit.getTotalOrthologs()>0) {
					unitsPositionMap.put(unit.getUniqueKey(), unitsWithOrthologs.size());
					unitsWithOrthologs.add(unit);
				}

			}
		}
		log.info("Total units with orthologs: "+unitsWithOrthologs.size());
		int [] group = new int [unitsWithOrthologs.size()];
		Arrays.fill(group, -1);

		Distribution distClusterSizes = new Distribution(0, 50, 1);
		//Build connected components
		for(int i=0;i<group.length;i++) {
			if(group[i]!=-1) continue;
			List<OrthologyUnit> cluster = new ArrayList<>();
			int groupNumber = orthologyUnitClusters.size();
			Queue<OrthologyUnit> agenda = new LinkedList<>();
			agenda.add(unitsWithOrthologs.get(i));
			while(agenda.size()>0) {
				OrthologyUnit unit = agenda.remove();
				int pos = unitsPositionMap.get(unit.getUniqueKey());
				int unitGroup = group[pos];
				if(unitGroup==-1) {
					cluster.add(unit);
					group[pos]=groupNumber;
					agenda.addAll(unit.getOrthologsAllGenomes());
				} else if(unitGroup!=groupNumber) log.warning("Possible connection between clusters "+unitGroup + " and "+groupNumber+" Unit: "+unit.getUniqueKey());

			}
			if(cluster.size()>1) {
				orthologyUnitClusters.add(cluster);
				distClusterSizes.processDatapoint(cluster.size());
			} else if (cluster.size()==0) {
				log.warning("Empty cluster from unit: "+unitsWithOrthologs.get(i).getUniqueKey()+" clusters: "+orthologyUnitClusters.size()+" groupNumber: "+groupNumber);
			}

		}
		log.info("Number of clusters: "+orthologyUnitClusters.size());
		//TODO: Report it better
		distClusterSizes.printDistributionInt(System.out);
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
			OrthologyUnit unitG2 = lcsResult.getUniqueOrtholog(genome2Id);
			lcsResult.setMateInLCS(unitG2);
			unitG2.setMateInLCS(lcsResult);
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

	private void completeLCS(int genomeId1, List<OrthologyUnit> chrUnits, int genomeId2) {
		int i1=-1;
		OrthologyUnit mate1=null;
		for(int i=0;i<chrUnits.size();i++) {
			OrthologyUnit chrUnit = chrUnits.get(i);
			OrthologyUnit mateInLCS = chrUnit.getLCSMate(genomeId2);
			if(mateInLCS==null) continue;

			if(i1==-1) {
				i1 = i;
				mate1 = mateInLCS;
				continue;
			}
			boolean reverse = false;
			String seqName2 = mate1.getSequenceName();
			//In yeast some genes intersect
			int firstG2 = mate1.getFirst()+1;
			int lastG2 = mateInLCS.getLast()-1;
			if(firstG2>lastG2) {
				reverse = true;
				firstG2 = mateInLCS.getFirst()+1;
				lastG2 = mate1.getLast()-1;
			}
			for(int j=i1+1;j<i;j++) {
				OrthologyUnit chrUnitB = chrUnits.get(j);
				Collection<OrthologyUnit> orthologs2 = chrUnitB.getOrthologs(genomeId2);
				if(chrUnitB.getId().equals("YJR023C_EC1118")) System.out.println("Orthologs for "+chrUnitB.getId()+" "+orthologs2.size()+" range G2: "+firstG2+"-"+lastG2);
				if(orthologs2.size()==0) continue;
				List<OrthologyUnit> matesInRange = new ArrayList<>();
				for(OrthologyUnit ortholog2:orthologs2) {
					if(seqName2.equals(ortholog2.getSequenceName()) && ortholog2.getFirst()>=firstG2 && ortholog2.getLast()<=lastG2) {
						matesInRange.add(ortholog2);
					}
				}
				if(chrUnitB.getId().equals("YJR023C_EC1118")) System.out.println("Mates in range for "+chrUnitB.getId()+" "+matesInRange.size());
				if(matesInRange.size()!=1) continue;
				OrthologyUnit mate = matesInRange.get(0);
				chrUnitB.setMateInLCS(mate);
				mate.setMateInLCS(chrUnitB);
				if(reverse) {
					lastG2 = mate.getLast()-1;
				} else {
					firstG2 = mate.getFirst()+1; 
				}
			}
			i1=i;
			mate1 = mateInLCS;
		}
	}

	public void printAlignmentResults() throws IOException {
		 // Create Javascript file for visualization variables
        try (PrintStream outJS = new PrintStream(jsFilename);) {
        	outJS.print("");
        }

		for(int i=0;i<genomes.size();i++) {
			AnnotatedReferenceGenome genome = genomes.get(i);
			int id = genome.getId();
			//Print metadata
			printGenomeMetadata(id, genome.getSequencesMetadata());
			// Print paralogs
			try (PrintStream outParalogs = new PrintStream(outputPrefix+"_paralogsG"+id+".tsv");) {
				outParalogs.println("geneId\tchromosome\tgeneStart\tgeneEnd\tparalogsCount\tgenomeId\tparalogId\tparalogChr\tparalogStart\tparalogEnd");
				for(OrthologyUnit unit:genome.getOrthologyUnits()) {
					Collection<OrthologyUnit> paralogs = unit.getParalogs(); 
					for(OrthologyUnit paralog: paralogs) {
						outParalogs.print(unit.getId()+"\t"+unit.getSequenceName()+"\t"+unit.getFirst()+"\t"+unit.getLast()+"\t"+paralogs.size());
						outParalogs.println("\t"+id+"\t"+paralog.getId()+"\t"+paralog.getSequenceName()+"\t"+paralog.getFirst()+"\t"+paralog.getLast());
					}
				}
			}
			try (PrintStream outParalogsJS = new PrintStream(new FileOutputStream(jsFilename, true));) {
				outParalogsJS.println("const paralogsG" + id + " = [");
				for(OrthologyUnit unit:genome.getOrthologyUnits()) {
					for(OrthologyUnit paralog: unit.getParalogs()) {
						outParalogsJS.println("{geneId: '"+unit.getId()
						+"', chromosome: '"+unit.getSequenceName()
						+"', geneStart: "+unit.getFirst()
						+", geneEnd: "+unit.getLast()
						+", paralogId: '"+paralog.getId()
						+"', paralogChr: '"+paralog.getSequenceName()
						+"', paralogStart: "+paralog.getFirst()
						+", paralogEnd: "+paralog.getLast()+"},");
					}
				}
				outParalogsJS.println("];");
			}
		}

		//Print ortholog clusters
		try (PrintStream outClusters = new PrintStream(outputPrefix+"_clusters.txt");) {
			for(List<OrthologyUnit> cluster:orthologyUnitClusters) {
				outClusters.print(cluster.get(0).getId());
				for(int i=1;i<cluster.size();i++) {
					OrthologyUnit unit = cluster.get(i);
					outClusters.print("\t"+unit.getId());
				}
				outClusters.println();
			}
		}

		try (PrintStream outD3Paralogs = new PrintStream(outputPrefix+"_circularParalogView.html");) {
			printD3Visualization(outD3Paralogs,"GenomesAlignerCircularParalogVisualizer.js", 5);
		}
		
		if(genomes.size()>1) {
			for(int i=0;i<genomes.size();i++) {
				AnnotatedReferenceGenome genome = genomes.get(i);
				int id = genome.getId();
				// Print orthologs
				try (PrintStream outOrthologs = new PrintStream(outputPrefix+"_orthologsG"+id+".tsv");) {
					try (PrintStream outOrthologsJS = new PrintStream(new FileOutputStream(jsFilename, true));) {
						outOrthologs.println("geneId\tchromosome\tgeneStart\tgeneEnd\tparalogsCount\tgenomeId2\tgeneIdG2\tchromosomeG2\tgeneStartG2\tgeneEndG2\ttype");
						outOrthologsJS.println("const orthologsG" + id + " = [");
						for(OrthologyUnit unit:genome.getOrthologyUnits()) {
							printOrthologyUnit(unit, outOrthologs, outOrthologsJS);
						}
						outOrthologsJS.println("];");
					}
				}
			}
			//Print D3 visualizations
			try (PrintStream outD3Linear = new PrintStream(outputPrefix+"_linearOrthologView.html");) {
				printD3Visualization(outD3Linear,"GenomesAlignerLinearOrthologVisualizer.js", 5);
			}
			
			try (PrintStream outD3Circular = new PrintStream(outputPrefix+"_circularOrthologView.html");) {
				printD3Visualization(outD3Circular,"GenomesAlignerCircularOrthologVisualizer.js", 5);
			}
		}



	}
	private void printGenomeMetadata(int id, QualifiedSequenceList sequencesMetadata) throws IOException {
		String outFilename = outputPrefix+"_genome"+id+".tsv";
		try (PrintStream out = new PrintStream(outFilename)) {
			try (PrintStream outJS = new PrintStream(new FileOutputStream(jsFilename, true))) {
				outJS.println("const genome" + id + " = [");
				out.println("Name\tLength");
				for(QualifiedSequence seq:sequencesMetadata) {
					out.println(""+seq.getName()+"\t"+seq.getLength());
					outJS.println("{Name: '"+seq.getName()+"', Length: "+seq.getLength()+"},");
				}
				outJS.println("];");
			}
		}
	}


	private void printOrthologyUnit(OrthologyUnit unit, PrintStream out, PrintStream outJS) {
		List<OrthologyUnit> orthologs = unit.getOrthologsOtherGenomes();
		if(orthologs.size()==0) return;
		char type = 'U';
		if(orthologs.size()>1) type = 'M';
		OrthologyUnit mateInLCS = unit.getLCSMate(orthologs.get(0).getGenomeId());
		for(OrthologyUnit ortholog:orthologs) {
			out.print(unit.getId()+"\t"+unit.getSequenceName()+"\t"+unit.getFirst()+"\t"+unit.getLast()+"\t"+unit.getParalogs().size());
			char typePrint = type;
			if (ortholog==mateInLCS) typePrint = 'L';
			out.println("\t"+ortholog.getGenomeId()+"\t"+ortholog.getId()+"\t"+ortholog.getSequenceName()+"\t"+ortholog.getFirst()+"\t"+ortholog.getLast()+"\t"+typePrint);
			outJS.println("{geneId: '"+unit.getId()
			+"', chromosome: '"+unit.getSequenceName()
			+"', geneStart: "+unit.getFirst()
			+", geneEnd: "+unit.getLast()
			+", genomeId2: '"+ortholog.getGenomeId()
			+"', geneIdG2: '"+ortholog.getId()
			+"', chromosomeG2: '"+ortholog.getSequenceName()
			+"', geneStartG2: "+ortholog.getFirst()
			+", geneEndG2: "+ortholog.getLast()
			+", type: '"+typePrint+"'},");
		}
	}

	private void printD3Visualization(PrintStream outD3, String jsFile, int preferredD3Version) throws IOException {
		outD3.println("<!DOCTYPE html>");
		outD3.println("<meta charset=\"utf-8\">");
		outD3.println("<head>");
		outD3.println("</head>");
		outD3.println(htmlStyleCode());
		outD3.println("<body>");
		//Adds the buttons for lcs, multiple and uniques
		outD3.println("<div id=\"option\"></div>\n");
		if (preferredD3Version == 3) {
			outD3.println("<script src=\"http://d3js.org/d3.v3.min.js\"></script>");
		}
		else if(preferredD3Version == 4) {
			outD3.println("<script src=\"http://d3js.org/d3.v4.min.js\"></script>");
		}
		else if(preferredD3Version == 5){
			outD3.println("<script src=\"http://d3js.org/d3.v5.min.js\"></script>");
		}
		outD3.println("<script src="+jsFilename+"></script>");
		outD3.println("<script>");
		//Print D3 script
		Class<? extends GenomesAligner> c = this.getClass();
		String resource = "/ngsep/genome/"+jsFile;
		log.info("Loading resource: "+resource);
		try (InputStream is = c.getResourceAsStream(resource);
				BufferedReader in = new BufferedReader(new InputStreamReader(is))) {
			String line=in.readLine();
			outD3.println("const MAX_HOMOLOGS_UNIT = "+maxHomologsUnit+";");

			while(line!=null) {
				outD3.println(line);
				line=in.readLine();
			}
		}
		outD3.println("</script>");
		outD3.println("</body>");
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
