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
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;

/**
 * @author Daniel Tello
 * @author Laura Gonzalez
 * @author Jorge Duitama
 */
public class GenomesAligner {

	// Constants for default values
	public static final String DEF_OUT_PREFIX = "genomesAlignment";
	public static final byte DEF_KMER_LENGTH = HomologRelationshipsFinder.DEF_KMER_LENGTH;
	public static final int DEF_MIN_PCT_KMERS = HomologRelationshipsFinder.DEF_MIN_PCT_KMERS;
	public static final int DEF_MAX_HOMOLOGS_UNIT = 3;
	public static final int DEF_MIN_HOMOLOGY_UNITS_BLOCK = PairwiseSyntenyBlocksFinder.DEF_MIN_HOMOLOGY_UNITS_BLOCK;
	public static final int DEF_MIN_BLOCK_LENGTH = PairwiseSyntenyBlocksFinder.DEF_MIN_BLOCK_LENGTH;
	public static final int DEF_MAX_DISTANCE_BETWEEN_UNITS = PairwiseSyntenyBlocksFinder.DEF_MAX_DISTANCE_BETWEEN_UNITS;
	public static final double DEF_MIN_FREQUENCY_SOFT_CORE = 0.9;

	// Logging and progress
	private Logger log = Logger.getLogger(GenomesAligner.class.getName());
	private ProgressNotifier progressNotifier=null;

	// Parameters
	private List<AnnotatedReferenceGenome> genomes = new ArrayList<>();
	private String outputPrefix = DEF_OUT_PREFIX;
	private int maxHomologsUnit = DEF_MAX_HOMOLOGS_UNIT;
	private boolean skipMCL= false;
	private int minHomologUnitsBlock = DEF_MIN_HOMOLOGY_UNITS_BLOCK;
	private int maxDistanceBetweenUnits = DEF_MAX_DISTANCE_BETWEEN_UNITS;
	private double minFrequencySoftCore = DEF_MIN_FREQUENCY_SOFT_CORE;
	private String inputFile = null;
	private String inputDirectory = null;

	
	// Model attributes
	private HomologRelationshipsFinder homologRelationshipsFinder = new HomologRelationshipsFinder();
	//private List<HomologyEdge> homologyEdges = new ArrayList<HomologyEdge>();
	
	private List<HomologyCluster> homologyClusters = new ArrayList<>();
	private List<PairwiseSyntenyBlock> orthologsSyntenyBlocks = new ArrayList<>();
	private int[][] paMatrix;

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
	}
	
	public boolean getSkipMCL() {
		return skipMCL;
	}
	public void setSkipMCL(boolean skipMCL) {
		this.skipMCL = skipMCL;
	}
	public void setSkipMCL(Boolean value) {
		setSkipMCL(value.booleanValue());
	}

	public byte getKmerLength() {
		return homologRelationshipsFinder.getKmerLength();
	}
	public void setKmerLength(byte kmerLength) {
		homologRelationshipsFinder.setKmerLength(kmerLength);
	}
	public void setKmerLength(String value) {
		setKmerLength((byte)OptionValuesDecoder.decode(value, Byte.class));
	}

	public int getMinPctKmers() {
		return homologRelationshipsFinder.getMinPctKmers();
	}
	public void setMinPctKmers(int minPctKmers) {
		homologRelationshipsFinder.setMinPctKmers(minPctKmers);
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
	
	public double getMinFrequencySoftCore() {
		return minFrequencySoftCore;
	}
	public void setMinFrequencySoftCore(double minFrequencySoftCore) {
		this.minFrequencySoftCore = minFrequencySoftCore;
	}
	public void setMinFrequencySoftCore(String value) {
		setMinFrequencySoftCore((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	public String getInputDirectory() {
		return inputDirectory;
	}
	public void setInputDirectory(String inputDirectory) {
		this.inputDirectory = inputDirectory;
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
		
		if(getInputFile()!= null) loadGenomesFromFile();
		if(genomes.size()==0) throw new IOException("At least one genome and its annotation should be provided");
		if(outputPrefix==null) throw new IOException("A prefix for output files is required");
		inferOrthologs();
		identifyHomologyClusters();
		if(genomes.size()>1) {
			alignGenomes();
			buildPAMatrix();
		}
		printResults();
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
		log.info("Genome: "+annGenome.getId()+" has "+annGenome.getTotalHomologyUnits()+" total homology units.");
		genomes.add(annGenome);
	}

	public void loadGenomesFromFile() throws IOException {
		if(inputDirectory == null) throw new IOException("You must specify the input folder");
		if(outputPrefix==null) throw new IOException("A prefix for output files is required");
		String fileSeparator = File.separator;
		try (FileReader reader = new FileReader(inputDirectory + fileSeparator + inputFile);
			 BufferedReader in = new BufferedReader(reader)) {
			String line = in.readLine();
			log.info("Loading genomes in " + inputDirectory + fileSeparator + inputFile);
			while(line!=null) {
				log.info("Loading genome " + line);
				String pathFile = inputDirectory + fileSeparator + line;
				String fastaFile = checkExistingFasta(pathFile);
				String gffFile = checkExistingGFF(pathFile);
				if(fastaFile == null || gffFile == null)
					log.warning("FASTA or GFF3 for the genome " + line + " did not exist. The genome was removed from the analysis.");
				else loadGenome(fastaFile,gffFile);
				line = in.readLine();
			}
		}
	}
	
	/**
	 * Check if a FASTA file exists given a prefix
	 * @param pathFile
	 * @return
	 */
	public String checkExistingFasta(String pathFile) 
	{
		List<String> extensions = new ArrayList<>();
		extensions.add(".fna");
		extensions.add(".fna.gz");
		extensions.add(".fa");
		extensions.add(".fa.gz");
		extensions.add(".fas");
		extensions.add(".fas.gz");
		extensions.add(".fasta");
		extensions.add(".fasta.gz");
		for (String string : extensions) {
			String pathFastaString = pathFile+string;
			File file = new File(pathFastaString);
			if(file.exists()) return pathFastaString;
		}	
		return null;
	}
	
	/**
	 * Check if a GFF3 annotation file exists given a prefix
	 * @param pathFile
	 * @return
	 */
	public String checkExistingGFF(String pathFile) 
	{
		List<String> extensions = new ArrayList<>();
		extensions.add(".gff");
		extensions.add(".gff.gz");
		extensions.add(".gff3");
		extensions.add(".gff3.gz");
		for (String string : extensions) {
			String pathGffString = pathFile+string;
			File file = new File(pathGffString);
			if(file.exists()) return pathGffString;
		}	
		return null;
	}
	
	private void inferOrthologs() {
		genomesDescription();
		// Identify paralogs
		for(int i=0;i<genomes.size();i++) {
			AnnotatedReferenceGenome genome = genomes.get(i);
			List<HomologyEdge> edges = homologRelationshipsFinder.calculateParalogs(genome.getHomologyUnits());
			log.info(String.format("Paralogs found for Genome #%d: %d", i+1, edges.size()));
		}
		
		//Identify orthology relationships between pairs of genomes
		for(int i=0;i<genomes.size();i++) {
			AnnotatedReferenceGenome genome1 = genomes.get(i);
			for (int j=0;j<genomes.size();j++) {
				AnnotatedReferenceGenome genome2 = genomes.get(j);
				if(i!=j) {
					List<HomologyEdge> edges = homologRelationshipsFinder.calculateOrthologs(genome1.getHomologyUnits(), genome2.getHomologyUnits());
					log.info(String.format("Orthologs found for Genome #%d #%d: %d", i+1, j+1, edges.size()));
				}
			}
		}
	}
	
	private void genomesDescription() {
		log.info("Total number of genomes: " + genomes.size());
		for(int i = 0; i < genomes.size(); i++) {
			AnnotatedReferenceGenome genome = genomes.get(i);
			log.info(String.format("Genome #%d has %d genes.", i+1, genome.getHomologyUnits().size()));
		}
	}
	
	public void identifyHomologyClusters () {
		HomologClustersCalculator calculator = new HomologClustersCalculator(homologRelationshipsFinder,skipMCL);
		calculator.setLog(log);
		homologyClusters = calculator.clusterHomologs(genomes);
	}
	
	private List<HomologyEdge> getAllHomologyEdges() {
		List<HomologyEdge> answer = new ArrayList<HomologyEdge>();
		for(AnnotatedReferenceGenome genome:genomes) {
			List<HomologyUnit> units = genome.getHomologyUnits();
			for(HomologyUnit unit:units) answer.addAll(unit.getAllHomologyRelationships());
		}
		return answer;
	}
	
	public void alignGenomes() {
		
		//Select finder:
		
		//PairwiseSyntenyBlocksFinder finder = new LCSMainPairwiseSyntenyBlocksFinder();
		//HalSyntenyPairwiseSyntenyBlocksFinder finder = new HalSyntenyPairwiseSyntenyBlocksFinder();
		DAGChainerPairwiseSyntenyBlocksFinder finder = new DAGChainerPairwiseSyntenyBlocksFinder();
		
		//Set parameters according to the user input
		finder.setMaxDistance(maxDistanceBetweenUnits);
		finder.setMinHomologUnitsBlock(minHomologUnitsBlock);
		
		for(int i=0;i<genomes.size();i++) {
			for(int j=i+1;j<genomes.size();j++) {
				AnnotatedReferenceGenome genome1 = genomes.get(i);
				AnnotatedReferenceGenome genome2 = genomes.get(j);
				log.info("Aligning genome "+i+" with genome "+j);
				orthologsSyntenyBlocks = finder.findSyntenyBlocks(genome1, genome2, homologyClusters);
			}
		}
	}
	
	/**
	 * Build P/A matrix based on Homology Units
	 */
	public void buildPAMatrix()
	{
		log.info("Building P/A matrix");
		int numGenomes = genomes.size();
		
		homologyClusters.addAll(getPrivateGeneFamilies());
		
		int numGeneFamilies = homologyClusters.size();
		
		paMatrix = new int[numGeneFamilies][numGenomes];
		
		for(int i=0; i<numGeneFamilies;i++)
		{
			List<HomologyUnit> cluster = homologyClusters.get(i).getHomologyUnitsCluster();
			for(HomologyUnit hom : cluster)
			{
				paMatrix[i][hom.getGenomeId()-1] ++;
			}
		}
		calculateFrequencies(minFrequencySoftCore);
		printPAMatrix();
			
		log.info("Genomes loaded: " + numGenomes);
		log.info("Gene Families loaded: " + numGeneFamilies);
		
	}
	
	/**
	 * Get private gene families to build the cloud/unique genome
	 * @return
	 */
	public List<HomologyCluster> getPrivateGeneFamilies()
	{
		List<HomologyCluster> privateGeneFamilies = new ArrayList<>();
		int count = homologyClusters.size();
	
		//Build a set with genes included in orthologs
		Set<String> genesIncluded = new HashSet<String>();
		for(int i=0; i<homologyClusters.size();i++)
		{
			List<HomologyUnit> cluster = homologyClusters.get(i).getHomologyUnitsCluster();
			for(HomologyUnit hom: cluster)
				genesIncluded.add(hom.getId());
		}
		
		for(AnnotatedReferenceGenome genome : genomes)
		{
			List<HomologyUnit> allGenes = genome.getHomologyUnits();
			
			for(HomologyUnit hom2 : allGenes)
			{
				if(!genesIncluded.contains(hom2.getId())) 
				{
					List<HomologyUnit> homologsList= new ArrayList<>();
					homologsList.add(hom2);
					
					HomologyCluster homclus = new HomologyCluster(count,homologsList);
					privateGeneFamilies.add(homclus);
					genesIncluded.add(hom2.getId());
					//log.info("Private size: " + privateGeneFamilies.size() + " Count value: " + count);
					count ++;
				}
			}
		}
		log.info("Detected "+privateGeneFamilies.size()+" private genes. Total clusters: " + count);
		return  privateGeneFamilies;
	}
	
	/**
	 * Calculate the frequency of each gene family within the genomes
	 * Categorize each family according to exact and soft threshold
	 */
	public void calculateFrequencies(double freqSoft)
	{
		for(int i=0;i<homologyClusters.size();i++)
		{
			double countFreq = 0;
			int totGenomes = paMatrix[i].length;
			for(int j=0;j<totGenomes;j++)
			{
				if(paMatrix[i][j]>0) countFreq ++;
			}
			double freq = countFreq/totGenomes;
			
			HomologyCluster cluster = homologyClusters.get(i);
			
			cluster.setFrequency(freq);
			String exact = (freq == 1) ? HomologyCluster.ECORE: HomologyCluster.EACCESORY;
			cluster.setExactCategory(exact);
			String soft = (freq >= freqSoft) ? HomologyCluster.SCORE: HomologyCluster.SACCESORY;
			cluster.setSoftCategory(soft);
		}
	}
	


	public void printResults() throws IOException {
		String jsFilename = outputPrefix + "_vizVariables.js";
		 // Create Javascript file for visualization variables
        try (PrintStream outJS = new PrintStream(jsFilename);) {
        	outJS.print("");
        }
        Map<String,Integer> clusterIdsByHomologyUnit = new HashMap<String, Integer>();
        for(HomologyCluster cluster:homologyClusters) {
        	for(HomologyUnit unit:cluster.getHomologyUnitsCluster()) clusterIdsByHomologyUnit.put(unit.getUniqueKey(), cluster.getClusterId());
        }
        Map<String,Integer> syntenyBlockIdByPair = new HashMap<String, Integer>();
        for(int i=0;i<orthologsSyntenyBlocks.size();i++) {
        	PairwiseSyntenyBlock block = orthologsSyntenyBlocks.get(i);
        	for(SyntenyVertex vertex:block.getHomologies()) {
        		List<HomologyUnit> units1 = vertex.getLocalRegion1().getHomologyUnitsCluster();
        		List<HomologyUnit> units2 = vertex.getLocalRegion2().getHomologyUnitsCluster();
        		for(HomologyUnit u1:units1) {
        			for(HomologyUnit u2:units2) {
        				syntenyBlockIdByPair.put(u1.getUniqueKey()+"-"+u2.getUniqueKey(), i+1);
        				syntenyBlockIdByPair.put(u2.getUniqueKey()+"-"+u1.getUniqueKey(), i+1);
        			}
        		}
        	}
        	
        }
		for(int i=0;i<genomes.size();i++) {
			AnnotatedReferenceGenome genome = genomes.get(i);
			int id = genome.getId();
			//Print metadata
			printGenomeMetadata(id, genome.getSequencesMetadata(), jsFilename);
			// Print paralogs
			try (PrintStream outParalogs = new PrintStream(outputPrefix+"_paralogsG"+id+".tsv");
				 PrintStream outParalogsJS = new PrintStream(new FileOutputStream(jsFilename, true))) {
				outParalogs.println("geneId\tchromosome\tgeneStart\tgeneEnd\torthogroup\tgenomeId\tparalogId\tparalogChr\tparalogStart\tparalogEnd\torthogroup2\tscore");
				outParalogsJS.println("const paralogsG" + id + " = [");
				for(HomologyUnit unit:genome.getHomologyUnits()) {
					int clusterId1 = clusterIdsByHomologyUnit.getOrDefault(unit.getUniqueKey(), -1);
					Collection<HomologyEdge> paralogRelationships = unit.getParalogRelationships(); 
					for(HomologyEdge paralogRelationship: paralogRelationships) {
						HomologyUnit paralog = paralogRelationship.getSubjectUnit();
						int clusterId2 = clusterIdsByHomologyUnit.getOrDefault(paralog.getUniqueKey(), -1);
						outParalogs.print(unit.getId()+"\t"+unit.getSequenceName()+"\t"+unit.getFirst()+"\t"+unit.getLast()+"\t"+clusterId1);
						outParalogs.print("\t"+id+"\t"+paralog.getId()+"\t"+paralog.getSequenceName()+"\t"+paralog.getFirst()+"\t"+paralog.getLast()+"\t"+clusterId2);
						outParalogs.println("\t"+paralogRelationship.getScore());
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

		/*try (PrintStream outD3Paralogs = new PrintStream(outputPrefix+"_circularParalogView.html");) {
			printD3Visualization(outD3Paralogs,"GenomesAlignerCircularParalogVisualizer.js", jsFilename, 5);
		}*/
		
		if(genomes.size()>1) {
			for(int i=0;i<genomes.size();i++) {
				AnnotatedReferenceGenome genome = genomes.get(i);
				int id = genome.getId();
				// Print orthologs
				try (PrintStream outOrthologs = new PrintStream(outputPrefix+"_orthologsG"+id+".tsv");
					PrintStream outOrthologsJS = new PrintStream(new FileOutputStream(jsFilename, true));) {
						outOrthologs.println("geneId\tchromosome\tgeneStart\tgeneEnd\torthogroup\tgenomeId2\tgeneIdG2\tchromosomeG2\tgeneStartG2\tgeneEndG2\torthogroup2\tscore\tblock");
						outOrthologsJS.println("const orthologsG" + id + " = [");
						for(HomologyUnit unit:genome.getHomologyUnits()) {
							printOrthologyUnit(unit, outOrthologs, outOrthologsJS,clusterIdsByHomologyUnit,syntenyBlockIdByPair);
						}
						outOrthologsJS.println("];");
				}
			}
			//Print D3 visualizations
			try (PrintStream outD3Linear = new PrintStream(outputPrefix+"_linearOrthologView.html");) {
				printD3Visualization(outD3Linear,"GenomesAlignerLinearOrthologVisualizer.js", jsFilename, 5);
			}
			
			/*try (PrintStream outD3Circular = new PrintStream(outputPrefix+"_circularOrthologView.html");) {
				printD3Visualization(outD3Circular,"GenomesAlignerCircularOrthologVisualizer.js", jsFilename, 5);
			}*/
		}

		//Print ortholog clusters
		CDNACatalogAligner.printResults(outputPrefix, homologyClusters);

		printSyntenyBlocks(outputPrefix+"_syntenyBlocks.txt", jsFilename);
		
	}
	
	/**
	 * Print synteny blocks
	 */
	private void printSyntenyBlocks(String outFilename, String jsFilename) throws IOException {
		try (PrintStream outSynteny = new PrintStream(outFilename);
			 PrintStream outSyntenyJS = new PrintStream(new FileOutputStream(jsFilename, true));){
			String headers = "SequenceName1\tStart1\tEnd1\tSequenceName2\tStart2\tEnd2";
			outSynteny.println(headers);
			outSyntenyJS.println("const syntenyBlocks = [");
			for (PairwiseSyntenyBlock sb : orthologsSyntenyBlocks) {
				GenomicRegion r1 = sb.getRegionGenome1();
				GenomicRegion r2 = sb.getRegionGenome2();
				int length1 = genomes.get(sb.getGenomeId1()-1).getSequenceByName(r1.getSequenceName()).getLength();
				int length2 = genomes.get(sb.getGenomeId2()-1).getSequenceByName(r2.getSequenceName()).getLength();
				char orientation = r2.isNegativeStrand()?'-':'+';
				String line = r1.getSequenceName() + "\t" + length1 + "\t" + r1.getFirst() +  "\t" + r1.getLast();
				line+= "\t"+r2.getSequenceName() + "\t" + length2 + "\t" + r2.getFirst() +  "\t" + r2.getLast()+ "\t"+orientation;
				outSynteny.println(line);
//				Printing of homology units that form the synteny block. 
				
				for (SyntenyVertex vertex : sb.getHomologies()) {
					LocalHomologyCluster c1 = vertex.getLocalRegion1();
					LocalHomologyCluster c2 = vertex.getLocalRegion2();
					line="\t";
					for(HomologyUnit u1:c1.getHomologyUnitsCluster()) line+=u1.getId()+",";
					line+="\t";
					for(HomologyUnit u2:c2.getHomologyUnitsCluster()) line+=u2.getId()+",";
					outSynteny.println(line);
				}
				outSyntenyJS.println("{genomeIdG1: "+sb.getGenomeId1()
				+", chromosomeG1: '"+r1.getSequenceName()+"'"
				+", regionStartG1: "+r1.getFirst()
				+", regionEndG1: "+r1.getLast()
				+", genomeIdG2: "+sb.getGenomeId2()
				+", chromosomeG2: '"+r2.getSequenceName() + "'" 
				+", regionStartG2: "+r2.getFirst()
				+", regionEndG2: "+r2.getLast()
				+", negativeG2: "+(r2.isNegativeStrand()?1:0)
				+"},");
			}
			outSyntenyJS.println("];");
		}
	}
	
	/**
	 * Print PA matrix and table of frequencies
	 */
	private void printPAMatrix()
	{
		//Print PA matrix and frequencies
		try (PrintStream outPA = new PrintStream(outputPrefix+"_paMatrix.txt");PrintStream outFreq = new PrintStream(outputPrefix+"_gfFreqs.txt")) 
		{
			outPA.print("");
			for(AnnotatedReferenceGenome genome : genomes)
			{
				outPA.print("\t"+genome.getId());
			}
			
			outFreq.print("GeneFamily\tFrequency\tExactGroup\tSoftGroup\n");
			for(int i=0;i<paMatrix.length;i++)
			{
				outPA.println();
				outPA.print("gf-"+i);
				int totGenomes = paMatrix[i].length;
				
				for(int j=0;j<totGenomes;j++)
				{
					outPA.print("\t"+paMatrix[i][j]);

				}
			
				HomologyCluster cluster = homologyClusters.get(i);			
				outFreq.print(i);
				outFreq.print("\t" + cluster.getFrequency());
				outFreq.print("\t" + cluster.getExactCategory());
				outFreq.print("\t" + cluster.getSoftCategory());
				outFreq.println();
				
			}
		}
		
		catch (Exception e) {
			log.warning ("" + e);
		}
	}
	
	private void printGenomeMetadata(int id, QualifiedSequenceList sequencesMetadata, String jsFilename) throws IOException {
		String outFilename = outputPrefix+"_genome"+id+".tsv";
		try (PrintStream out = new PrintStream(outFilename);
			 PrintStream outJS = new PrintStream(new FileOutputStream(jsFilename, true))) {
			outJS.println("const genome" + id + " = [");
			out.println("Name\tLength");
			for(QualifiedSequence seq:sequencesMetadata) {
				out.println(""+seq.getName()+"\t"+seq.getLength());
				outJS.println("{Name: '"+seq.getName()+"', Length: "+seq.getLength()+"},");
			}
			outJS.println("];");
		}
	}


	private void printOrthologyUnit(HomologyUnit unit, PrintStream out, PrintStream outJS, Map<String,Integer> clusterIdsByHomologyUnit, Map<String,Integer> syntenyBlockIdByPair) {
		int clusterId1 = clusterIdsByHomologyUnit.getOrDefault(unit.getUniqueKey(), -1);
		for(int i=0;i<genomes.size();i++) {
			AnnotatedReferenceGenome genome = genomes.get(i);
			int genomeId = genome.getId();
			if(genomeId==unit.getGenomeId()) continue;
			Collection<HomologyEdge> orthologRelationships = unit.getOrthologRelationships(genomeId);
			if(orthologRelationships.size()==0) return;
			for(HomologyEdge edge:orthologRelationships) {
				HomologyUnit ortholog = edge.getSubjectUnit();
				int clusterId2 = clusterIdsByHomologyUnit.getOrDefault(ortholog.getUniqueKey(), -1);
				int syntenyBlock = syntenyBlockIdByPair.getOrDefault(unit.getUniqueKey()+"-"+ortholog.getUniqueKey(), -1);
				out.print(unit.getId()+"\t"+unit.getSequenceName()+"\t"+unit.getFirst()+"\t"+unit.getLast()+"\t"+clusterId1);
				out.print("\t"+ortholog.getGenomeId()+"\t"+ortholog.getId()+"\t"+ortholog.getSequenceName()+"\t"+ortholog.getFirst()+"\t"+ortholog.getLast()+"\t"+clusterId2);
				out.println("\t"+edge.getScore()+"\t"+syntenyBlock);
				outJS.println("{geneId: '"+unit.getId()
				+"', chromosome: '"+unit.getSequenceName()
				+"', geneStart: "+unit.getFirst()
				+", geneEnd: "+unit.getLast()
				+", genomeId2: '"+ortholog.getGenomeId()
				+"', geneIdG2: '"+ortholog.getId()
				+"', chromosomeG2: '"+ortholog.getSequenceName()
				+"', geneStartG2: "+ortholog.getFirst()
				+", geneEndG2: "+ortholog.getLast()
				+", score: "+edge.getScore()
				+", block: '"+syntenyBlock+"'},");
			}
		}
		
	}

	private void printD3Visualization(PrintStream outD3, String jsFile, String dataJsFilename, int preferredD3Version) throws IOException {
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
		outD3.println("<script src="+dataJsFilename+"></script>");
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
