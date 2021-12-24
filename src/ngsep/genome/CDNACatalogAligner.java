package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.ProteinTranslator;

public class CDNACatalogAligner {
	// Constants for default values
	public static final int INPUT_TYPE_CDNA = HomologyCatalog.INPUT_TYPE_CDNA;
	public static final int INPUT_TYPE_PROTEIN = HomologyCatalog.INPUT_TYPE_PROTEIN;
	public static final String DEF_OUT_PREFIX = "catalogsAlignment";
	public static final byte DEF_KMER_LENGTH = HomologRelationshipsFinder.DEF_KMER_LENGTH;
	public static final int DEF_MIN_PCT_KMERS = HomologRelationshipsFinder.DEF_MIN_PCT_KMERS;
	public static final int DEF_MAX_HOMOLOGS_UNIT = 3;
	
	// Logging and progress
	private Logger log = Logger.getLogger(CDNACatalogAligner.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String outputPrefix = DEF_OUT_PREFIX;
	private boolean skipMCL= false;
	private int inputType = INPUT_TYPE_CDNA;
	
	// Model attributes
	private HomologRelationshipsFinder homologRelationshipsFinder = new HomologRelationshipsFinder();
	private List<HomologyCatalog> homologyCatalogs = new ArrayList<>();
	private List<HomologyCluster> orthologyUnitClusters=new ArrayList<>();
	
	//logging
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
	
	//program arguments
	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	//program arguments
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
	public int getInputType() {
		return inputType;
	}
	public void setInputType(int inputType) {
		this.inputType = inputType;
	}
	public void setInputType(String value) {
		setInputType((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public static void main(String[] args) throws Exception {
		CDNACatalogAligner instance = new CDNACatalogAligner();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		while(i<args.length) {
			String fileOrganism = args[i++];
			instance.loadFile(fileOrganism);
		}
		instance.run();
	}
	
	private void loadFile(String fileName) throws IOException {
		HomologyCatalog catalog = HomologyCatalog.loadFromFasta(fileName,homologyCatalogs.size()+1,inputType);
		homologyCatalogs.add(catalog);
	}
	
	public void run () throws IOException {
		logParameters();
		if(homologyCatalogs.size()==0) throw new IOException("At least one organism's data should be provided");
		if(outputPrefix==null) throw new IOException("A prefix for output files is required");
		generateOrthologs();
		printRelationships(outputPrefix+"_rawHomologs.txt");
		generateClusters();
		printRelationships(outputPrefix+"_finalHomologs.txt");
		printResults(outputPrefix,orthologyUnitClusters);
		log.info("Process finished");
	}
	
	private void generateOrthologs() {
		catalogsDescription();
		int n = homologyCatalogs.size();
		for(int i=0;i<n;i++) {
			HomologyCatalog catalog = homologyCatalogs.get(i);
			List<HomologyEdge> edges = homologRelationshipsFinder.calculateParalogsOrganism(catalog);
			log.info(String.format("Paralogs found for Organism #%d: %d", i+1, edges.size()));
		}
		
		
		for(int i=0;i<n;i++) {
			HomologyCatalog catalog1 = homologyCatalogs.get(i);
			for (int j=0;j<n;j++) {
				HomologyCatalog catalog2 = homologyCatalogs.get(j);
				if(i!=j) {
					List<HomologyEdge> edges = homologRelationshipsFinder.calculateOrthologs(catalog1, catalog2);
					log.info(String.format("Orthologs found for Organisms #%d #%d: %d", i+1, j+1, edges.size()));
				}
			}
		}
	}
	
	private void catalogsDescription() {
		log.info("Total number of catalogs: " + homologyCatalogs.size());
		for(int i = 0; i < homologyCatalogs.size(); i++) {
			HomologyCatalog catalog = homologyCatalogs.get(i);
			log.info(String.format("Catalog #%d has %d genes.", i+1, catalog.getHomologyUnits().size()));
		}
	}
	
	public void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Loaded: "+ homologyCatalogs.size()+" annotated genomes");
		out.println("Output prefix:"+ outputPrefix);
		out.println("K-mer length: "+ getKmerLength());
		out.println("Minimum percentage of k-mers to call orthologs: "+ getMinPctKmers());
		log.info(os.toString());
	}
	
	private void generateClusters() {
		HomologClustersCalculator calculator = new HomologClustersCalculator(skipMCL);
		calculator.setLog(log);
		calculator.setKmerLength(homologRelationshipsFinder.getKmerLength());
		calculator.setMinPctKmers(homologRelationshipsFinder.getMinPctKmers());
		orthologyUnitClusters = calculator.clusterHomologsCatalogs(homologyCatalogs);
	}
	
	private void printRelationships(String filename) throws IOException {
		try (PrintStream out = new PrintStream(filename);) {
			out.println("Unit1\tUnit2\tscore");
			for(HomologyCatalog catalog:homologyCatalogs) {
				List<HomologyUnit> units = catalog.getHomologyUnits();
				for(HomologyUnit unit:units) {
					Collection<HomologyEdge> edges = unit.getAllHomologyRelationships();
					for(HomologyEdge edge:edges) {
						HomologyUnit target = edge.getSubjectUnit();
						out.println(unit.getId()+"\t"+target.getId()+"\t"+edge.getScore());
					}
				}
			}
		}
	}
	public static void printResults(String outputPrefix, List<HomologyCluster> orthologyUnitClusters) throws FileNotFoundException {
		//Print ortholog clusters
		try (PrintStream outClusters = new PrintStream(outputPrefix+"_clusters.txt");) {
			for(HomologyCluster cluster:orthologyUnitClusters) {
				List<HomologyUnit> memCluster = cluster.getHomologyUnitsCluster();
				outClusters.print(cluster.getClusterId());
				for(int i=0;i<memCluster.size();i++) {
					HomologyUnit unit = memCluster.get(i);
					outClusters.print("\t"+unit.getId());
				}
				outClusters.println();
			}
		}
	}
}
