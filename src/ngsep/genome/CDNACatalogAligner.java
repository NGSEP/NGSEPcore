package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.main.ThreadPoolManager;

public class CDNACatalogAligner {
	// Constants for default values
	public static final int INPUT_TYPE_CDNA = HomologyCatalog.INPUT_TYPE_CDNA;
	public static final int INPUT_TYPE_PROTEIN = HomologyCatalog.INPUT_TYPE_PROTEIN;
	public static final String DEF_OUT_PREFIX = "catalogsAlignment";
	public static final byte DEF_KMER_LENGTH = HomologRelationshipsFinder.DEF_KMER_LENGTH;
	public static final int DEF_MIN_PCT_KMERS = HomologRelationshipsFinder.DEF_MIN_PCT_KMERS;
	public static final int DEF_MAX_HOMOLOGS_UNIT = 3;
	public static final int DEF_NUM_THREADS = 1;
	
	// Logging and progress
	private Logger log = Logger.getLogger(CDNACatalogAligner.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String outputPrefix = DEF_OUT_PREFIX;
	private boolean skipMCL= false;
	private int inputType = INPUT_TYPE_CDNA;
	private boolean calculateNucleotideEvolutionStatistics = false;
	private int numThreads = 1;
	
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
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	public void setNumThreads(String value) {
		setNumThreads((int)OptionValuesDecoder.decode(value, Integer.class));
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
		if(calculateNucleotideEvolutionStatistics) calculateNucleotideEvolutionStatistics();
		printResults(outputPrefix,orthologyUnitClusters);
		log.info("Process finished");
	}
	
	
	private void generateOrthologs() {
		catalogsDescription();
		int n = homologyCatalogs.size();
		ThreadPoolManager poolParalogs = new ThreadPoolManager(numThreads, numThreads);
		for(int i=0;i<n;i++) {
			final int index = i;
			try {
				poolParalogs.queueTask(()->calculateParalogs(index));
			} catch (InterruptedException e) {
				e.printStackTrace();
				throw new RuntimeException("Concurrence error calculating paralogs for genome "+i,e);
			} 
		}
		try {
			poolParalogs.terminatePool();
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException("Concurrence error calculating paralogs",e);
		}
		ThreadPoolManager poolOrthologs = new ThreadPoolManager(numThreads, numThreads);
		for(int i=0;i<n;i++) {
			for (int j=0;j<n;j++) {
				if(i==j) continue;
				try {
					final int index1=i;
					final int index2=j;
					poolOrthologs.queueTask(()->calculateOrthologs(index1, index2));
				} catch (InterruptedException e) {
					e.printStackTrace();
					throw new RuntimeException("Concurrence error calculating orthologs between "+i+" "+j,e);
				}
			}
		}
		try {
			poolOrthologs.terminatePool();
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException("Concurrence error calculating orthologs",e);
		}
	}
	private void calculateOrthologs(int i, int j) {
		HomologyCatalog catalog1 = homologyCatalogs.get(i);
		HomologyCatalog catalog2 = homologyCatalogs.get(j);
		List<HomologyEdge> edges = homologRelationshipsFinder.calculateOrthologs(catalog1.getHomologyUnits(), catalog2.getHomologyUnits());
		log.info(String.format("Orthologs found for Organisms #%d #%d: %d", i+1, j+1, edges.size()));
	}
	private void calculateParalogs(int i) {
		HomologyCatalog catalog = homologyCatalogs.get(i);
		List<HomologyEdge> edges = homologRelationshipsFinder.calculateParalogs(catalog.getHomologyUnits());
		log.info(String.format("Paralogs found for Organism #%d: %d", i+1, edges.size()));
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
		if(inputType==INPUT_TYPE_CDNA)  out.println("Loaded: "+ homologyCatalogs.size()+" transcriptomes");
		if(inputType==INPUT_TYPE_PROTEIN)  out.println("Loaded: "+ homologyCatalogs.size()+" proteomes");
		out.println("Output prefix:"+ outputPrefix);
		out.println("K-mer length: "+ getKmerLength());
		out.println("Minimum percentage of k-mers to call orthologs: "+ getMinPctKmers());
		if(skipMCL) out.println("Skip the MCL algorithm. Use connected components");
		out.println("Number of threads: "+ getNumThreads());
		log.info(os.toString());
	}
	
	private void generateClusters() {
		HomologClustersCalculator calculator = new HomologClustersCalculator(homologRelationshipsFinder,skipMCL);
		calculator.setLog(log);
		orthologyUnitClusters = calculator.clusterHomologsCatalogs(homologyCatalogs);
	}
	
	private void printRelationships(String filename) throws IOException {
		Map<String,Integer> clusterAssignment = new HashMap<String,Integer>();
		for(HomologyCluster cluster:orthologyUnitClusters) {
			for(HomologyUnit unit:cluster.getHomologyUnitsCluster()) {
				clusterAssignment.put(unit.getUniqueKey(), cluster.getClusterId());
			}
		}
		try (PrintStream out = new PrintStream(filename);) {
			out.println("Unit1\tUnit2\tscore\tClusterUnit1\tClusterUnit2");
			for(HomologyCatalog catalog:homologyCatalogs) {
				List<HomologyUnit> units = catalog.getHomologyUnits();
				for(HomologyUnit unit:units) {
					int clusterId1 = clusterAssignment.getOrDefault(unit.getUniqueKey(), -1);
					Collection<HomologyEdge> edges = unit.getAllHomologyRelationships();
					for(HomologyEdge edge:edges) {
						HomologyUnit target = edge.getSubjectUnit();
						int clusterId2 = clusterAssignment.getOrDefault(target.getUniqueKey(), -1);
						out.println(unit.getId()+"\t"+target.getId()+"\t"+edge.getScore()+"\t"+clusterId1+"\t"+clusterId2);
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
	private void calculateNucleotideEvolutionStatistics() {
		for(HomologyCluster cluster:orthologyUnitClusters) {
			calculateNucleotideEvolutionStatistics(cluster);
		}
		
	}
	private void calculateNucleotideEvolutionStatistics(HomologyCluster cluster) {
		for(HomologyUnit unit:cluster.getHomologyUnitsCluster()) {
			Collection<HomologyEdge> edges = unit.getAllHomologyRelationships();
			for(HomologyEdge edge:edges) {
				HomologyUnit u1 = edge.getQueryUnit();
				HomologyUnit u2 = edge.getSubjectUnit();
				CharSequence s1 = u1.getCdsSequence();
				CharSequence s2 = u2.getCdsSequence();
				//TODO: 1. Align sequences
				//2. Calculate statistics
				//3. Update edge with statistics
			}
		}
		
	}
}
