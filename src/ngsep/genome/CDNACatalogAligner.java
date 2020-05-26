package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaFileReader;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.ProteinTranslator;
import ngsep.transcriptome.Transcript;

public class CDNACatalogAligner {
	// Constants for default values
	public static final String DEF_OUT_PREFIX = "catalogsAlignment";
	public static final byte DEF_KMER_LENGTH = HomologRelationshipsFinder.DEF_KMER_LENGTH;
	public static final int DEF_MIN_PCT_KMERS = HomologRelationshipsFinder.DEF_MIN_PCT_KMERS;
	public static final int DEF_MAX_HOMOLOGS_UNIT = 3;
	
	// Logging and progress
	private Logger log = Logger.getLogger(CDNACatalogAligner.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private ProteinTranslator translator = new ProteinTranslator();
	private String outputPrefix = DEF_OUT_PREFIX;
	private boolean skipMCL= false;
	
	// Model attributes
	private HomologRelationshipsFinder homologRelationshipsFinder = new HomologRelationshipsFinder();
	private List<HomologyCatalog> cdnaCatalogs = new ArrayList<>();
	private List<HomologyEdge> homologyEdges = new ArrayList<HomologyEdge>();
	private List<List<HomologyUnit>> orthologyUnitClusters=new ArrayList<>();
	
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
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> sequences = handler.loadSequences(fileName);
		
		List<HomologyUnit> units = new ArrayList<>();
		for(QualifiedSequence seq : sequences) {
			HomologyUnit unit = new HomologyUnit(cdnaCatalogs.size()+1, seq.getName(), translator.getProteinSequence(seq.getCharacters()));
			units.add(unit);
		}
		
		HomologyCatalog catalog = new HomologyCatalog(units);
		cdnaCatalogs.add(catalog);
	}
	
	public void run () throws IOException {
		logParameters();
		if(cdnaCatalogs.size()==0) throw new IOException("At least one organism's data should be provided");
		if(outputPrefix==null) throw new IOException("A prefix for output files is required");
		generateOrthologs();
		printPartialResults();
		generateClusters();
		printResults();
		log.info("Process finished");
	}
	
	private void generateOrthologs() {
		catalogsDescription();
		
		for(int i=0;i<cdnaCatalogs.size();i++) {
			HomologyCatalog catalog = cdnaCatalogs.get(i);
			List<HomologyEdge> edges = homologRelationshipsFinder.calculateParalogsOrganism(catalog);
			homologyEdges.addAll(edges);
			log.info(String.format("Paralogs found for Organism #%d: %d", i+1, edges.size()));
		}
		
		
		for(int i=0;i<cdnaCatalogs.size();i++) {
			HomologyCatalog catalog1 = cdnaCatalogs.get(i);
			for (int j=0;j<cdnaCatalogs.size();j++) {
				HomologyCatalog catalog2 = cdnaCatalogs.get(j);
				if(i!=j) {
					List<HomologyEdge> edges = homologRelationshipsFinder.calculateOrthologs(catalog1, catalog2);
					homologyEdges.addAll(edges);
					log.info(String.format("Orthologs found for Organisms #%d #%d: %d", i+1, j+1, edges.size()));
				}
			}
		}
	}
	
	private void catalogsDescription() {
		log.info("Total number of catalogs: " + cdnaCatalogs.size());
		for(int i = 0; i < cdnaCatalogs.size(); i++) {
			HomologyCatalog catalog = cdnaCatalogs.get(i);
			log.info(String.format("Catalog #%d has %d genes.", i+1, catalog.getHomologyUnits().size()));
		}
	}
	
	public void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Loaded: "+ cdnaCatalogs.size()+" annotated genomes");
		out.println("Output prefix:"+ outputPrefix);
		out.println("K-mer length: "+ getKmerLength());
		out.println("Minimum percentage of k-mers to call orthologs: "+ getMinPctKmers());
		log.info(os.toString());
	}
	
	private void generateClusters() {
		HomologClustersCalculator calculator = new HomologClustersCalculator(skipMCL);
		calculator.setLog(log);
		orthologyUnitClusters = calculator.clusterHomologsCatalogs(cdnaCatalogs, homologyEdges);
	}
	
	private void printPartialResults() throws FileNotFoundException {
		//Print orthology relationships
		try (PrintStream outOrthologs = new PrintStream(outputPrefix+"_rawOrthologs.txt");) {
			for(HomologyEdge edge : homologyEdges) {
				outOrthologs.print(String.format("%s\t%s\t%f", edge.getQueryUnit().getId(), edge.getSubjectUnit().getId(), edge.getScore()));
				outOrthologs.println();
			}
		}
	}
	
	public void printResults() throws FileNotFoundException {	
		//Print ortholog clusters
		try (PrintStream outClusters = new PrintStream(outputPrefix+"_clusters.txt");) {
			for(List<HomologyUnit> cluster:orthologyUnitClusters) {
				outClusters.print(cluster.get(0).getId());
				for(int i=1;i<cluster.size();i++) {
					HomologyUnit unit = cluster.get(i);
					outClusters.print("\t"+unit.getId());
				}
				outClusters.println();
			}
		}
	}
}
