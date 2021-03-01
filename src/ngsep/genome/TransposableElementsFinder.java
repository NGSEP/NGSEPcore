package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;

public class TransposableElementsFinder {

	// Logging and progress
	private Logger log = Logger.getLogger(TransposableElementsFinder.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String outputFile = null;
	private String transposonsDatabaseFile = null;
	
	// Model attributes
	
	// Get and set methods
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
	
	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	public String getTransposonsDatabaseFile() {
		return transposonsDatabaseFile;
	}
	public void setTransposonsDatabaseFile(String transposonsDatabaseFile) {
		this.transposonsDatabaseFile = transposonsDatabaseFile;
	}
	
	public static void main(String[] args) throws Exception {
		TransposableElementsFinder instance = new TransposableElementsFinder();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	public void run() throws IOException {
		long time = System.currentTimeMillis();
		logParameters ();
		if(inputFile==null) throw new IOException("The input genome is a required parameter");
		if(outputFile==null) throw new IOException("The output file is a required parameter");
		ReferenceGenome genome = new ReferenceGenome(inputFile);
		List<TransposableElementAnnotation> transposonAnnotations = findTransposons(genome);
		saveTransposons(transposonAnnotations,outputFile);
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		log.info("Process finished in "+seconds+" seconds");
	}
	public void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		if(transposonsDatabaseFile!=null) out.println("Database of transposable elements: "+ transposonsDatabaseFile);
		log.info(os.toString());
	}
	public List<TransposableElementAnnotation> findTransposons(ReferenceGenome genome) {
		GenomicRegionSortedCollection<TransposableElementAnnotation> annotations = new GenomicRegionSortedCollection<TransposableElementAnnotation>(genome.getSequencesMetadata());
		annotations.addAll(findTransposonsDeNovo(genome));
		if(transposonsDatabaseFile!=null) annotations.addAll(findTransposonsBySimilarity(genome));
		return removeRedundantAnnotations(annotations);
	}
	
	private List<TransposableElementAnnotation> findTransposonsDeNovo(ReferenceGenome genome) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		// TODO implement
		return answer;
	}
	private List<TransposableElementAnnotation> findTransposonsBySimilarity(ReferenceGenome genome) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		// TODO implement
		return answer;
	}
	private List<TransposableElementAnnotation> removeRedundantAnnotations(GenomicRegionSortedCollection<TransposableElementAnnotation> annotations) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		// TODO implement
		return answer;
	}
	public void saveTransposons(List<TransposableElementAnnotation> transposonAnnotations, String outputFile) {
		// TODO Auto-generated method stub
		
	}
	

}
