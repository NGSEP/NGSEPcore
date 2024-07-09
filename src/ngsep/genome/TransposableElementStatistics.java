package ngsep.genome;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;

public class TransposableElementStatistics {
	// Logging and progress
	private Logger log = Logger.getLogger(TransposableElementsFinder.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String outputFile = null;
	
	// Model attributes
	private Map<String,TEFamilyStatistics> familyStatistics = new TreeMap<>();
	
	
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
	
	public static void main(String[] args) throws Exception {
		TransposableElementStatistics instance = new TransposableElementStatistics();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	public void run() throws IOException {
		logParameters ();
		if(inputFile==null) throw new IOException("The input genome is a required parameter");
		if(outputFile==null) throw new IOException("The output file is a required parameter");
		calculateStatistics(inputFile);
		try (PrintStream out = new PrintStream(outputFile)) {
			printStatistics(out);
		}
		
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		log.info(os.toString());
	}
	
	public void calculateStatistics(String filename)  throws IOException {
		familyStatistics = new TreeMap<>();
		try (FileReader reader = new FileReader(filename);
			 BufferedReader in = new BufferedReader(reader)) {
			String line = in.readLine();
			for (int i = 0; line != null; i++) {
				String[] items = line.split("\t|#");
				if(items.length<6) {
					line = in.readLine();
					continue;
				}
				GenomicRegion event = new TransposableElementAnnotation(items[0], Integer.parseInt(items[1]), Integer.parseInt(items[2]));
				TEFamilyStatistics stats = familyStatistics.computeIfAbsent(items[5], v->new TEFamilyStatistics(items[5]));
				stats.addEvent(event);
				line = in.readLine();
			}
		}
	}
	public void printStatistics(PrintStream out) {
		for(Map.Entry<String,TEFamilyStatistics> entry:familyStatistics.entrySet() ) {
			TEFamilyStatistics stats = entry.getValue();
			out.println(entry.getKey()+"\t"+stats.getNumEvents()+"\t"+stats.getTotalSpan());
		}	
	}

}
class TEFamilyStatistics {
	private String familyName;
	private GenomicRegionSortedCollection<GenomicRegion> events = new GenomicRegionSortedCollection<>();
	
	public TEFamilyStatistics(String familyName) {
		super();
		this.familyName = familyName;
	}
	public void addEvent(GenomicRegion event) {
		events.add(event);
	}
	
	public String getFamilyName() {
		return familyName;
	}
	public int getNumEvents() {
		return events.size();
	}
	public int getTotalSpan() {
		return events.getTotalSpan();
	}
}
