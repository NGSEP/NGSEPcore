package ngsep.sequences;

import java.io.ByteArrayOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import ngsep.assembly.Assembler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.io.FastqFileReader;
import ngsep.sequences.io.SimpleSequenceListLoader;

public class FastqFileFilter {

	// Constants for default values
	public static final int DEF_MIN_READ_LENGTH = 0;
	public static final int DEF_MIN_READ_AVG_QUAL = 0;
	
	// Logging and progress
	private Logger log = Logger.getLogger(Assembler.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	// Parameters
	private String inputFile = null;
	private String outputFile = null;
	private int minReadLength = DEF_MIN_READ_LENGTH;
	private int minReadAverageQuality = DEF_MIN_READ_AVG_QUAL;
	private String selectReadIdsFile = null;
	
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
	public int getMinReadLength() {
		return minReadLength;
	}
	public void setMinReadLength(int minReadLength) {
		this.minReadLength = minReadLength;
	}
	public void setMinReadLength(String value) {
		setMinReadLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public int getMinReadAverageQuality() {
		return minReadAverageQuality;
	}
	public void setMinReadAverageQuality(int minReadAverageQuality) {
		this.minReadAverageQuality = minReadAverageQuality;
	}
	public void setMinReadAverageQuality(String value) {
		setMinReadAverageQuality((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public String getSelectReadIdsFile() {
		return selectReadIdsFile;
	}
	public void setSelectReadIdsFile(String selectReadIdsFile) {
		this.selectReadIdsFile = selectReadIdsFile;
	}
	
	public static void main(String[] args) throws Exception {
		FastqFileFilter instance = new FastqFileFilter();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	public void run() throws IOException {
		logParameters();
		if(inputFile==null) throw new IOException("The input file with raw reads is required");
		if(outputFile==null) throw new IOException("An output file is required");
		run (inputFile, outputFile);
		log.info("Process finished");
		
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		out.println("Minimum read length: "+minReadLength);
		out.println("Minimum average read quality: "+minReadAverageQuality);
		if(selectReadIdsFile!=null) out.println("File with read ids to select: "+selectReadIdsFile);
		log.info(os.toString());
	}
	public void run(String inputFile, String outputFile) throws IOException {
		int n = 0;
		int np = 0;
		Set<String> readIdsToSelect = null;
		if(selectReadIdsFile!=null) {
			readIdsToSelect = loadReadIdsFile(selectReadIdsFile);
		}
		if(!outputFile.toLowerCase().endsWith(".gz")) outputFile=outputFile+".gz";
		try (FastqFileReader reader = new FastqFileReader(inputFile);
			 OutputStream os = new GZIPOutputStream(new FileOutputStream(outputFile));
			 PrintStream out = new PrintStream(os)) {
			reader.setSequenceType(DNASequence.class);
			reader.setMinAverageQuality(minReadAverageQuality);
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				boolean pass = read.getLength()>=minReadLength;
				if(readIdsToSelect!=null) pass = pass && readIdsToSelect.contains(read.getName());
				if(pass) {
					read.save(out);
					np++;
				}
				n++;
			}
		}
		log.info("Total reads: "+n+". PASS reads: "+np);
	}
	private Set<String> loadReadIdsFile(String filename) throws IOException {
		SimpleSequenceListLoader loader = new SimpleSequenceListLoader();
		QualifiedSequenceList seqIds = loader.loadSequences(filename);
		Set<String> answer = new HashSet<String>();
		for(QualifiedSequence seq:seqIds) answer.add(seq.getName());
		return answer;
	}
}
