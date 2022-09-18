package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

public class GenomeAssemblyMask {
	
	// Logging and progress
	private Logger log = Logger.getLogger(TransposableElementsFinder.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	
	// Parameters 
	private String inputFile = null;
	private String outputFile = null;
	private String regionsFile = null;
	private boolean hardMask = false;
	
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

	public String getRegionsFile() {
		return regionsFile;
	}
	public void setRegionsFile(String regionsFile) {
		this.regionsFile = regionsFile;
	}
	
	public boolean isHardMask() {
		return hardMask;
	}
	public void setHardMask(boolean hardMask) {
		this.hardMask = hardMask;
	}
	public void setHardMask(Boolean hardMask) {
		this.setHardMask(hardMask.booleanValue());
	}
	
	public static void main(String[] args) throws Exception {
		GenomeAssemblyMask instance = new GenomeAssemblyMask();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	public void run() throws IOException {
		long time = System.currentTimeMillis();
		logParameters ();
		if(inputFile==null) throw new IOException("The input genome is a required parameter");
		if(outputFile==null) throw new IOException("The output file is a required parameter");
		if(regionsFile==null) throw new IOException("The regions file is a required parameter");
		ReferenceGenome genome = new ReferenceGenome(inputFile);
		Map<String,List<GenomicRegion>> regions = loadRegions(regionsFile);
		ReferenceGenome maskedGenome = maskGenome(genome,regions);
		saveGenome(maskedGenome,outputFile);
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		log.info("Process finished in "+seconds+" seconds");
		
	}
	private Map<String,List<GenomicRegion>> loadRegions(String regionsFile) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		Map<String,List<GenomicRegion>> regions = handler.loadRegionsAsMap(regionsFile);
		return regions;
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		if(regionsFile!=null) out.println("Regions file: "+ regionsFile);
		log.info(os.toString());
		
	}
	private ReferenceGenome maskGenome(ReferenceGenome genome, Map<String,List<GenomicRegion>> regions) throws IOException {
		QualifiedSequenceList sequences = genome.getSequencesList();
		QualifiedSequenceList newSequences = new QualifiedSequenceList();
		for(QualifiedSequence qseq:sequences) {
			String name = qseq.getName();
			if(qseq.getCharacters()==null) System.err.println("Null sequence for name: "+name);
			String seq = qseq.getCharacters().toString();
			StringBuilder newSeq = new StringBuilder();
			int nextPos = 0;
			List<GenomicRegion> regionsChr = regions.get(name);
			if(regionsChr==null) {
				log.info("No regions found for sequence: "+name);
				newSequences.add(qseq);
				continue;
			}
			log.info("Masking "+regionsChr.size()+" regions for sequence "+name);
			for(GenomicRegion r:regionsChr) {
				int firstZeroBased = r.getFirst()-1;
				int endZeroBased = r.getLast();
				if(nextPos<firstZeroBased) {
					String seqBefore = seq.substring(nextPos,firstZeroBased);
					newSeq.append(seqBefore.toUpperCase());
				}
				if(hardMask) {
					String nString = "N".repeat(endZeroBased-firstZeroBased);
					newSeq.append(nString);
				} else {
					String regionSeq = seq.substring(firstZeroBased,endZeroBased);
					newSeq.append(regionSeq.toLowerCase());
				}
				nextPos = endZeroBased;
			}
			String seqAfter = seq.substring(nextPos);
			newSeq.append(seqAfter.toUpperCase());
			newSequences.add(new QualifiedSequence(name,new DNAMaskedSequence(newSeq)));
		}
		return new ReferenceGenome(newSequences);
	}
	private void saveGenome(ReferenceGenome maskedGenome, String outputFile) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		try (PrintStream out = new PrintStream(outputFile)) {
			handler.saveSequences(maskedGenome.getSequencesList(), out, 100);
		}
		
		
	}

}
