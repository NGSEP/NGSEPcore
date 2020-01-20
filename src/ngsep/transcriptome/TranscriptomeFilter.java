package ngsep.transcriptome;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;
import ngsep.transcriptome.io.GFF3TranscriptomeWriter;

public class TranscriptomeFilter {

	// Constants for default values
	public static final int DEF_MIN_PROTEIN_LENGTH=0;
	public static final byte FORMAT_GFF = 0;
	public static final byte FORMAT_GENE_LIST = 1;
	public static final byte FORMAT_GENE_REGIONS = 2;
	public static final byte FORMAT_TRANSCRIPT_LIST = 3;
	public static final byte FORMAT_TRANSCRIPT_REGIONS = 4;
	
	// Logging and progress
	private Logger log = Logger.getLogger(TranscriptomeFilter.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String outputFile = null;
	private String genomeFile = null;
	private ReferenceGenome genome;
	private byte outputFormat = FORMAT_GFF;
	private boolean selectCompleteProteins = false;
    private int minProteinLength=DEF_MIN_PROTEIN_LENGTH;
	private GenomicRegionSortedCollection<GenomicRegion> regionsToFilter = null;
    private GenomicRegionSortedCollection<GenomicRegion> regionsToSelect = null;
	
	private ProteinTranslator translator=  new ProteinTranslator();
	
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
	
	public String getGenomeFile() {
		return genomeFile;
	}
	public void setGenomeFile(String genomeFile) {
		this.genomeFile = genomeFile;
	}
	public ReferenceGenome getGenome() {
		return genome;
	}

	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	
	public byte getOutputFormat() {
		return outputFormat;
	}
	public void setOutputFormat(byte outputFormat) {
		if(outputFormat<FORMAT_GFF || outputFormat>FORMAT_TRANSCRIPT_REGIONS) throw new IllegalArgumentException("Invalid output format: "+outputFormat);
		this.outputFormat = outputFormat;
	}
	public void setOutputFormat(String value) {
		this.setOutputFormat((byte)OptionValuesDecoder.decode(value, Byte.class));
	}

	public boolean isSelectCompleteProteins() {
		return selectCompleteProteins;
	}
	public void setSelectCompleteProteins(boolean selectCompleteProteins) {
		this.selectCompleteProteins = selectCompleteProteins;
	}
	public void setSelectCompleteProteins(Boolean selectCompleteProteins) {
		this.setSelectCompleteProteins(selectCompleteProteins.booleanValue());
	}
	
	public int getMinProteinLength() {
		return minProteinLength;
	}
	public void setMinProteinLength(int minProteinLength) {
		this.minProteinLength = minProteinLength;
	}
	public void setMinProteinLength(String value) {
		this.setMinProteinLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public List<GenomicRegion> getRegionsToFilter() {
		return regionsToFilter.asList();
	}
	public void setRegionsToFilter(List<GenomicRegion> regions) {
		this.regionsToFilter = new GenomicRegionSortedCollection<GenomicRegion>(regions);
	}
	public void setRegionsToFilter(String regionsFile) throws IOException {
		if(regionsFile==null || regionsFile.length()==0) {
			this.regionsToFilter = null;
			return;
		}
		SimpleGenomicRegionFileHandler regionFileHandler = new SimpleGenomicRegionFileHandler();
		List<GenomicRegion> regions = regionFileHandler.loadRegions(regionsFile);
		this.regionsToFilter = new GenomicRegionSortedCollection<GenomicRegion>(regions);
	}

	public List<GenomicRegion> getRegionsToSelect() {
		return regionsToSelect.asList();
	}
	public void setRegionsToSelect(List<GenomicRegion> regions) {
		this.regionsToSelect = new GenomicRegionSortedCollection<GenomicRegion>(regions);
	}
	public void setRegionsToSelect(String regionsFile) throws IOException {
		if(regionsFile==null || regionsFile.length()==0) {
			this.regionsToSelect = null;
			return;
		}
		SimpleGenomicRegionFileHandler regionFileHandler = new SimpleGenomicRegionFileHandler();
		List<GenomicRegion> regions = regionFileHandler.loadRegions(regionsFile);
		this.regionsToSelect = new GenomicRegionSortedCollection<GenomicRegion>(regions);
	}
	
	public static void main(String[] args) throws Exception {
		TranscriptomeFilter instance = new TranscriptomeFilter();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}

	public void run() throws IOException {
		logParameters();
		if (genomeFile!=null ) {
			log.info("Loading genome from: "+genomeFile);
			setGenome(new ReferenceGenome(genomeFile));
			log.info("Loaded genome with: "+genome.getNumSequences()+" sequences. Total length: "+genome.getTotalLength());
		} else if (genome != null) {
			log.info("Running with loaded genome from: "+genome.getFilename()+" number of sequences: "+genome.getNumSequences()+". Total length: "+genome.getTotalLength());
		} else {
			throw new IOException("The file with the reference genome is a required parameter");
		}
		if(inputFile==null) throw new IOException("The input transcriptome file is a required parameter");
		if(outputFile==null) filterTranscriptome(inputFile, System.out);
		else {
			try (PrintStream out = new PrintStream(outputFile)) {
				filterTranscriptome(inputFile, System.out);
			}
		}
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		if(outputFile!=null) out.println("Output file:"+ outputFile);
		else out.println("Write to standard output");
		if (genomeFile!=null) out.println("Given reference genome file: "+genomeFile);
		else if (genome!=null) out.println("Reference genome with "+genome.getNumSequences()+" and total length: "+genome.getTotalLength()+" previously loaded from: "+genome.getFilename());
		if (outputFormat == FORMAT_GFF) out.println("Output GFF");
		if (outputFormat == FORMAT_GENE_LIST) out.println("Output list of genes");
		if (outputFormat == FORMAT_GENE_REGIONS) out.println("Output list of gene annotated genomic regions");
		if (outputFormat == FORMAT_TRANSCRIPT_LIST) out.println("Output list of transcripts");
		if (outputFormat == FORMAT_TRANSCRIPT_REGIONS) out.println("Output list of transcript annotated genomic regions");
		if(selectCompleteProteins) out.println("Select only proteins with valid start and stop codons");
		out.println("Minimum protein length: "+ minProteinLength);
		if(regionsToFilter!=null) out.println("Loaded "+ regionsToFilter.size()+" regions to filter");
		if(regionsToSelect!=null) out.println("Loaded "+ regionsToSelect.size()+" regions to select");
		log.info(os.toString());
		
	}
	/**
	 * Filters the genes and transcripts in given transcriptome
	 * @param transcriptomeFile in GFF format with transcripts to filter
	 * @param out Stream for the output file
	 * @throws IOException 
	 */
	public void filterTranscriptome(String transcriptomeFile, PrintStream out) throws IOException {
		QualifiedSequenceList sequences = genome.getSequencesMetadata();
		GFF3TranscriptomeHandler handler = new GFF3TranscriptomeHandler(sequences);
		Transcriptome transcriptome = handler.loadMap(transcriptomeFile);
		transcriptome.fillSequenceTranscripts(genome, log);
		List<Transcript> transcripts = transcriptome.getAllTranscripts();
		Set<String> genes = new HashSet<>();
		int n=0;
		Transcriptome filteredTranscriptome = new Transcriptome(sequences);
		for(Transcript transcript:transcripts) {
			//if(n<100) System.out.println("Filtering transcript: "+transcript.getId());
			if(passFilters(transcript)) {
				//if(n<100) System.out.println("Transcript: "+transcript.getId()+" passed. Out format: "+outputFormat);
				if(outputFormat==FORMAT_TRANSCRIPT_LIST ) {
					out.println(transcript.getId());	
				} else if (outputFormat==FORMAT_TRANSCRIPT_REGIONS) {
					out.println(transcript.getSequenceName()+"\t"+transcript.getFirst()+"\t"+transcript.getLast()+"\t"+transcript.getId());
				} else if(outputFormat==FORMAT_GENE_LIST ) {
					String geneId = transcript.getGeneId();
					if(!genes.contains(geneId)) {
						out.println(geneId);
						genes.add(geneId);
					}
				} else if(outputFormat==FORMAT_GENE_LIST ) {
					String geneId = transcript.getGeneId();
					if(!genes.contains(geneId)) {
						Gene gene = transcript.getGene();
						out.println(gene.getSequenceName()+"\t"+gene.getFirst()+"\t"+gene.getLast()+"\t"+gene.getId());
						genes.add(geneId);
					}
				} else if (outputFormat==FORMAT_GFF ) {
					filteredTranscriptome.addTranscript(transcript);
				}
			}
			n++;
			if (progressNotifier!=null && n%100==0) {
				int progress = n/100;
				if (!progressNotifier.keepRunning(progress)) {
					return;
				}
			}
		}
		if (outputFormat==FORMAT_GFF ) {
			GFF3TranscriptomeWriter writer = new GFF3TranscriptomeWriter();
			writer.printTranscriptome(filteredTranscriptome, out);
		}
		
	}

	private boolean passFilters(Transcript transcript) {
		if (regionsToFilter!=null) {
			if(regionsToFilter.findSpanningRegions(transcript).size()>0) return false;
		}
		if (regionsToSelect!=null) {
			if(regionsToSelect.findSpanningRegions(transcript).size()==0) return false;
		}
		if(selectCompleteProteins || minProteinLength>0) {
			if(!transcript.isCoding()) return false;
			String protein = null;
			Codon startCodon = transcript.getStartCodon();
			Codon stopCodon = transcript.getStopCodon();	
			
			protein = transcript.getProteinSequence(translator);
			if(protein == null) {
				return false;
			}
			if(protein.length()<minProteinLength) return false;
			boolean complete = startCodon!=null && startCodon.isStart() && stopCodon!=null && stopCodon.isStop();
			if (selectCompleteProteins && !complete) return false;
		}
		return true;
	}

}
