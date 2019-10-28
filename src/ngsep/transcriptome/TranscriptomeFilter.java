package ngsep.transcriptome;

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

	public static final int DEF_MIN_PROTEIN_LENGTH=0;
	public static final byte FORMAT_GFF = 0;
	public static final byte FORMAT_GENE_LIST = 1;
	public static final byte FORMAT_GENE_REGIONS = 2;
	public static final byte FORMAT_TRANSCRIPT_LIST = 3;
	public static final byte FORMAT_TRANSCRIPT_REGIONS = 4;
	private Logger log = Logger.getLogger(TranscriptomeFilter.class.getName());
	
	private ProgressNotifier progressNotifier=null;
	
	private ReferenceGenome genome;
	
	private GenomicRegionSortedCollection<GenomicRegion> regionsToFilter = null;
    private GenomicRegionSortedCollection<GenomicRegion> regionsToSelect = null;
    private boolean selectCompleteProteins = false;
    private int minProteinLength=DEF_MIN_PROTEIN_LENGTH;
    
	private byte outputFormat = FORMAT_GFF;
	
	private ProteinTranslator translator=  new ProteinTranslator();
	
	
	public static void main(String[] args) throws Exception {
		TranscriptomeFilter instance = new TranscriptomeFilter();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.genome = new ReferenceGenome(args[i++]);
		String transcriptomeFile = args[i++];
		instance.filterTranscriptome(transcriptomeFile, System.out);
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}

	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
    
    public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	
	public ReferenceGenome getGenome() {
		return genome;
	}

	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
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
	
	/**
	 * @return the selectCompleteProteins
	 */
	public boolean isSelectCompleteProteins() {
		return selectCompleteProteins;
	}

	/**
	 * @param selectCompleteProteins the selectCompleteProteins to set
	 */
	public void setSelectCompleteProteins(boolean selectCompleteProteins) {
		this.selectCompleteProteins = selectCompleteProteins;
	}
	public void setSelectCompleteProteins(Boolean selectCompleteProteins) {
		this.setSelectCompleteProteins(selectCompleteProteins.booleanValue());
	}
	
	/**
	 * @return the minProteinLength
	 */
	public int getMinProteinLength() {
		return minProteinLength;
	}

	/**
	 * @param minProteinLength the minProteinLength to set
	 */
	public void setMinProteinLength(int minProteinLength) {
		this.minProteinLength = minProteinLength;
	}
	
	public void setMinProteinLength(String value) {
		this.setMinProteinLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	

	/**
	 * @return the outputFormat
	 */
	public byte getOutputFormat() {
		return outputFormat;
	}

	/**
	 * @param outputFormat the outputFormat to set
	 */
	public void setOutputFormat(byte outputFormat) {
		if(outputFormat<FORMAT_GFF || outputFormat>FORMAT_TRANSCRIPT_REGIONS) throw new IllegalArgumentException("Invalid output format: "+outputFormat);
		this.outputFormat = outputFormat;
	}
	
	public void setOutputFormat(String value) {
		this.setOutputFormat((byte)OptionValuesDecoder.decode(value, Byte.class));
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
