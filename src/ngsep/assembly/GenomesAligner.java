package ngsep.assembly;

import java.io.IOException;
import java.io.PrintStream;

import ngsep.genome.GenomicRegion;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;

public class GenomesAligner {

	private ReferenceGenome genome1;
	private ReferenceGenome genome2;
	private Transcriptome transcriptome1;
	private Transcriptome transcriptome2;
	
	
	
	public static void main(String[] args) throws Exception {
		GenomesAligner instance = new GenomesAligner();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String fileGenome1 = args[i++];
		String fileTranscriptome1 = args[i++];
		String fileGenome2 = args[i++];
		String fileTranscriptome2 = args[i++];
		String outPrefix = args[i++]; 
		instance.loadGenomes(fileGenome1,fileTranscriptome1, fileGenome2, fileTranscriptome2);
		instance.findUniqueCopyOrthologs();
		instance.alignOrthologUnits();
		try (PrintStream outAlignmnent = new PrintStream(outPrefix+"_lcs.txt");) {
			instance.printAlignmentOrthologUnits (outAlignmnent);
		}
	}

	public void loadGenomes(String fileGenome1, String fileTranscriptome1, String fileGenome2, String fileTranscriptome2) throws IOException {
		genome1 = new ReferenceGenome(fileGenome1);
		genome2 = new ReferenceGenome(fileGenome2);
		GFF3TranscriptomeHandler transcriptomeHandler = new GFF3TranscriptomeHandler();
		transcriptome1 = transcriptomeHandler.loadMap(fileTranscriptome1);
		transcriptome1.fillSequenceTranscripts(genome1);
		transcriptome2 = transcriptomeHandler.loadMap(fileTranscriptome2);
		transcriptome2.fillSequenceTranscripts(genome2);
		
	}



	public void findUniqueCopyOrthologs() {
		//Create orthology units based on transcripts
		
		//Build two FMIndex objects for proteins in the two genomes
		
		
		//Query each orthology unit against its own genome and select the units that are unique
		
		//Query selected units against the other genome to assign mates
		
	}



	public void alignOrthologUnits() {
		// TODO: Implement LCS algorithm to align ortholog units
		
	}



	private void printAlignmentOrthologUnits(PrintStream outAlignmnent) {
		// TODO Auto-generated method stub
		
	}

}
class OrthologyUnit implements GenomicRegion {
	private String id;
	private String sequenceName;
	private int first;
	private int last;
	private boolean negativeStrand = false;
	private String proteinSequence;
	private String mateId;
	
	public OrthologyUnit(String id, String sequenceName, int first, int last) {
		super();
		this.id = id;
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
	}
	@Override
	public String getSequenceName() {
		return sequenceName;
	}
	@Override
	public int getFirst() {
		return first;
	}
	@Override
	public int getLast() {
		return last;
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return !negativeStrand;
	}
	@Override
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	/**
	 * @param negativeStrand the negativeStrand to set
	 */
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
	/**
	 * @return the proteinSequence
	 */
	public String getProteinSequence() {
		return proteinSequence;
	}
	/**
	 * @param proteinSequence the proteinSequence to set
	 */
	public void setProteinSequence(String proteinSequence) {
		this.proteinSequence = proteinSequence;
	}
	/**
	 * @return the mateId
	 */
	public String getMateId() {
		return mateId;
	}
	/**
	 * @param mateId the mateId to set
	 */
	public void setMateId(String mateId) {
		this.mateId = mateId;
	}
}
