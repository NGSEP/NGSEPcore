package ngsep.transcriptome;

import ngsep.variants.GenomicVariant;


public class VariantFunctionalAnnotation {
	
	
	public static final String ANNOTATION_SPLICE_DONOR="splice_donor_variant";
	public static final String ANNOTATION_SPLICE_ACCEPTOR="splice_acceptor_variant";
	public static final String ANNOTATION_SPLICE_REGION="splice_region_variant";
	public static final String ANNOTATION_FRAMESHIFT="frameshift_variant";
	public static final String ANNOTATION_NONSENSE="stop_gained";
	public static final String ANNOTATION_START_LOSS="start_lost";
	public static final String ANNOTATION_INFRAME_DEL="inframe_deletion";
	public static final String ANNOTATION_INFRAME_INS="inframe_insertion";
	public static final String ANNOTATION_STOP_LOSS="stop_lost";
	public static final String ANNOTATION_MISSENSE="missense_variant";
	public static final String ANNOTATION_SYNONYMOUS="synonymous_variant";
	
	public static final String ANNOTATION_CODING="coding_sequence_variant";
	public static final String ANNOTATION_5P_UTR="5_prime_UTR_variant";
	public static final String ANNOTATION_3P_UTR="3_prime_UTR_variant";
	public static final String ANNOTATION_NONCODINGRNA="non_coding_transcript_exon_variant";
	
	public static final String ANNOTATION_UPSTREAM="upstream_gene_variant";
	public static final String ANNOTATION_DOWNSTREAM="downstream_gene_variant";
	public static final String ANNOTATION_INTRON="intron_variant";
	public static final String ANNOTATION_INTERGENIC="intergenic_variant";
	
	
	private GenomicVariant variant;
	private String annotation;
	private Transcript transcript;
	private byte altAlleleIdx=-1;
	private int codonNumber=0;
	private byte codonPosition=0;
	private String aminoacidChange;
	public VariantFunctionalAnnotation(GenomicVariant variant, String annotation) {
		super();
		this.variant = variant;
		this.annotation = annotation;
	}
	/**
	 * @return the variant
	 */
	public GenomicVariant getVariant() {
		return variant;
	}
	/**
	 * @return the annotation
	 */
	public String getAnnotation() {
		return annotation;
	}
	/**
	 * @return the transcript
	 */
	public Transcript getTranscript() {
		return transcript;
	}
	/**
	 * @param transcript the transcript to set
	 */
	public void setTranscript(Transcript transcript) {
		this.transcript = transcript;
	}

	/**
	 * @param altAlleleIdx the altAlleleIdx to set
	 */
	public void setAltAlleleIdx(byte altAlleleIdx) {
		this.altAlleleIdx = altAlleleIdx;
	}
	/**
	 * @return the altAlleleIdx
	 */
	public int getAltAlleleIdx() {
		return altAlleleIdx;
	}
	/**
	 * @return the codonNumber
	 */
	public int getCodonNumber() {
		return codonNumber;
	}
	/**
	 * @param codonNumber the codonNumber to set
	 */
	public void setCodonNumber(int codonNumber) {
		this.codonNumber = codonNumber;
	}
	/**
	 * @return the codonPosition
	 */
	public byte getCodonPosition() {
		return codonPosition;
	}
	/**
	 * @param codonPosition the codonPosition to set
	 */
	public void setCodonPosition(byte codonPosition) {
		this.codonPosition = codonPosition;
	}
	/**
	 * @return the aminoacidChange
	 */
	public String getAminoacidChange() {
		return aminoacidChange;
	}
	/**
	 * @param aminoacidChange the aminoacidChange to set
	 */
	public void setAminoacidChange(String aminoacidChange) {
		this.aminoacidChange = aminoacidChange;
	}
}
