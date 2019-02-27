/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.variants;

public class GenomicVariantAnnotation {
	//Reserved 1000g info attributes
	public static final String ATTRIBUTE_ANCESTRAL_ALLELE = "AA";
	public static final String ATTRIBUTE_PRIOR_ALLELE_FREQ = "AF";
	public static final String ATTRIBUTE_NUMBER_ALLELES = "AN";
	public static final String ATTRIBUTE_CIGAR = "CIGAR";
	public static final String ATTRIBUTE_DBSNP = "DB";
	public static final String ATTRIBUTE_END = "END";
	public static final String ATTRIBUTE_EXPECTED_ALLELE_COUNTS = "EC";
	public static final String ATTRIBUTE_HAPMAP2 = "H2";
	public static final String ATTRIBUTE_HAPMAP3 = "H3";
	public static final String ATTRIBUTE_SAMPLES_MQ_ZERO = "MQ0";
	public static final String ATTRIBUTE_SAMPLES_GENOTYPED = "NS";
	public static final String ATTRIBUTE_SOMATIC = "SOMATIC";
	public static final String ATTRIBUTE_VALIDATED = "VALIDATED";
	public static final String ATTRIBUTE_1000G = "1000G";
	//New annotation attributes
	/**
	 * Code of the annotation for closeness to a transcript
	 */
	public static final String ATTRIBUTE_TRANSCRIPT_ANNOTATION = "TA";
	/**
	 * Code of the id of the source transcript for the TA annotation 
	 */
	public static final String ATTRIBUTE_TRANSCRIPT_ID = "TID";
	/**
	 * Name of the gene related with the variant annotation 
	 */
	public static final String ATTRIBUTE_GENE_NAME = "TGN";
	/**
	 * For coding variants, first codon affected by the variant, and start codon position
	 */
	public static final String ATTRIBUTE_TRANSCRIPT_CODON = "TCO";
	/**
	 * For non-synonymous coding variants, aminoacid change
	 */
	public static final String ATTRIBUTE_TRANSCRIPT_AMINOACID_CHANGE = "TACH";
	/**
	 * Tells the number of samples with CNVs spanning this region
	 */
	public static final String ATTRIBUTE_IN_CNV = "CNV";
	/**
	 * Tells the minor allele frequency
	 */
	public static final String ATTRIBUTE_MAF = "MAF";
	/**
	 * Type of variant. See ngsep.variants.GenomicVariant interface for possible types
	 */
	public static final String ATTRIBUTE_TYPE = "TYPE";
	/**
	 * Counts of alleles observed in the population, including the reference allele
	 */
	public static final String ATTRIBUTE_ALLELE_FREQUENCY_SPECTRUM = "AFS";
	
	/**
	 * Phred scaled fisher strand bias 
	 */
	public static final String ATTRIBUTE_FISHER_STRAND_BIAS = "FS";
	
	private GenomicVariant genomicVariant;
	private String attribute;
	private Object value;
	public GenomicVariantAnnotation(GenomicVariant genomicVariant, String attribute, Object value) {
		super();
		this.genomicVariant = genomicVariant;
		//TODO: check if it is a known attribute and replace the key to save memory
		this.attribute = attribute;
		this.value = value;
	}
	public GenomicVariant getGenomicVariant() {
		return genomicVariant;
	}
	public String getAttribute() {
		return attribute;
	}
	public Object getValue() {
		return value;
	}
	public void setValue(Object value) {
		this.value = value;
	}
}
