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
package ngsep.vcf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import ngsep.genome.GenomicRegion;
import ngsep.transcriptome.VariantFunctionalAnnotation;
import ngsep.transcriptome.VariantFunctionalAnnotationType;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.DiversityStatistics;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;

public class VCFRecord implements GenomicRegion {
	
	
	public static final String FORMAT_GENOTYPE = "GT";
	public static final String FORMAT_GENOTYPE_PHRED_LIKELIHOOD = "PL";
	public static final String FORMAT_GENOTYPE_LIKELIHOOD = "GL";
	public static final String FORMAT_GENOTYPE_QUALITY = "GQ";
	public static final String FORMAT_DEPTH = "DP";
	
	
	
	//New genotype attributes
	
	/**
	 * Number of base calls (depth) for the 4 nucleotides in called SNVs sorted as A,C,G,T
	 */
	public static final String FORMAT_ALL_BASES_SNP_DEPTH = "BSDP";
	
	/**
	 * Number of copies of each of the alleles at this site, taking into account copy number variation events surrounding the variant
	 */
	public static final String FORMAT_ALLELE_COPY_NUMBER = "ACN";
	
	/**
	 * Number of base calls (depth) for each allele, including the reference
	 */
	public static final String FORMAT_ALLELE_DEPTH = "ADP";
	/**
	 * Number of fragments supporting a structural event
	 */
	public static final String FORMAT_NUMBER_SUPPORTING_FRAGMENTS = "NSF";
	/**
	 * Number of split reads supporting a large indel
	 */
	public static final String FORMAT_NUMBER_SPLIT_READS = "NSR";
	/**
	 * Number of copies of a copy number variant as a real number
	 */
	public static final String FORMAT_REAL_NUMBER_COPIES = "RNC";
	/**
	 * Number of fragments showing evidence of tandem duplication
	 */
	public static final String FORMAT_NUMBER_TANDEM_DUP_FRAGS = "NTADF";
	/**
	 * Number of fragments showing evidence of non-tandem (trans) duplication
	 */
	public static final String FORMAT_NUMBER_TRANS_DUP_FRAGS = "NTRDF";
	/**
	 * Text describing the genotype for copy number variants
	 */
	public static final String FORMAT_TEXT_GENOTYPE = "TGEN";
	
	
	
	
	//Constants to keep arrays with format fields that can be understood by NGSEP
	public static final int FORMAT_IDX_GT = 0;
	public static final int FORMAT_IDX_PL = 1;
	public static final int FORMAT_IDX_GL = 2;
	public static final int FORMAT_IDX_GQ = 3;
	public static final int FORMAT_IDX_DP = 4;
	public static final int FORMAT_IDX_ADP = 5;
	public static final int FORMAT_IDX_BSDP = 6;
	public static final int FORMAT_IDX_ACN = 7;
	public static final int FORMAT_IDX_NSF = 8;
	public static final int FORMAT_IDX_NSR = 9;
	public static final int FORMAT_IDX_RNC = 10;
	public static final int FORMAT_IDX_NTADF = 11;
	public static final int FORMAT_IDX_NTRDF = 12;
	public static final int FORMAT_IDX_TGEN = 13;
	
	public static final int [] DEF_FORMAT_ARRAY_NONE = {};
	public static final int [] DEF_FORMAT_ARRAY_MINIMAL = {FORMAT_IDX_GT};
	public static final int [] DEF_FORMAT_ARRAY_QUALITY = {FORMAT_IDX_GT,FORMAT_IDX_GQ};
	public static final int [] DEF_FORMAT_ARRAY_NGSEP_SNV = {FORMAT_IDX_GT,FORMAT_IDX_PL,FORMAT_IDX_GQ,FORMAT_IDX_DP,FORMAT_IDX_BSDP,FORMAT_IDX_ACN};
	public static final int [] DEF_FORMAT_ARRAY_NGSEP_NOSNV = {FORMAT_IDX_GT,FORMAT_IDX_PL,FORMAT_IDX_GQ,FORMAT_IDX_DP,FORMAT_IDX_ADP,FORMAT_IDX_ACN};
	public static final Map<String, Integer> KNOWN_FORMAT_FIELDS_MAP = loadMap();
	private static Map<String, Integer> loadMap() {
		Map<String, Integer> answer = new TreeMap<String, Integer>();
		answer.put(FORMAT_GENOTYPE, FORMAT_IDX_GT);
		answer.put(FORMAT_GENOTYPE_PHRED_LIKELIHOOD, FORMAT_IDX_PL);
		answer.put(FORMAT_GENOTYPE_LIKELIHOOD, FORMAT_IDX_GL);
		answer.put(FORMAT_GENOTYPE_QUALITY, FORMAT_IDX_GQ);
		answer.put(FORMAT_DEPTH, FORMAT_IDX_DP);
		answer.put(FORMAT_ALLELE_DEPTH, FORMAT_IDX_ADP);
		answer.put(FORMAT_ALL_BASES_SNP_DEPTH, FORMAT_IDX_BSDP);
		answer.put(FORMAT_ALLELE_COPY_NUMBER, FORMAT_IDX_ACN);
		answer.put(FORMAT_NUMBER_SUPPORTING_FRAGMENTS, FORMAT_IDX_NSF);
		answer.put(FORMAT_NUMBER_SPLIT_READS, FORMAT_IDX_NSR);
		answer.put(FORMAT_REAL_NUMBER_COPIES, FORMAT_IDX_RNC);
		answer.put(FORMAT_NUMBER_TANDEM_DUP_FRAGS, FORMAT_IDX_NTADF);
		answer.put(FORMAT_NUMBER_TRANS_DUP_FRAGS, FORMAT_IDX_NTRDF);
		answer.put(FORMAT_TEXT_GENOTYPE, FORMAT_IDX_TGEN);
		
		//Kept for backwards compatibility
		answer.put("AC", FORMAT_IDX_ADP);
		answer.put("AAC", FORMAT_IDX_BSDP);
		
		//Added to load counts from GATK VCF files
		answer.put("AD", FORMAT_IDX_ADP);
		
		return answer;
	}
	
	public static final String [] KNOWN_FORMAT_FIELDS_ARRAY = loadFormatNames();
	
	private static String[] loadFormatNames() {
		String [] answer = new String[14];
		answer[FORMAT_IDX_GT] = FORMAT_GENOTYPE;
		answer[FORMAT_IDX_PL] = FORMAT_GENOTYPE_PHRED_LIKELIHOOD;
		answer[FORMAT_IDX_GL] = FORMAT_GENOTYPE_LIKELIHOOD;
		answer[FORMAT_IDX_GQ] = FORMAT_GENOTYPE_QUALITY;
		answer[FORMAT_IDX_DP] = FORMAT_DEPTH;
		answer[FORMAT_IDX_ADP] = FORMAT_ALLELE_DEPTH;
		answer[FORMAT_IDX_BSDP] = FORMAT_ALL_BASES_SNP_DEPTH;
		answer[FORMAT_IDX_ACN] = FORMAT_ALLELE_COPY_NUMBER;
		answer[FORMAT_IDX_NSF] = FORMAT_NUMBER_SUPPORTING_FRAGMENTS;
		answer[FORMAT_IDX_NSR] = FORMAT_NUMBER_SPLIT_READS;
		answer[FORMAT_IDX_RNC] = FORMAT_REAL_NUMBER_COPIES;
		answer[FORMAT_IDX_NTADF] = FORMAT_NUMBER_TANDEM_DUP_FRAGS;
		answer[FORMAT_IDX_NTRDF] = FORMAT_NUMBER_TRANS_DUP_FRAGS;
		answer[FORMAT_IDX_TGEN] = FORMAT_TEXT_GENOTYPE;
		return answer;
	}
	
	private GenomicVariant variant;
	private Set<String> filters = new TreeSet<String>();
	private Map<String,GenomicVariantAnnotation> infoFields = new LinkedHashMap<String, GenomicVariantAnnotation>();
	private List<CalledGenomicVariant> calls;
	private int [] fieldsFormat;
	private VCFFileHeader header;
	public VCFRecord(GenomicVariant variant, List<String> filters,List<GenomicVariantAnnotation> infoFields, int [] format, List<CalledGenomicVariant> calls, VCFFileHeader header) {
		this.variant = variant;
		this.filters.addAll(filters);
		for(GenomicVariantAnnotation ann:infoFields) addAnnotation(ann);
		fieldsFormat = Arrays.copyOf(format, format.length);
		this.calls = calls;
		this.header = header;
	}

	public VCFRecord(GenomicVariant variant, int [] format, List<CalledGenomicVariant> calls, VCFFileHeader header) {
		super();
		this.variant = variant;
		fieldsFormat = Arrays.copyOf(format, format.length);
		this.calls = calls;
		this.header = header;
	}
	
	public VCFRecord(GenomicVariant variant, int [] format, CalledGenomicVariant call, VCFFileHeader header) {
		super();
		this.variant = variant;
		fieldsFormat = Arrays.copyOf(format, format.length);
		this.calls = new ArrayList<CalledGenomicVariant>();
		this.calls.add(call);
		this.header = header;
	}
	
	
	public String getSequenceName() {
		return variant.getSequenceName();
	}

	public int getFirst() {
		return variant.getFirst();
	}

	public int getLast() {
		return variant.getLast();
	}

	public int length() {
		return variant.length();
	}
	

	public boolean isPositiveStrand() {
		return variant.isPositiveStrand();
	}

	public boolean isNegativeStrand() {
		return variant.isNegativeStrand();
	}

	public GenomicVariant getVariant() {
		return variant;
	}
	public List<String> getFilters() {
		return new ArrayList<String>(filters);
	}
	public List<GenomicVariantAnnotation> getInfoFields() {
		return new ArrayList<GenomicVariantAnnotation>(infoFields.values());
	}
	public GenomicVariantAnnotation getInfoField(String key) {
		return infoFields.get(key);
	}
	public List<CalledGenomicVariant> getCalls() {
		return calls;
	}
	public void addFilter(String filter) {
		filters.add(filter);
	}
	public void removeFilter(String filter) {
		filters.remove(filter);
	}
	public void addAnnotation(GenomicVariantAnnotation annotation) {
		infoFields.put(annotation.getAttribute(), annotation);
	}
	public void removeAnnotation (String key) {
		infoFields.remove(key);
	}
	public int[] getFieldsFormat() {
		return fieldsFormat;
	}
	public void setFieldsFormat(int[] fieldsFormat) {
		this.fieldsFormat = fieldsFormat;
	}
	public VCFFileHeader getHeader() {
		return header;
	}
	public void setHeader(VCFFileHeader header) {
		this.header = header;
	}
	public VariantFunctionalAnnotation getNGSEPFunctionalAnnotation () {
		VariantFunctionalAnnotation answer = null;
		GenomicVariantAnnotation ann = getInfoField(GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ANNOTATION);
		if(ann == null) return null;
		VariantFunctionalAnnotationType type = VariantFunctionalAnnotationType.getTypeBySearchKey(ann.getValue().toString());
		if(type == null) return null;
		answer = new VariantFunctionalAnnotation(variant, type);
		return answer;
	}
	public static VCFRecord createDefaultPopulationVCFRecord(GenomicVariant variant, List<CalledGenomicVariant> calls, VCFFileHeader header) {
		DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(calls, false);
		int [] format = variant.isSNV()?VCFRecord.DEF_FORMAT_ARRAY_NGSEP_SNV:VCFRecord.DEF_FORMAT_ARRAY_NGSEP_NOSNV;
		VCFRecord record = new VCFRecord(variant, format, calls, header);
		record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_SAMPLES_GENOTYPED, divStats.getNumSamplesGenotyped()));
		record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_NUMBER_ALLELES, divStats.getNumCalledAlleles()));
		record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_ALLELE_FREQUENCY_SPECTRUM, format(divStats.getAlleleCounts())));
		if(variant.isBiallelic()) record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_MAF, divStats.getMaf()));
		return record;
	}
	private static String format(int[] alleleCounts) {
		StringBuilder answer = new StringBuilder(""+alleleCounts[0]);
		for(int i=1;i<alleleCounts.length;i++) answer.append(","+alleleCounts[i]);
		return answer.toString();
	}
}
