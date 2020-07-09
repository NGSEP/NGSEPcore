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

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.main.io.ParseUtils;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.Sample;

public class VCFFileHeader {
	private List<Sample> samples = new ArrayList<Sample>();
	private Map<String,Integer> samplesIndexMap = new TreeMap<String, Integer>();
	private Map<String, Sample> samplesWithHeaderLine = new TreeMap<String, Sample>();
	private List<VCFHeaderLine> idHeaderLines = new ArrayList<VCFHeaderLine>();
	private List<String> unstructuredHeaderLines = new ArrayList<String>();
	private String version;
	
	public VCFFileHeader () {
		
	}
	public VCFFileHeader cloneEmpty () {
		VCFFileHeader answer = new VCFFileHeader();
		answer.idHeaderLines = idHeaderLines;
		answer.unstructuredHeaderLines = unstructuredHeaderLines;
		answer.version = version;
		return answer;
	}
	public static VCFFileHeader makeDefaultEmptyHeader () {
		VCFFileHeader header = new VCFFileHeader();
		header.version = "VCFv4.2";
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_IN_CNV, "\"Number of samples with CNVs around this variant\"", "1", "Integer"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ANNOTATION, "\"Variant annotation based on a gene model\"", "1", "String"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ID, "\"Id of the transcript related to the variant annotation\"","1","String"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_GENE_NAME,"\"Name of the gene related to the variant annotation\"","1","String"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_CODON, "\"One based codon position of the start of the variant. The decimal is the codon position\"","1","Float"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_AMINOACID_CHANGE, "\"Description of the aminoacid change produced by a non-synonymous mutation. String encoded as reference aminoacid, position and mutated aminoacid\"","1","String"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_SAMPLES_GENOTYPED, "\"Number of samples genotyped\"","1","Integer"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_MAF,"\"Minor allele frequency\"","1","Float"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_NUMBER_ALLELES,"\"Number of alleles in called genotypes\"","1","Integer"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_ALLELE_FREQUENCY_SPECTRUM,"\"Allele counts over the population for all alleles, including the reference\"","R","Integer"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_TYPE,"\"Type of variant\"", "1", "String"));
		header.idHeaderLines.add(new VCFHeaderLine("INFO", GenomicVariantAnnotation.ATTRIBUTE_FISHER_STRAND_BIAS,"\"Phred-scaled p-value using Fisher's exact test to detect strand bias\"","1","Float"));
		header.idHeaderLines.add(new VCFHeaderLine("FORMAT", VCFRecord.FORMAT_GENOTYPE,"\"Genotype\"","1","String"));
		header.idHeaderLines.add(new VCFHeaderLine("FORMAT", VCFRecord.FORMAT_GENOTYPE_PHRED_LIKELIHOOD,"\"Phred-scaled genotype likelihoods rounded to the closest integer\"","G","Integer"));
		header.idHeaderLines.add(new VCFHeaderLine("FORMAT", VCFRecord.FORMAT_GENOTYPE_QUALITY,"\"Genotype quality\"","1","Integer"));
		header.idHeaderLines.add(new VCFHeaderLine("FORMAT", VCFRecord.FORMAT_DEPTH,"\"Read depth\"","1","Integer"));
		header.idHeaderLines.add(new VCFHeaderLine("FORMAT", VCFRecord.FORMAT_ALLELE_DEPTH,"\"Counts for observed alleles, including the reference allele\"","R","Integer"));
		header.idHeaderLines.add(new VCFHeaderLine("FORMAT", VCFRecord.FORMAT_ALL_BASES_SNP_DEPTH,"\"Number of base calls (depth) for the 4 nucleotides in called SNVs sorted as A,C,G,T\"","4","Integer"));
		header.idHeaderLines.add(new VCFHeaderLine("FORMAT", VCFRecord.FORMAT_ALLELE_COPY_NUMBER,"\"Predicted copy number of each allele taking into account the prediction of number of copies of the region surrounding the variant\"","R","Integer"));
		return header;
	}
	public void loadHeaderLine(String line) throws IOException {
		if(line.startsWith("##SAMPLE=")) {
			loadSampleHeader(line.substring(10,line.length()-1));
		} else if(line.startsWith("##fileformat=")){
			version = line.substring(13);
		} else {
			int equalI = line.indexOf("=<ID=");
			if(equalI>0) loadIdHeaderLine(line.substring(2,equalI),line.substring(equalI+2,line.length()-1));
			else unstructuredHeaderLines.add(line.substring(2));
		}
		
	}
	private void loadIdHeaderLine(String headerType, String data) throws IOException {
		String [] items = ParseUtils.parseStringWithText(data, ',', '\"');
		VCFHeaderLine line = null;
		for(int i=0;i<items.length;i++) {
			int eqIdx = items[i].indexOf('=');
			if(eqIdx<=0 || eqIdx ==items[i].length()-1) throw new IOException("Can not load header line: "+data+". Wrong format for item: "+items[i]);
			String property = items[i].substring(0, eqIdx);
			String value = items[i].substring(eqIdx+1);
			if(i==0) {
				//Id
				line = new VCFHeaderLine(headerType, value);
			} else if ("Type".equals(property)) {
				line.setDataType(value);
			} else if ("Number".equals(property)) {
				line.setNumber(value);
			} else if ("Description".equals(property)) {
				line.setDescription(value);
			}
			line.setAttribute(property, value);
		}
		idHeaderLines.add(line);
	}
	private void loadSampleHeader(String line) {
		String [] items = line.split("=|,");
		String id = "";
		byte ploidy = GenomicVariant.DEFAULT_PLOIDY;
		for(int i=0;i<items.length;i++) {
			if("ID".equals(items[i]) && i<items.length-1) {
				id = items[i+1];
			} else if ("PL".equals(items[i]) && i<items.length-1) {
				ploidy = Byte.parseByte(items[i+1]);
			}
		}
		if(id!=null) {
			Sample s = new Sample(id);
			if(ploidy>0) s.setNormalPloidy(ploidy);
			samplesWithHeaderLine.put(id, s);
		}	
	}
	
	public void loadSampleIds(String samplesLine) throws IOException {
		samples = new ArrayList<Sample>();
		if(samplesLine!=null && samplesLine.startsWith("#")) {
			String [] items = ParseUtils.parseString(samplesLine, '\t');
			for(int i=9;i<items.length;i++) {
				String sampleId = items[i];
				Sample s = samplesWithHeaderLine.get(sampleId);
				if(s == null) s = new Sample(sampleId);
				samples.add(s);
			}
		} else {
			throw new IOException("VCF file does not have line with sample ids");
		}
		samplesIndexMap = new TreeMap<String, Integer>();
		for(int i=0;i<samples.size();i++) {
			String sampleId = samples.get(i).getId(); 
			samplesIndexMap.put(sampleId, i);
		}
	}
	
	public void addSamples (List<Sample> samples) {
		for(Sample sample:samples) {
			addSample(sample, false);
		}
	}
	public void addSample(Sample s, boolean headerLine) {
		samples.add(s);
		samplesIndexMap.put(s.getId(), samples.size()-1);
		if(headerLine) samplesWithHeaderLine.put(s.getId(), s);
	}
	
	public void addDefaultSample(String sampleId) {
		addSample(new Sample(sampleId), false);
	}
	
	public List<Sample> getSamples() {
		return samples;
	}
	
	public int getIndexSampleId(String sampleId) {
		Integer i = samplesIndexMap.get(sampleId);
		if(i!=null) return i;
		return -1;
	}
	
	public int[] calculateIndexes(String [] selectedSampleIds) {
		int [] indexes = new int[selectedSampleIds.length];
		Arrays.fill(indexes, -1);
		for(int i=0;i<selectedSampleIds.length;i++) {
			indexes[i] = getIndexSampleId(selectedSampleIds[i]);
		}
		return indexes;
	}
	public Map<String, Sample> getSamplesWithHeaderLine() {
		return samplesWithHeaderLine;
	}
	public List<String> getSampleIds() {
		List<String> sampleIds = new ArrayList<String>();
		for(Sample s:samples) sampleIds.add(s.getId());
		return sampleIds;
	}
	public void print(PrintStream out) {
		out.println("##fileformat="+version);
		for(VCFHeaderLine line:idHeaderLines) {
			out.print("##"+line.getHeaderType()+"=<");
			int i=0;
			for(String property:line.getAttributes().keySet()) {
				String value = line.getAttribute(property);
				if(i>0) out.print(",");
				out.print(property+"="+value);
				i++;
			}
			out.println(">");
		}
		for(String line:unstructuredHeaderLines) {
			out.println("##"+line);
		}
		for(Sample sample:samplesWithHeaderLine.values()) out.println("##SAMPLE=<ID="+sample.getId()+",PL="+sample.getNormalPloidy()+">");
		printColumnNamesLine(out);
	}
	private void printColumnNamesLine(PrintStream out) {
		out.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		if(samples.size()>0) out.print("\tFORMAT");
		for(int i=0;i<samples.size();i++) {
			out.print("\t"+samples.get(i).getId());
		}
		out.println();
	}
}
