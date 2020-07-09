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

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.main.io.ParseUtils;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledCNV;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.variants.VariantCallReport;

public class VCFFileReader implements Iterable<VCFRecord>,Closeable {
	private Logger log = Logger.getLogger(VCFFileReader.class.getName());
	
	public static final String NO_INFO_CHAR = ".";
	public static final int LOAD_MODE_CALLINFO = 0;
	public static final int LOAD_MODE_QUALITY = 1;
	public static final int LOAD_MODE_COPY_NUMBER = 2;
	public static final int LOAD_MODE_MINIMAL = 3;
	
	
	
	private BufferedReader in;
	private QualifiedSequenceList sequences = new QualifiedSequenceList();
	private VCFFileHeader header = new VCFFileHeader();
	
	private VCFFileIterator currentIterator = null;
	
	private int loadMode = LOAD_MODE_CALLINFO;
	
	public VCFFileReader (String filename) throws IOException {
		init(null,new File(filename));
	}
	public VCFFileReader (File file) throws IOException {
		init(null,file);
	}
	public VCFFileReader (InputStream stream) throws IOException {
		init(stream,null);
	}
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		if (log == null) throw new NullPointerException("Log can not be null");
		this.log = log;
	}
	
	
	public int getLoadMode() {
		return loadMode;
	}
	public void setLoadMode(int loadMode) {
		this.loadMode = loadMode;
	}
	
	public QualifiedSequenceList getSequences() {
		return sequences;
	}
	public void setSequences(QualifiedSequenceList sequences) {
		this.sequences = sequences;
	}
	public List<String> getSampleIds() {
		return header.getSampleIds();
	}
	

	@Override
	public void close() throws IOException {
		in.close();		
	}
	
	

	@Override
	public Iterator<VCFRecord> iterator() {
		if (in == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (currentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        currentIterator = new VCFFileIterator(); 
		return currentIterator;
	}
	
	private void init (InputStream stream, File file) throws IOException {
		if (stream != null && file != null) throw new IllegalArgumentException("Stream and file are mutually exclusive");
		if(file!=null) {
			stream = new FileInputStream(file);
			if(file.getName().endsWith(".gz")) {
				stream = new ConcatGZIPInputStream(stream);
			}
		}
		in = new BufferedReader(new InputStreamReader(stream));
		String samplesLine = loadHeader();
		header.loadSampleIds(samplesLine);
	}
	
	private String loadHeader() throws IOException {
		String line = in.readLine();
		while(line!=null && line.startsWith("##")) {
			header.loadHeaderLine(line);
			line = in.readLine(); 
		}
		return line;
	}
	
	
	private VCFRecord loadVCFRecord (String line) {
		String [] items = ParseUtils.parseString(line,'\t');
		if(items.length<8) {
			log.severe("Could not load line: "+line+". VCF records must have at least 8 columns");
			return null;
		}
		GenomicVariant variant = loadGenomicVariant(items);
		if(variant == null) return null;
		List<String> filters = loadFilters(items[6]);
		List<GenomicVariantAnnotation> infoFields = loadInfoField(variant, items[7]);
		
		//if(variant.getType() == GenomicVariant.TYPE_UNDETERMINED && items[3].length()>1) variant.setType(GenomicVariant.TYPE_INDEL);
		List<CalledGenomicVariant> calls = new ArrayList<CalledGenomicVariant>();
		List<Sample> samples = header.getSamples();
		
		if(items.length==8) {
			if (samples.size()>0) {
				log.severe("Can not load genomic variant at "+items[0]+":"+items[1]+". Number of genotyped samples does not coincide with number of samples in the header");
				return null;
			}
			return new VCFRecord(variant, filters, infoFields, new int [0], calls, header);
		}
		//If genotype information is present
		int[] formatInput = loadInputFormat(items[8]);
		if(items.length-9!=samples.size()) {
			log.severe("Can not load genomic variant at "+items[0]+":"+items[1]+". Number of genotyped samples does not coincide with number of samples in the header");
			return null;
		}
		boolean nonDefaultCN = false;
		for(int i=9;i<items.length;i++) {
			String [] itemsSample = ParseUtils.parseString(items[i], ':');
			Sample s = samples.get(i-9);
			CalledGenomicVariant call = loadCalledVariant(variant,formatInput,itemsSample,s); 
			calls.add(call);
			if(call.getCopyNumber()!=CalledGenomicVariant.DEFAULT_PLOIDY) nonDefaultCN = true;
		}
		int [] formatLoad = makeLoadFormat (formatInput,loadMode!=LOAD_MODE_MINIMAL && nonDefaultCN && variant.getType()<=GenomicVariant.TYPE_STR);	
		return new VCFRecord(variant, filters, infoFields, formatLoad, calls, header);
	}
	
	private GenomicVariant loadGenomicVariant(String[] items) {
		QualifiedSequence seq;
		try {
			seq = sequences.addOrLookupName(items[0]);
		} catch (RuntimeException e) {
			log.severe("Can not load genomic variant at "+items[0]+":"+items[1]+". Unrecognized sequence name. "+e.getMessage());
			return null;
		}
		int position;
		try {
			position = Integer.parseInt(items[1]);
			if(position<=0) throw new NumberFormatException("Negative position: "+position);
		} catch(NumberFormatException e) {
			log.severe("Can not load genomic variant at "+items[0]+":"+items[1]+". Position must be a positive number");
			return null;
		}
		String id = null;
		if(!NO_INFO_CHAR.equals(items[2])) {
			id = items[2];
		}
		if(items[3].length()==0) {
			log.severe("Can not load genomic variant at "+items[0]+":"+items[1]+". The field with the reference allele is empty");
			return null;
		}
		if(items[4].length()==0) {
			log.severe("Can not load genomic variant at "+items[0]+":"+items[1]+". The field with alternative alleles is empty");
			return null;
		}
		List<String> alleles = new ArrayList<String>();
		alleles.add(items[3]);
		String [] altAlleles = ParseUtils.parseString(items[4], ',');
		
		if(items[4].charAt(0)!='.') alleles.addAll(Arrays.asList(altAlleles));
		if(alleles.size()>GenomicVariant.MAX_NUM_ALLELES) {
			log.severe("Can not load genomic variant at "+items[0]+":"+items[1]+". Number of alleles "+alleles.size()+" is greater than the maximum allowed "+GenomicVariant.MAX_NUM_ALLELES);
			return null;
		}
		short variantQS=0;
		if(items[5].length()> 0 && items[5].charAt(0)!='.') {
			try {
				//It is loaded as double for format compatibility. However it is treated internally as a short
				double qsD = Double.parseDouble(items[5]);
				if(qsD > Short.MAX_VALUE) qsD = 255;
				variantQS = (short) Math.round(qsD);
				if(variantQS<0) throw new NumberFormatException("Negative variant QS: "+variantQS);
			} catch (NumberFormatException e) {
				log.severe ("Can not load genomic variant at "+items[0]+":"+items[1]+". Quality score must be a non negative number");
				return null;
			}
		}
		char c1 = items[3].charAt(0);
		char c2 = items[4].charAt(0);
		if(alleles.size()==2 && items[3].length()==1 && items[4].length()==1 && DNASequence.isInAlphabeth(c1) && DNASequence.isInAlphabeth(c2)) {
			SNV snv = new SNV(seq.getName(), position, c1, c2);
			snv.setId(id);
			snv.setVariantQS(variantQS);
			return snv;
		} else {
			GenomicVariantImpl variant = new GenomicVariantImpl(seq.getName(), position, alleles);
			variant.setId(id);
			variant.setVariantQS(variantQS);
			return variant;
		}
	}
	private List<String> loadFilters(String filtersField) {
		//TODO: Improve memory usage
		if(NO_INFO_CHAR.equals(filtersField)) return new ArrayList<String>();
 		return Arrays.asList(ParseUtils.parseString(filtersField, ';'));
	}
	private List<GenomicVariantAnnotation> loadInfoField (GenomicVariant variant,String infoField) {
		List<GenomicVariantAnnotation> annotations = new ArrayList<GenomicVariantAnnotation>();
		if(NO_INFO_CHAR.equals(infoField)) return annotations;
		String [] infoItems = ParseUtils.parseStringWithText(infoField, ';','\"');
		for(int i=0;i<infoItems.length;i++) {
			int idx = infoItems[i].indexOf("=");
			if(idx <0) annotations.add(new GenomicVariantAnnotation(variant, infoItems[i], true));
			else {
				String attribute = infoItems[i].substring(0,idx);
				String value = infoItems[i].substring(idx+1);
				GenomicVariantAnnotation ann = new GenomicVariantAnnotation(variant, attribute, value); 
				annotations.add(ann);
				if(GenomicVariantAnnotation.ATTRIBUTE_TYPE.equals(attribute)) {
					byte type = GenomicVariantImpl.getVariantTypeId(value);
					if(type>0 && type<=GenomicVariant.TYPE_INVERSION) {
						variant.setType(type);
					} /*else {
						log.severe("Can not load type for genomic variant at "+variant.getSequenceName()+":"+variant.getFirst()+". Invalid type: "+value);
					}*/
				}
				if(GenomicVariantAnnotation.ATTRIBUTE_END.equals(attribute) && (variant instanceof GenomicVariantImpl)) {
					//Imprecise variant
					int last = Integer.parseInt(value);
					GenomicVariantImpl impl = (GenomicVariantImpl)variant; 
					impl.setLast(last);
					//TODO: only change length if the length attribute is not present
					impl.setLength(last-variant.getFirst()+1);
				}
			}
		}
		return annotations;
	}
	
	private int[] loadInputFormat(String formatStr) {
		String [] itemsFormat = ParseUtils.parseString(formatStr, ':');
		int [] answer = new int [itemsFormat.length];
		Arrays.fill(answer, -1);
		for(int i=0;i<itemsFormat.length;i++) {
			Integer index = VCFRecord.KNOWN_FORMAT_FIELDS_MAP.get(itemsFormat[i]);
			if(index!=null) {
				if(loadMode == LOAD_MODE_MINIMAL && index != VCFRecord.FORMAT_IDX_GT) continue;
				if(loadMode == LOAD_MODE_COPY_NUMBER && index != VCFRecord.FORMAT_IDX_GT && index != VCFRecord.FORMAT_IDX_ACN) continue;
				if(loadMode == LOAD_MODE_QUALITY && index != VCFRecord.FORMAT_IDX_GT && index != VCFRecord.FORMAT_IDX_ACN && index != VCFRecord.FORMAT_IDX_GQ) continue;
				answer[i] = index;
			}
		}
		return answer;
	}
	private CalledGenomicVariant loadCalledVariant(GenomicVariant variant,int [] format, String[] itemsSample, Sample sample) {
		String sampleId = sample.getId();
		if(itemsSample.length>format.length) {
			log.severe("Can not load genotype of sample "+sampleId+" for genomic variant at "+variant.getSequenceName()+":"+variant.getFirst()+". Sample information does not match format");
			CalledGenomicVariantImpl answer = new CalledGenomicVariantImpl(variant,new byte[0]);
			answer.setSampleId(sampleId);
			answer.updateAllelesCopyNumberFromCounts(sample.getNormalPloidy());
			return answer;
		}
		String [] knownItemsSample = new String [VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY.length];
		Arrays.fill(knownItemsSample, null);
		for(int i=0;i<itemsSample.length;i++) {
			if(format[i]>=0)knownItemsSample[format[i]] = itemsSample[i];
		}
		//Load genotype field
		String [] alleles = variant.getAlleles();
		int numAlleles = alleles.length;
		//Load called alleles with copy number if present
		String genotypeStr=knownItemsSample[VCFRecord.FORMAT_IDX_GT];
		if(genotypeStr==null) genotypeStr = ".";
		String [] callItems = ParseUtils.parseString(genotypeStr, '|', '/');
		boolean phased = callItems.length>1 && genotypeStr.charAt(callItems[0].length())=='|';
		short [] allelesCNG = new short[numAlleles];
		byte [] phasedAlleles = new byte[callItems.length];
		short totalCNG = (short)Math.min(CalledGenomicVariant.MAX_PLOIDY_SAMPLE, callItems.length);
		Arrays.fill(allelesCNG, (short)0);
		Set<Byte> uniqueAlleles = new TreeSet<Byte>();
		for(int j=0;j<callItems.length;j++)  {
			if(callItems[j].length()>0 && callItems[j].charAt(0)!='.') {
				byte nextAlleleId;
				try {
					nextAlleleId = Byte.parseByte(callItems[j]);
				} catch (NumberFormatException e) {
					log.severe("Can not load genotype of sample "+sampleId+" for genomic variant at "+variant.getSequenceName()+":"+variant.getFirst()+". Called allele "+callItems[j]+" is not a number");
					uniqueAlleles.clear();
					break;
				}
				if(nextAlleleId<0 || nextAlleleId>=numAlleles) {
					log.severe("Can not load genotype of sample "+sampleId+" for genomic variant at "+variant.getSequenceName()+":"+variant.getFirst()+". Inconsistent called allele "+nextAlleleId+" for the total number of alleles: "+numAlleles);
					uniqueAlleles.clear();
					break;
				}
				uniqueAlleles.add(nextAlleleId);
				allelesCNG[nextAlleleId]++;
				phasedAlleles [j] = nextAlleleId;
			}
		}
		byte [] calledAlleleIds = new byte[uniqueAlleles.size()];
		if(calledAlleleIds.length==0) phased = false;
		Iterator<Byte> it = uniqueAlleles.iterator();
		for(int j=0;it.hasNext();j++)  {
			byte nextAlleleId = it.next();
			calledAlleleIds[j] = nextAlleleId;
		}
		
		if(calledAlleleIds.length>numAlleles) {
			log.severe("Can not load genotype of sample "+sampleId+" for genomic variant at "+variant.getSequenceName()+":"+variant.getFirst()+". More called alleles than total alleles");
			calledAlleleIds = new byte[0];
		}
		//Load variant-specific optional information 
		int [] allCounts = loadCounts(knownItemsSample[VCFRecord.FORMAT_IDX_BSDP],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_BSDP],sampleId,variant,4);
		int [] counts = loadCounts(knownItemsSample[VCFRecord.FORMAT_IDX_ADP],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_ADP],sampleId,variant,alleles.length);
		double [][] logConditionals = loadConditionals (numAlleles,knownItemsSample[VCFRecord.FORMAT_IDX_PL],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_PL],sampleId,variant,callItems.length==1,true);
		if(logConditionals==null) logConditionals = loadConditionals (numAlleles,knownItemsSample[VCFRecord.FORMAT_IDX_GL],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_GL],sampleId,variant,callItems.length==1,false);
		//Create object consistent with the variant information
		CalledGenomicVariant answer = null;
		if(variant instanceof SNV) {
			byte genotype = CalledSNV.GENOTYPE_HOMOREF;
			if(calledAlleleIds.length==0) genotype = CalledSNV.GENOTYPE_UNDECIDED;
			else if(calledAlleleIds.length>1) genotype=CalledSNV.GENOTYPE_HETERO;
			else if (calledAlleleIds[0]>0) genotype = CalledSNV.GENOTYPE_HOMOALT;
			CalledSNV snv = new CalledSNV((SNV)variant, genotype);
			answer = snv;
			if(allCounts!=null) snv.setAllBaseCounts(allCounts);
			else if (counts!=null && counts.length==2) {
				snv.setCountReference(counts[0]);
				snv.setCountAlternative(counts[1]);
			}
			if(logConditionals!=null) snv.setRefAltGenotypeLogConditionals(logConditionals);
		} else if (variant.getType() == GenomicVariant.TYPE_CNV) {
			CalledCNV cnv;
			if (calledAlleleIds.length==1) {
				//If the genotype field is not empty, then it is interpreted as the number of copies
				cnv = new CalledCNV(variant, (byte)calledAlleleIds[0]);
			} else {
				cnv = new CalledCNV(variant);
			}
			answer = cnv;
			Double v = loadSingleNumber(knownItemsSample[VCFRecord.FORMAT_IDX_RNC],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_RNC],sampleId,variant, false);
			if(v!=null) cnv.setNumCopies(v.floatValue(),false);
			v = loadSingleNumber(knownItemsSample[VCFRecord.FORMAT_IDX_NTADF],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_NTADF],sampleId,variant, true);
			if(v!=null) cnv.setTandemFragments(v.intValue());
			v = loadSingleNumber(knownItemsSample[VCFRecord.FORMAT_IDX_NTRDF],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_NTRDF],sampleId,variant, true);
			if(v!=null) cnv.setTransDupFragments(v.intValue());
			String textGen = knownItemsSample[VCFRecord.FORMAT_IDX_TGEN];
			if(textGen!=null) cnv.setTextGenotype(textGen);
		} else {
			CalledGenomicVariantImpl cv = new CalledGenomicVariantImpl(variant, calledAlleleIds);
			answer = cv;
			if(allCounts!=null) cv.setAllCounts(allCounts);
			cv.setCallReport(new VariantCallReport(alleles, counts, logConditionals));
		}
		answer.setSampleId(sampleId);
		//Allowed real numbers in this field to be able to load freebayes GQ fields
		Double v = loadSingleNumber(knownItemsSample[VCFRecord.FORMAT_IDX_GQ],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_GQ],sampleId,variant, false);
		if(v!=null) answer.setGenotypeQuality(v.shortValue());
		v = loadSingleNumber(knownItemsSample[VCFRecord.FORMAT_IDX_DP],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_DP],sampleId,variant, true);
		if(v!=null) answer.setTotalReadDepth(v.intValue());
		//Load alleles copy number
		if(variant.getType() != GenomicVariant.TYPE_CNV) {
			short [] allelesCN = null;
			int [] countsA = loadCounts(knownItemsSample[VCFRecord.FORMAT_IDX_ACN],VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[VCFRecord.FORMAT_IDX_ACN],sampleId,variant,alleles.length);
			int totalCopyNumber = 0;
			if(countsA!=null) {
				allelesCN = new short [countsA.length];		
				for(int j=0;j<countsA.length;j++) {
					allelesCN[j] = (short)countsA[j];
					totalCopyNumber += countsA[j];
				}
				if(totalCopyNumber>CalledGenomicVariant.MAX_PLOIDY_SAMPLE) {
					log.severe("Can not load alleles copy number for sample "+sampleId+" at genomic variant "+variant.getSequenceName()+":"+variant.getFirst()+". Total copy number is larger than the maximum allowed value "+CalledGenomicVariant.MAX_PLOIDY_SAMPLE);
					allelesCN = null;
					totalCopyNumber = 0;
				}
			}
			if(allelesCN==null) {	
				allelesCN = allelesCNG;
				totalCopyNumber = totalCNG;
			}
			try {
				if(answer.isUndecided()) answer.updateAllelesCopyNumberFromCounts((short)totalCopyNumber);
				else answer.setAllelesCopyNumber(allelesCN);
			} catch (IllegalArgumentException e) {
				log.severe("Can not load alleles copy number for sample "+sampleId+" at genomic variant "+variant.getSequenceName()+":"+variant.getFirst()+". "+e.getMessage());
				answer.updateAllelesCopyNumberFromCounts((short)totalCopyNumber);
			}
		}
		//Load phasing
		if(phased) {
			int copyNumber = answer.getCopyNumber();
			if(copyNumber == phasedAlleles.length) {
				//TODO: Load SNVs with  
				if(answer instanceof CalledSNV && copyNumber == 2) {
					CalledSNV csnv = (CalledSNV)answer;
					csnv.setPhasingCN2(phasedAlleles[0]==1);
				} else if (answer instanceof CalledGenomicVariantImpl) {
					CalledGenomicVariantImpl call = (CalledGenomicVariantImpl)answer;
					call.setIndexesPhasedAlleles(phasedAlleles);
				} else {
					log.severe("Can not load phasing information for sample "+sampleId+" at genomic variant "+variant.getSequenceName()+":"+variant.getFirst()+". Phasing of SNVs with high copy number still not supported");
				}
				
			} else {
				log.severe("Can not load phasing information for sample "+sampleId+" at genomic variant "+variant.getSequenceName()+":"+variant.getFirst()+". Number of phased alleles inconsistent with total copy number");
			}
		}
		
		return answer;
		
	}
	private Double loadSingleNumber (String value, String formatField, String sampleId, GenomicVariant var, boolean integer) {
		if(value==null || NO_INFO_CHAR.equals(value)) return null;
		try {
			if(integer) {
				return 0.0+Integer.parseInt(value);
			}
			return Double.parseDouble(value);
		} catch (NumberFormatException e) {
			log.severe("Can not load value of format field "+formatField+" for sample "+sampleId+" at genomic variant at "+var.getSequenceName()+":"+var.getFirst()+". Error parsing value: "+value);
			return null;
		}
	}
	private int[] loadCounts(String countsStr, String formatField, String sampleId, GenomicVariant var, int expectedCounts) {
		if(countsStr==null || NO_INFO_CHAR.equals(countsStr)) return null; 
		String [] countsItems = ParseUtils.parseString(countsStr, ',');
		if(countsItems.length!=expectedCounts) {
			log.severe("Can not load counts of format field "+formatField+" for sample "+sampleId+" at genomic variant at "+var.getSequenceName()+":"+var.getFirst()+". Number of depths: "+countsItems.length+" different than the expected number: "+expectedCounts);
			return null;
		}
		int [] counts = new int[countsItems.length];
		for(int j=0;j<countsItems.length;j++) {
			try {
				counts[j] = Integer.parseInt(countsItems[j]);
			} catch (NumberFormatException e) {
				log.severe("Can not load counts of format field "+formatField+" for sample "+sampleId+" at genomic variant at "+var.getSequenceName()+":"+var.getFirst()+". Error parsing count: "+countsItems[j]);
				return null;
			}
		}
		return counts;
	}
	private double[][] loadConditionals(int numAlleles, String dataStr, String formatField, String sampleId, GenomicVariant var, boolean haploidGT, boolean phredScaled) {
		if(dataStr==null || NO_INFO_CHAR.equals(dataStr)) return null;
		double [][] answer = new double [numAlleles][numAlleles];
		String [] dataItems = ParseUtils.parseString(dataStr, ',');
		if(haploidGT && numAlleles==dataItems.length) {
			for(int i=0;i<dataItems.length;i++) {
				Arrays.fill(answer[i], -100);
				double next;
				try {
					if(phredScaled) {
						next = Integer.parseInt(dataItems[i]);
						next = -next/10;
					} else {
						next = Double.parseDouble(dataItems[i]);
					}
				} catch (NumberFormatException e) {
					log.severe("Can not load values of format field "+formatField+" for sample "+sampleId+" at genomic variant at "+var.getSequenceName()+":"+var.getFirst()+". Error parsing value: "+dataItems[i]);
					return null;
				}
				answer[i][i] = next;
			}
			return answer;
		}
		int k=0;
		for(int j=0;j<numAlleles;j++) {
			for(int i=0;i<=j;i++) {
				//Not enough fields found
				if(k>=dataItems.length) {
					log.severe("Can not load genotype data of format field "+formatField+" for sample "+sampleId+" at genomic variant at "+var.getSequenceName()+":"+var.getFirst()+". Unexpected number of values: "+dataItems.length+" for "+numAlleles+" alleles");
					return null;
				}
				double next;
				try {
					if(phredScaled) {
						next = Integer.parseInt(dataItems[k]);
						next = -next/10;
					} else {
						next = Double.parseDouble(dataItems[k]);
					}
				} catch (NumberFormatException e) {
					log.severe("Can not load values of format field "+formatField+" for sample "+sampleId+" at genomic variant at "+var.getSequenceName()+":"+var.getFirst()+". Error parsing value: "+dataItems[k]);
					return null;
				}
				answer[i][j] = next;
				if(i!=j) answer [j][i] = next;
				k++;
			}
		}
		return answer;
	}
	
	private int[] makeLoadFormat(int[] formatInput, boolean forceCN) {
		List<Integer> loadFields = new ArrayList<Integer>();
		for(int i=0;i<formatInput.length;i++) {
			int idx = formatInput[i];
			if(idx>=0) loadFields.add(idx);
		}
		//Force copy number if it was originally on the genotype information
		if(forceCN && !loadFields.contains(VCFRecord.FORMAT_IDX_ACN)) loadFields.add(VCFRecord.FORMAT_IDX_ACN);
		int [] answer = new int [loadFields.size()];
		for(int i=0;i<answer.length;i++) {
			answer[i] = loadFields.get(i);
		}
		return answer;
	}
	
	/**
	 * Loads only the variants basic information
	 * @param filename Name for the VCF file
	 * @return List<GenomicVariant> List of variants within the given file
	 * @throws IOException If the file can not be read
	 */
	public static List<GenomicVariant> loadVariants(String filename, boolean filterReferenceSitesGVCF) throws IOException {
		List<GenomicVariant> answer = new ArrayList<GenomicVariant>();
		
		try (VCFFileReader in = new VCFFileReader(filename)) {
			in.setLoadMode(LOAD_MODE_MINIMAL);
			Iterator<VCFRecord> it = in.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				GenomicVariant var = record.getVariant();
				if(filterReferenceSitesGVCF && var.getAlleles().length<2) continue;
				answer.add(var);
			}
		}
		return answer;
	}

	/**
	 * Loads all calls for the VCF of one individual
	 * @param filename Name for the VCF file
	 * @return List<CAlledGenomicVariant> List of calls within the given file
	 * @throws IOException If the file can not be read
	 */
	public static List<CalledGenomicVariant> loadCalledVariantsSingleIndividualVCF(String filename) throws IOException {
		List<CalledGenomicVariant> answer = new ArrayList<CalledGenomicVariant>();
		try (VCFFileReader in = new VCFFileReader(filename)) {
			Iterator<VCFRecord> it = in.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				List<CalledGenomicVariant> calls = record.getCalls();
				if(calls.size()>=1) answer.add(calls.get(0));
				if(answer.size()%100000==0) System.out.println("Loaded "+answer.size()+" calls");
			}
		}
		return answer;
	}

	
	
	public VCFFileHeader getHeader() {
		return header;
	}

	private class VCFFileIterator implements Iterator<VCFRecord> {
		private VCFRecord nextRecord;
		public VCFFileIterator() {
			nextRecord = loadRecord();
		}
		@Override
		public boolean hasNext() {
			return nextRecord!=null;
		}

		@Override
		public VCFRecord next() {
			if(nextRecord==null) throw new NoSuchElementException();
			VCFRecord answer = nextRecord;
			nextRecord = loadRecord();
			return answer;
		}

		private VCFRecord loadRecord() {
			String line;
			while(true) {
				try {
					line = in.readLine();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
				if(line==null) return null;
				VCFRecord answer = loadVCFRecord(line);
				if(answer !=null) return answer;
			} 
		}
		@Override
		public void remove() {
			throw new UnsupportedOperationException("Remove not supported by VCFFileIterator");
		}
	}

}
