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

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import ngsep.main.io.ParseUtils;
import ngsep.variants.CalledCNV;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.Sample;
import ngsep.variants.VariantCallReport;

public class VCFFileWriter {
	
	public void printVCFRecords (List<VCFRecord> records, PrintStream out) {
		for(VCFRecord record:records) {
			printVCFRecord(record, out);
		}
	}
	public void printVCFRecord (VCFRecord record, PrintStream out) {
		GenomicVariant var = record.getVariant();
		//Add type as annotation if still not added
		byte type = var.getType();
		String typeName = GenomicVariantImpl.getVariantTypeName(var.getType());
		if(type!=GenomicVariant.TYPE_UNDETERMINED && type!=GenomicVariant.TYPE_BIALLELIC_SNV && typeName!=null) record.addAnnotation(new GenomicVariantAnnotation(var, GenomicVariantAnnotation.ATTRIBUTE_TYPE, typeName));
		printBasicVariantInfo(var, out);
		printFilters(record.getFilters(),out);
		printInfoField(record.getInfoFields(), out);
		List<CalledGenomicVariant> calls = record.getCalls();
		if(calls.size()>0) {
			int [] outFormat = record.getFieldsFormat();
			printGenotypeFormat(out,outFormat);
			//Genotype
			List<Sample> samples = null;
			if(record.getHeader()!=null) samples = record.getHeader().getSamples();
			for(int i=0;i<calls.size();i++) {
				short ploidy = GenomicVariant.DEFAULT_PLOIDY;
				if(samples!=null) ploidy = samples.get(i).getNormalPloidy();
				printGenotypeInfo(calls.get(i), out, outFormat,ploidy);
			}
		}
		
		out.println();
	}
	
	private void printFilters(List<String> filters, PrintStream out) {
		out.print("\t");
		if(filters==null || filters.size()==0) {
			out.print(VCFFileReader.NO_INFO_CHAR);
			return;
		}
		boolean printed = false;
		for(String filter:filters) {
			if(printed) out.print(";");
			printed = true;
			out.print(filter);
		}
	}
	private void printBasicVariantInfo(GenomicVariant var,PrintStream out) {
		out.print(var.getSequenceName()+"\t");
		out.print(var.getFirst()+"\t");
		String id = var.getId();
		if(id==null) {
			id = VCFFileReader.NO_INFO_CHAR;
		}
		out.print(id+"\t");
		
		String [] alleles = var.getAlleles();
		out.print(alleles[0]+"\t");
		if(alleles.length==1) out.print(VCFFileReader.NO_INFO_CHAR);
		else {
			//Starts at 1 to ignore the reference allele
			for(int i=1;i<alleles.length;i++) {
				if(i>1) out.print(",");
				out.print(alleles[i]);
			}
		}
		out.print("\t");
		out.print(var.getVariantQS());
	}
	private void printInfoField(List<GenomicVariantAnnotation> info, PrintStream out) {
		out.print("\t");
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
		boolean printed = false;
		for(GenomicVariantAnnotation ann:info) {
			Object value = ann.getValue();
			if(value == null) continue;
			if(value instanceof Boolean) {
				if((Boolean) value) {
					if(printed) out.print(";");
					printed = true;
					out.print(ann.getAttribute());
				}
			} else if (value instanceof Collection<?>) {
				Collection<?> values = (Collection<?>)value;
				if(values.size()>0) {
					if(printed) out.print(";");
					printed = true;
					out.print(ann.getAttribute()+"=");
				}
				Iterator<?> it = values.iterator();
				for(int i=0;it.hasNext();i++) {
					Object o = it.next();
					if(i>0) out.print(",");
					if(o instanceof Double) {
						out.print(fmt.format(o));
					} else {
						out.print(o.toString());
					}
					
				}
			} else if (value instanceof Double) {
				if(printed) out.print(";");
				printed = true;
				double valN = (Double)value;
				out.print(ann.getAttribute()+"="+fmt.format(valN));
			} else {
				if(printed) out.print(";");
				printed = true;
				out.print(ann.getAttribute()+"="+ann.getValue().toString());
			}
		}
		if(!printed) out.print(VCFFileReader.NO_INFO_CHAR);
	}
	private void printGenotypeFormat(PrintStream out, int [] format) {
		//Genotype format
		out.print("\t");
		for(int f=0;f<format.length;f++) {
			if(f>0) out.print(":");
			int formatIdx = format[f];
			String formatName = VCFRecord.KNOWN_FORMAT_FIELDS_ARRAY[formatIdx];
			out.print(formatName);
		}
	}
	private void printGenotypeInfo(CalledGenomicVariant var, PrintStream out, int [] format, short ploidy) {
		out.print("\t");
		String[] alleles = var.getAlleles();
		VariantCallReport report = var.getCallReport();
		CalledCNV cnv = null;
		if (var instanceof CalledCNV) {
			cnv = (CalledCNV) var;
		}
		byte [] idxsCalledAlleles = var.getIndexesCalledAlleles();
		for(int f=0;f<format.length;f++) {
			if(f>0) out.print(":");
			int formatIdx = format[f];
			if(formatIdx == VCFRecord.FORMAT_IDX_GT) {
				boolean phased = var.isPhased();
				if (idxsCalledAlleles.length == 0) {
					//Undecided call
					out.print(".");
					if(ploidy>1) out.print ("/."); 
					/*for(int i=1;i<ploidy;i++) {
						out.print ("/.");
					}*/
				} else if(idxsCalledAlleles.length == 1) {
					//Homozygous call
					int idAllele = idxsCalledAlleles[0];
					out.print(""+idAllele);
					if(phased) {
						for(int i=1;i<ploidy;i++) {
							out.print ("|"+idAllele);
						}
					} else if(ploidy>1) out.print ("/"+idAllele);
				} else {
					//Heterozygous call
					byte [] finalAlleles = idxsCalledAlleles;
					if(phased) finalAlleles = var.getIndexesPhasedAlleles();
					for(int i=0;i<finalAlleles.length;i++) {
						//Since v2.1.4, alleles are not explicitly written with copy number anymore. Allele copy numbers are saved in the new format field Local Allele Copy Numbers (ACN)
						int idAllele = finalAlleles[i];
						if(i>0) out.print((phased?"|":"/"));
						out.print(""+idAllele);
					}
				}
			} else if (formatIdx == VCFRecord.FORMAT_IDX_PL) {
				//Phred likelihoods
				
				for(int j=0;j<alleles.length;j++) {
					for(int i=0;i<=j;i++) {
						if(i>0 || j>0) out.print (",");
						int condPhred = 0;
						if(report!=null && report.logConditionalsPresent()) {
							condPhred = (int) Math.round(-10*report.getLogConditionalProbability(alleles[i], alleles[j]));
						}
						out.print(condPhred);
					}
				}
			} else if (formatIdx == VCFRecord.FORMAT_IDX_GL) {
				//Likelihoods not phred scaled
				
				for(int j=0;j<alleles.length;j++) {
					for(int i=0;i<=j;i++) {
						if(i>0 || j>0) out.print (",");
						double logCond = 0;
						if(report!=null && report.logConditionalsPresent()) {
							logCond = report.getLogConditionalProbability(alleles[i], alleles[j]);
						}
						out.print(ParseUtils.ENGLISHFMT.format(logCond));
					}
				}
			} else if (formatIdx == VCFRecord.FORMAT_IDX_GQ) {
				//Phred of the genotype posterior
				out.print(var.getGenotypeQuality());
			} else if (formatIdx == VCFRecord.FORMAT_IDX_ACN) {
				//Local alleles copy number
				short totalCopyNumber = var.getCopyNumber();
				if(totalCopyNumber == 0) {
					out.print(VCFFileReader.NO_INFO_CHAR);
					continue;
				}
				short [] varAllelesCopyNumber = var.getAllelesCopyNumber();
				if(var.isUndecided()) varAllelesCopyNumber[0] = totalCopyNumber;
				for(int j=0;j<varAllelesCopyNumber.length;j++) {
					if(j>0) out.print (",");
					out.print(""+varAllelesCopyNumber[j]);
				}
			} else if (formatIdx == VCFRecord.FORMAT_IDX_DP) {
				//Read depth
				out.print(var.getTotalReadDepth());
			} else if (formatIdx == VCFRecord.FORMAT_IDX_ADP) {
				if(report!=null && report.countsPresent()) {
					for(int i=0;i<alleles.length;i++) {
						if(i>0) out.print(",");
						out.print(report.getCount(alleles[i]));
					}
				} else {
					for(int i=0;i<alleles.length;i++) {
						if(i>0) out.print(",");
						out.print("0");
					}
				}
			} else if (formatIdx == VCFRecord.FORMAT_IDX_BSDP) {
				if(var==null) {
					out.print(VCFFileReader.NO_INFO_CHAR);
					continue;
				}
				int [] allCounts = var.getAllCounts();
				if(allCounts == null) {
					out.print("0,0,0,0");
					continue;
				}
				for(int i=0;i<allCounts.length;i++) {
					if(i>0) out.print(",");
					out.print(allCounts[i]);
				}
			} else if (formatIdx == VCFRecord.FORMAT_IDX_RNC) {
				//Num copies
				if(cnv==null) {
					out.print(VCFFileReader.NO_INFO_CHAR);
					continue;
				}
				out.print(ParseUtils.ENGLISHFMT.format(cnv.getNumCopies()));
			} else if (formatIdx == VCFRecord.FORMAT_IDX_NTADF) {
				//Num tandem duplication fragments
				if(cnv==null) {
					out.print(VCFFileReader.NO_INFO_CHAR);
					continue;
				}
				out.print(cnv.getTandemFragments());
			} else if (formatIdx == VCFRecord.FORMAT_IDX_NTRDF) {
				//Num trans duplication fragments
				if(cnv==null) {
					out.print(VCFFileReader.NO_INFO_CHAR);
					continue;
				}
				out.print(cnv.getTransDupFragments());
			} else if (formatIdx == VCFRecord.FORMAT_IDX_TGEN) {
				//Text genotype
				if(cnv==null) {
					out.print(VCFFileReader.NO_INFO_CHAR);
					continue;
				}
				out.print(cnv.getTextGenotype());
			} else if (formatIdx == VCFRecord.FORMAT_IDX_NSF) {
				//TODO: NSF is more for SVs than for only CNVs
				if(cnv==null) {
					out.print(VCFFileReader.NO_INFO_CHAR);
					continue;
				}
				out.print(cnv.getTotalReadDepth());
			} 
		}
	}
	public void printHeader(VCFFileHeader header, PrintStream out) {
		header.print(out);
	}
	
}
