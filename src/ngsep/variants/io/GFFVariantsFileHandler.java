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
package ngsep.variants.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import ngsep.main.io.ParseUtils;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledCNV;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.ReadPairCalledGenomicVariant;


public class GFFVariantsFileHandler {
	public static final String TAG_ID = "ID";
	public static final String TAG_GROUP_ID = "GROUP_ID";
	public static final String TAG_LENGTH = "LENGTH";
	public static final String TAG_SOURCE = "SOURCE";
	public static final String TAG_NSF = "NSF";
	public static final String TAG_NC = "NC";
	public static final String TAG_NUF = "NUF";
	public static final String TAG_HET = "HET";
	public static final String TAG_NTADF = "NTADF";
	public static final String TAG_NTRDF = "NTRDF";
	public static final String TAG_TGEN = "TGEN";
	public static final String TAG_NSR = "NSR";
	
	public List<CalledGenomicVariant> loadVariants(String filename) throws IOException {
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(filename);
			if(filename.endsWith(".gz")) {
				return loadVariants (new GZIPInputStream(fis));
			}
			return loadVariants(fis);
		} finally {
			if (fis!=null) fis.close();
		}
		
	}
	public List<CalledGenomicVariant> loadVariants(InputStream is) throws IOException {
		List<CalledGenomicVariant> answer = new ArrayList<CalledGenomicVariant>();
		
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(is));
			String line = in.readLine();
			QualifiedSequenceList seqNames = new QualifiedSequenceList();
			QualifiedSequenceList sources = new QualifiedSequenceList();
			
			while (line != null) {
				String[] items = line.split("\t");
				if(items.length<9) {
					line = in.readLine();
					continue;
				}
				String seqName = seqNames.addOrLookupName(items[0]).getName();
				int first = Integer.parseInt(items[3]);
				int last = Integer.parseInt(items[4]);
				byte type = GenomicVariantImpl.getVariantTypeId(items[2]);
				short genotypeQuality = Short.parseShort(items[5]);
				String[] items2 = items[8].split("=|;");
				GenomicVariantImpl var = new GenomicVariantImpl(seqName, first, last, type);
				if(type == GenomicVariant.TYPE_CNV || type == GenomicVariant.TYPE_REPEAT) {
					CalledCNV calledCNV = new CalledCNV(var);
					String textGen = null;
					for(int i=0;i<items2.length;i+=2) {
						String tag = items2[i];
						String value = items2[i+1];
						if(TAG_ID.equals(tag)) {
							var.setId(value);
						} else if (TAG_SOURCE.equals(tag)) {
							String source = sources.addOrLookupName(value).getName();
							calledCNV.setSource(source);
						} else if (TAG_NSF.equals(tag)) {
							if(type == GenomicVariant.TYPE_REPEAT) calledCNV.setNonUniqueAlns(Integer.parseInt(value));
							else calledCNV.setTotalReadDepth(Integer.parseInt(value));
						} else if (TAG_NC.equals(tag)) {
							calledCNV.setNumCopies(Float.parseFloat(value),true);
						} else if (TAG_NUF.equals(tag)) {
							//Backwards compatibility
							calledCNV.setUniqueAlns((int)Math.round(Double.parseDouble(value)));
						} else if (TAG_HET.equals(tag)) {
							calledCNV.setHeterozygousVariants(Integer.parseInt(value));
						} else if (TAG_NTADF.equals(tag)) {
							calledCNV.setTandemFragments(Integer.parseInt(value));
						} else if (TAG_NTRDF.equals(tag)) {
							calledCNV.setTransDupFragments(Integer.parseInt(value));
						} else if (TAG_TGEN.equals(tag)) {
							if(CalledCNV.TEXT_GEN_DEL.equals(value)) textGen = CalledCNV.TEXT_GEN_DEL;
							else if(CalledCNV.TEXT_GEN_TANDEMDUP.equals(value)) textGen = CalledCNV.TEXT_GEN_TANDEMDUP;
							else if(CalledCNV.TEXT_GEN_TRANSDUP.equals(value)) textGen = CalledCNV.TEXT_GEN_TRANSDUP;
							else System.err.println("Unrecognized text genotype "+value+" for CNV at "+var.getSequenceName()+":"+var.getFirst()+"-"+var.getLast());
							calledCNV.setTextGenotype(textGen);
						}
					}
					calledCNV.setGenotypeQuality(genotypeQuality);
					answer.add(calledCNV);
				} else if (type == GenomicVariant.TYPE_LARGEDEL || type == GenomicVariant.TYPE_LARGEINS || type == GenomicVariant.TYPE_INVERSION ) {
					int splitReads = 0;
					int predictedLength = 0;
					int numAlns = 0;
					for(int i=0;i<items2.length;i+=2) {
						String tag = items2[i];
						String value = items2[i+1];
						if(TAG_ID.equals(tag)) {
							var.setId(value);
						} else if (TAG_NSF.equals(tag)) {
							numAlns = Integer.parseInt(value);
						} else if (TAG_LENGTH.equals(tag)) {
							predictedLength = Integer.parseInt(value);
						} else if (TAG_NSR.equals(tag)) {
							splitReads = Integer.parseInt(value);
						}
					}
					ReadPairCalledGenomicVariant call = new ReadPairCalledGenomicVariant(var, CalledGenomicVariant.GENOTYPE_HOMOALT, predictedLength);
					call.setGenotypeQuality(genotypeQuality);
					call.setSupportingFragments(numAlns);
					call.setNumSplitReads(splitReads);
					answer.add(call);
				}
				line = in.readLine();
			}
		} finally {
			if(in!=null) in.close();
		}
		return answer;
	}
	public void saveVariants(List<? extends CalledGenomicVariant> variants, PrintStream out) {
		DecimalFormat df = ParseUtils.ENGLISHFMT;
		out.println("##gff-version	3");
		int index = 1;
		Set<String> usedIds = new TreeSet<String>();
		for(CalledGenomicVariant v:variants) if(v.getId()!=null)usedIds.add(v.getId());
		for(CalledGenomicVariant v:variants) {
			out.print(v.getSequenceName()+"\tNGSEP\t"+GenomicVariantImpl.getVariantTypeName(v.getType())+"\t"+v.getFirst()+"\t"+v.getLast());
			out.print("\t"+v.getGenotypeQuality()+"\t+\t.\t");
			String id = v.getId();
			while(id==null) {
				String idTest = "SV_"+index;
				if(!usedIds.contains(idTest)) id = idTest;
				index++;
			}
			out.print("ID="+id+";");
			
			out.print("LENGTH="+v.length());
			
			if(v instanceof CalledCNV) {
				CalledCNV cnv = (CalledCNV) v;
				if (cnv.getSource()!=null) out.print(";"+TAG_SOURCE+"="+cnv.getSource());
				if(v.getType() == GenomicVariant.TYPE_REPEAT) {
					out.print(";"+TAG_NSF+"="+cnv.getNonUniqueAlns());
				} else {
					out.print(";"+TAG_NSF+"="+cnv.getTotalReadDepth());
				}
				out.print(";"+TAG_NC+"="+df.format(cnv.getNumCopies()));
				if(cnv.getUniqueAlns()>0) out.print(";"+TAG_NUF+"="+cnv.getUniqueAlns());
				out.print(";"+TAG_HET+"="+cnv.getHeterozygousVariants());
				out.print(";"+TAG_NTADF+"="+cnv.getTandemFragments());
				out.print(";"+TAG_NTRDF+"="+cnv.getTransDupFragments());
				if(cnv.getTextGenotype()!=null) out.print(";"+TAG_TGEN+"="+cnv.getTextGenotype());
			} else if(v instanceof ReadPairCalledGenomicVariant) {
				ReadPairCalledGenomicVariant call = (ReadPairCalledGenomicVariant)v;
				if (call.getSource()!=null) out.print(";"+TAG_SOURCE+"="+call.getSource());
				out.print(";"+TAG_NSF+"="+call.getSupportingFragments());
				out.print(";"+TAG_NSR+"="+call.getNumSplitReads());
			}
			out.println();
			index++;
		}
	}
}
