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
package ngsep.benchmark;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.math.CountsRankHelper;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;
/**
 * @author Jorge Duitama
 */
public class BiparentalHaploidGoldStandardBuilder {
	private Map<String,String> files = new TreeMap<>();
	private int filesP1;
	private int filesP2;
	private ReferenceGenome genome;
	private QualifiedSequenceList sequenceNames;
	private GenomicRegionSortedCollection<CalledGenomicVariant> allCalls;
	private int posPrint = -1;
	public static void main(String[] args) throws Exception {
		BiparentalHaploidGoldStandardBuilder instance = new BiparentalHaploidGoldStandardBuilder();
		instance.loadGenome (args[0]);
		instance.loadDescriptor(args[1]);
		try (PrintStream outVCF = new PrintStream(args[2]+"_GS.vcf");
			 PrintStream outRegions = new PrintStream(args[2]+"_GS_regions.txt"	)) {
			instance.buildGoldStandardF1(outVCF, outRegions);
		}
	}
	public void loadGenome(String filename) throws IOException {	
		genome = new ReferenceGenome(filename);
		sequenceNames = genome.getSequencesMetadata();
		allCalls = new GenomicRegionSortedCollection<>(sequenceNames);
	}
	public void loadDescriptor(String filename) throws IOException {
		try (FileReader reader = new FileReader(filename); 
			 BufferedReader in = new BufferedReader(reader)) {
			int i=0;
			int j=0;
			String line = in.readLine();
			while (line != null) {
				String [] items = line.split("\t");
				String fileId=null;
				if("1".equals(items[0])) {
					fileId = "P1_"+i+"_"+items[1];
					i++;
				}
				if("2".equals(items[0])) {
					fileId = "P2_"+j+"_"+items[1];
					j++;
				}
				if(fileId!=null) files.put(fileId, items[1]);
				line = in.readLine();
			}
			filesP1 = i;
			filesP2 = j;
		}
	}
	public void buildGoldStandardF1(PrintStream outVCF, PrintStream outRegions) throws IOException {
		for(String fileId:files.keySet()) {
			loadFileCalls(files.get(fileId),fileId);
		}
		VCFFileWriter writer = new VCFFileWriter();
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		header.addDefaultSample("GSBiparentalHaploid");
		writer.printHeader(header, outVCF);
		for(QualifiedSequence sequence: sequenceNames) {
			String seqName = sequence.getName();
			List<CalledGenomicVariant> callsSequence = allCalls.getSequenceRegions(seqName).asList();
			int first = 0;
			int last = 0;
			int nextConfidence = 1;
			List<CalledGenomicVariant> clusterCalls=new ArrayList<>();
			for(CalledGenomicVariant call:callsSequence) {
				if(call.getFirst()>last) {
					VCFRecord record = processCluster(seqName,first,last,clusterCalls, header);
					if(record!=null) {
						if(record.getFirst()==posPrint) System.out.println("Printing record for variant "+record.getSequenceName()+": "+record.getFirst()+"-"+record.getLast()+" type: "+record.getVariant().getType());
						writer.printVCFRecord(record, outVCF);
					} else {
						saveConfidenceRegion (seqName, nextConfidence, first-1, outRegions);
						nextConfidence = last+1;
					}
					first = call.getFirst();
					clusterCalls.clear();
				}
				clusterCalls.add(call);
				last = Math.max(last, call.getLast());
			}
			VCFRecord record = processCluster(seqName,first,last,clusterCalls,header);
			if(record!=null) {
				writer.printVCFRecord(record, outVCF);
			} else {
				saveConfidenceRegion (seqName, nextConfidence, first-1, outRegions);
				nextConfidence = last+1;
			}
			saveConfidenceRegion (seqName, nextConfidence, sequence.getLength(), outRegions);
		}
	}
	
	
	private void loadFileCalls(String filename, String sampleId) throws IOException {
		System.out.println("Loading calls from file: "+filename);
		try (VCFFileReader reader = new VCFFileReader(filename)) {
			Iterator<VCFRecord> it = reader.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				CalledGenomicVariant call = record.getCalls().get(0);
				call.setSampleId(sampleId);
				allCalls.add(call);
			}
		}
	}
	private void saveConfidenceRegion(String sequenceName, int first, int last, PrintStream out) {
		//Avoid very small reference regions
		if(first == posPrint) System.out.println("Seq name: "+sequenceName+" first: "+first+" last: "+last);
		if(last-first<5) return;
		CharSequence refAllele = genome.getReference(sequenceName,first,last);
		if(first == posPrint) System.out.println("Seq name: "+sequenceName+" first: "+first+" last: "+last+" ref allele: "+refAllele);
		if(refAllele==null) return;
		out.println(sequenceName+"\t"+first+"\t"+last);
	}
	private VCFRecord processCluster(String sequenceName, int first, int last, List<CalledGenomicVariant> clusterCalls, VCFFileHeader header) {
		if(last==0) return null;
		//Extract alleles
		CharSequence refAlleleC = genome.getReference(sequenceName,first,last);
		if(refAlleleC==null) return null;
		String refAllele = refAlleleC.toString().toUpperCase();
		if(refAllele.length()==0) throw new RuntimeException("Empty reference at "+sequenceName+":"+first+"-"+last+" clusterCalls: "+serializeCalls(clusterCalls));
		List<String> alleles = new ArrayList<>();
		alleles.add(refAllele);
		Set<String> fileIds = files.keySet();
		CountsRankHelper<Byte> countsP1 = new CountsRankHelper<>();
		CountsRankHelper<Byte> countsP2 = new CountsRankHelper<>();
		
		for(String fileId:fileIds) {
			String call = makeCall(sequenceName, first, last, refAllele, clusterCalls, fileId);
			if(call==null) continue;
			int j = alleles.indexOf(call);
			if(j<0) {
				j = alleles.size();
				if(first == posPrint) System.out.println("For "+sequenceName+":"+first+"-"+last+". Adding allele: "+call+" from file: "+fileId);
				alleles.add(call);
			}
			if(fileId.startsWith("P1")) {
				countsP1.add((byte)j);
			} else {
				countsP2.add((byte)j);
			}
		}
		byte aP1 = -1;
		if(countsP1.getNumDifferent()>0) {
			Map.Entry<Byte, Integer> bestP1 = countsP1.selectBest(1).entrySet().iterator().next();
			if(bestP1.getValue()>=filesP1-2) aP1 = bestP1.getKey(); 
			if(first == posPrint) System.out.println("For "+sequenceName+":"+first+"-"+last+". Best allele P1: "+bestP1.getKey()+" votes: "+bestP1.getValue()+" aP1: "+aP1);
		}
		byte aP2 = -1;
		if(countsP2.getNumDifferent()>0) {
			Map.Entry<Byte, Integer> bestP2 = countsP2.selectBest(1).entrySet().iterator().next();
			if(bestP2.getValue()>=filesP2-2) aP2 = bestP2.getKey(); 
			if(first == posPrint) System.out.println("For "+sequenceName+":"+first+"-"+last+". Best allele P2: "+bestP2.getKey()+" votes: "+bestP2.getValue()+" aP2: "+aP2);
		}
		List<String> finalAlleles = new ArrayList<>();
		finalAlleles.add(refAllele);
		
		if(aP1>0) finalAlleles.add(alleles.get(aP1));
		if(aP2>0 && aP2!=aP1) finalAlleles.add(alleles.get(aP2));
		byte [] calledAlleles;
		byte [] phasedAlleles=new byte[2];
		if(aP1<0 || aP2<0) {
			//Undecided call
			return null;
		} else if (aP1==aP2) {
			// Only variant calls are generated by this voting procedure
			if(aP1==0) return null;
			calledAlleles = new byte[1];
			calledAlleles[0] = 1;
			phasedAlleles[0] = phasedAlleles[1]= 1;
		} else {
			calledAlleles = new byte[2];
			phasedAlleles[0] = calledAlleles[0] = (byte) (aP1==0?0:1);
			phasedAlleles[1] = calledAlleles[1] = (byte) (aP2==0?0:finalAlleles.size()-1);
			if(calledAlleles[0]> calledAlleles[1]) {
				byte tmp = calledAlleles[0];
				calledAlleles[0] = calledAlleles[1];
				calledAlleles[1] = tmp;
			}
		}
		
		GenomicVariantImpl outVariant = new GenomicVariantImpl(sequenceName, first, last, finalAlleles);
		CalledGenomicVariantImpl outCall = new CalledGenomicVariantImpl(outVariant, calledAlleles);
		int [] format = VCFRecord.DEF_FORMAT_ARRAY_MINIMAL;
		outVariant.setVariantQS((short) 255);
		byte type = GenomicVariant.TYPE_UNDETERMINED;
		for(CalledGenomicVariant call:clusterCalls) {
			type = GenomicVariantImpl.mergeType(type, call.getType());
		}
		if(type == GenomicVariant.TYPE_UNDETERMINED && first==last && finalAlleles.size()==2) type = GenomicVariant.TYPE_BIALLELIC_SNV;
		outVariant.setType(type);
		outCall.setGenotypeQuality((short)255);
		if(phasedAlleles!=null) outCall.setIndexesPhasedAlleles(phasedAlleles);
		format = VCFRecord.DEF_FORMAT_ARRAY_QUALITY;
		
		//if(outCall.getFirst()==3518) System.out.println("Phased: "+outCall.isPhased()+" heterozygous: "+outCall.isHeterozygous()+" phased alleles: "+outCall.getPhasedAlleles()[0]+" "+outCall.getPhasedAlleles()[1]);
		VCFRecord answer = new VCFRecord(outVariant, format, outCall, header);
		return answer;
	}
	private String serializeCalls(List<CalledGenomicVariant> clusterCalls) {
		StringBuilder out = new StringBuilder();
		for(CalledGenomicVariant call:clusterCalls) {
			out.append(call.getSampleId()+":"+call.getAlleles().length+".");
		}
		return out.toString();
	}
	private String makeCall(String sequenceName, int first, int last, String refAllele, List<CalledGenomicVariant> clusterCalls, String fileId) {
		
		List<CalledGenomicVariant> callsFile = new ArrayList<>();
		for(CalledGenomicVariant call:clusterCalls) {
			if(fileId.equals(call.getSampleId())) {
				if(call.isUndecided() || call.isHeterozygous()) return null;
				if(call.isHomozygousReference()) continue;
				callsFile.add(call);
			}
		}
		//if(first == posPrint) System.out.println("For "+sequenceName+":"+first+"-"+last+". File: "+fileId+" clusterCalls: "+callsFile.size());
		if(callsFile.size()==0) return refAllele;
		if(first == last) return callsFile.get(0).getCalledAlleles()[0];
		StringBuilder outCall = new StringBuilder();
		int nextPosCall=first;
		for(CalledGenomicVariant call:callsFile) {
			//if(first == posPrint) System.out.println("For "+sequenceName+":"+first+"-"+last+". File: "+fileId+" call allele: "+call.getCalledAlleles()[0]+" quality: "+call.getGenotypeQuality());
			if(nextPosCall<call.getFirst()) {
				CharSequence left = genome.getReference(sequenceName, nextPosCall, call.getFirst()-1);
				if(left == null ) return null;
				outCall.append(left.toString().toUpperCase());
			} else if (nextPosCall>call.getFirst()) {
				//Overlapping calls within a single file
				return null;
			}
			outCall.append(call.getCalledAlleles()[0]);
			nextPosCall = call.getLast()+1;
		}
		if(nextPosCall<=last) {
			CharSequence right = genome.getReference(sequenceName, nextPosCall, last);
			if(right == null ) return null;
			outCall.append(right.toString().toUpperCase());
		} else if (nextPosCall>last+1) {
			//Calls go over the last position
			return null;
		}
		if(first == posPrint) System.out.println("For "+sequenceName+":"+first+"-"+last+". File: "+fileId+" final call: "+outCall.toString());
		return outCall.toString();
	}
}
