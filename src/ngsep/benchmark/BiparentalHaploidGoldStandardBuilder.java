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

public class BiparentalHaploidGoldStandardBuilder {
	private Map<String,String> files = new TreeMap<>();
	private int filesP1;
	private int filesP2;
	private ReferenceGenome genome;
	private QualifiedSequenceList sequenceNames;
	private GenomicRegionSortedCollection<CalledGenomicVariant> allCalls;
	private int posPrint = 0;
	public static void main(String[] args) throws Exception {
		BiparentalHaploidGoldStandardBuilder instance = new BiparentalHaploidGoldStandardBuilder();
		instance.loadGenome (args[0]);
		instance.loadDescriptor(args[1]);
		try (PrintStream out = new PrintStream(args[2])) {
			instance.buildGoldStandardF1(out);
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
	public void buildGoldStandardF1(PrintStream out) throws IOException {
		for(String fileId:files.keySet()) {
			loadFileCalls(files.get(fileId),fileId);
		}
		VCFFileWriter writer = new VCFFileWriter();
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		header.addDefaultSample("GSBiparentalHaploid");
		writer.printHeader(header, out);
		for(QualifiedSequence sequence: sequenceNames) {
			String seqName = sequence.getName();
			List<CalledGenomicVariant> callsSequence = allCalls.getSequenceRegions(seqName).asList();
			int first = 0;
			int last = 0;
			List<CalledGenomicVariant> clusterCalls=new ArrayList<>();
			for(CalledGenomicVariant call:callsSequence) {
				if(call.getFirst()>last) {
					VCFRecord record = processCluster(seqName,first,last,clusterCalls, header);
					if(record!=null) writer.printVCFRecord(record, out);
					printReferenceCalls(seqName, last+1,call.getFirst()-1, header, writer,out);
					first = call.getFirst();
					clusterCalls.clear();
				}
				clusterCalls.add(call);
				last = Math.max(last, call.getLast());
			}
			VCFRecord record = processCluster(seqName,first,last,clusterCalls,header);
			if(record!=null) writer.printVCFRecord(record, out);
			printReferenceCalls(seqName,last+1,sequence.getLength(),header, writer,out);
			
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
	private void printReferenceCalls(String sequenceName, int first, int last, VCFFileHeader header, VCFFileWriter writer, PrintStream out) {
		if(last<first) return;
		CharSequence refAllele = genome.getReference(sequenceName,first,last);
		if(refAllele==null) return;
		if(refAllele.length()==0) throw new RuntimeException("Empty reference at "+sequenceName+":"+first+"-"+last);
		List<String> alleles = new ArrayList<>();
		alleles.add(refAllele.toString());
		GenomicVariant refVar = new GenomicVariantImpl(sequenceName, first, last, alleles);
		byte [] calledAlleles = {0};
		CalledGenomicVariantImpl call = new CalledGenomicVariantImpl(refVar, calledAlleles);
		writer.printVCFRecord(new VCFRecord(refVar, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, call, header), out);
	}
	private VCFRecord processCluster(String sequenceName, int first, int last, List<CalledGenomicVariant> clusterCalls, VCFFileHeader header) {
		if(last==0) return null;
		//Extract alleles
		CharSequence refAlleleC = genome.getReference(sequenceName,first,last);
		if(refAlleleC==null) return null;
		String refAllele = refAlleleC.toString();
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
		if(aP1<0 || aP2<0) {
			calledAlleles = new byte[0];
		} else if (aP1==aP2) {
			calledAlleles = new byte[1];
			calledAlleles[0] = (byte) (aP1==0?0:1);
		} else {
			calledAlleles = new byte[2];
			calledAlleles[0] = (byte) (aP1==0?0:1);
			calledAlleles[1] = (byte) (aP2==0?0:finalAlleles.size()-1);
			if(calledAlleles[0]> calledAlleles[1]) {
				byte tmp = calledAlleles[0];
				calledAlleles[0] = calledAlleles[1];
				calledAlleles[1] = tmp;
			}
		}
		
		GenomicVariantImpl outVariant = new GenomicVariantImpl(sequenceName, first, last, finalAlleles);
		CalledGenomicVariantImpl outCall = new CalledGenomicVariantImpl(outVariant, calledAlleles);
		int [] format = VCFRecord.DEF_FORMAT_ARRAY_MINIMAL;
		if(calledAlleles.length>0 && (calledAlleles.length>1 || calledAlleles[0]!=0)) {
			outVariant.setVariantQS((short) 255);
			byte type = GenomicVariant.TYPE_UNDETERMINED;
			for(CalledGenomicVariant call:clusterCalls) {
				type = GenomicVariantImpl.mergeType(type, call.getType());
			}
			if(type == GenomicVariant.TYPE_UNDETERMINED && first==last && finalAlleles.size()==2 && calledAlleles.length>0) type = GenomicVariant.TYPE_BIALLELIC_SNV;
			outVariant.setType(type);
			outCall.setGenotypeQuality((short)255);
			format = VCFRecord.DEF_FORMAT_ARRAY_QUALITY;
		}
		
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
				outCall.append(left);
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
			outCall.append(right);
		} else if (nextPosCall>last+1) {
			//Calls go over the last position
			return null;
		}
		if(first == posPrint) System.out.println("For "+sequenceName+":"+first+"-"+last+". File: "+fileId+" final call: "+outCall.toString());
		return outCall.toString();
	}
}
