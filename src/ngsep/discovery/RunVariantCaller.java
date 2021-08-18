package ngsep.discovery;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import JSci.maths.statistics.SampleStatistics;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegion;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.CliquesFinder;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class RunVariantCaller {
	
	public static final int DEF_MIN_MQ_UNIQUE_ALIGNMENT = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	public static final String MAX_CLIQUE_FINDER_ALGORITHM = "Clique";
	
	private Map<String, List<Signature>> signatures;
	private ReferenceGenome refGenome;
	private int minMQ = DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	RunVariantCaller(){
		signatures = new LinkedHashMap<>();
	}
	
	public void setRefGenome(String referenceFile) throws IOException {
		refGenome = new ReferenceGenome(referenceFile);
	}
	
	public void readAlignments(String alignmentFile) throws IOException{
		try(ReadAlignmentFileReader alignmentReader = new ReadAlignmentFileReader(alignmentFile, refGenome)){
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			alignmentReader.setFilterFlags(filterFlags);
			alignmentReader.setMinMQ(minMQ);
			Iterator<ReadAlignment> it = alignmentReader.iterator();
			ReadAlignment aln = it.next();
			while(it.hasNext()) {
				aln = it.next();
				findSignatures(aln);
				//System.out.println(" Alignment processed at: " + aln.getSequenceName() + " pos: " + aln.getFirst());
			}
		}
	}
	
	public Map<String, List<Signature>> getSignatures() {
		return signatures;
	}
	
	public ReferenceGenome getRefGenome() {
		return refGenome;
	}
	/**
	public boolean isAlnGap(char cigOp) {
		return (cigOp == 'N' || cigOp == 'D') || cigOp == 'I';
	}
	**/
	public void findSignatures(ReadAlignment aln) {
		/**
		String cigarStr = aln.getCigarString();
		String strAmount = "";
		int amount = 0;
		int currentPos = aln.getFirst();
		for(int i = 0; i < cigarStr.length(); i++) {
			char op = cigarStr.charAt(i);
			if(Character.isDigit(op)) strAmount += op;
			else {
				amount = Integer.parseInt(strAmount);
				if(amount >= 50 && isAlnGap(op)) {
					int begin = currentPos;
					int end = currentPos + amount;
					int length = end - begin + 1;
					if(op == 'I') end = begin + 1;
					String chr = aln.getSequenceName();
					Signature sig = new Signature(chr, begin, end, length, Signature.ALN_TYPENAME_INTRALIGNMENT);
					sig.setType(op);
					if(signatures.containsKey(chr)) {
						List<Signature> signs = signatures.get(chr);
						signs.add(sig);
						signatures.put(chr, signs);
					}else {
						List<Signature> signs = new ArrayList<>();
						signs.add(sig);
						signatures.put(chr, signs);
					}
				}
				currentPos += amount;
				strAmount = "";
				amount = 0;
			}
		}
		**/
		Map<Integer, GenomicVariant> calls = aln.getIndelCalls();
		if(calls != null) {
			for (Map.Entry<Integer, GenomicVariant> call : calls.entrySet()) {
				GenomicVariant indel = call.getValue();
				if(indel.length()<50) continue;
				int begin = indel.getFirst();
				int end = indel.getLast();
				String seqName = indel.getSequenceName();
				int length = indel.length();
				byte varType = indel.getType();
				Signature sig = new Signature(seqName, begin, end, length, Signature.ALN_TYPENAME_INTRALIGNMENT);
				sig.setType(GenomicVariantImpl.getVariantTypeName(varType));
				System.out.println(" Indel first: " + indel.getFirst() + " last: " + indel.getLast() + " chr: " + indel.getSequenceName() + " with length: " + indel.length());
				if(signatures.containsKey(seqName)) {
					List<Signature> signs = signatures.get(seqName);
					signs.add(sig);
					signatures.put(seqName, signs);
				}else {
					List<Signature> signs = new ArrayList<>();
					signs.add(sig);
					signatures.put(seqName, signs);
				}
			}
		}
	}
	
	public List<GenomicVariantAnnotation> annotateStructuralVariant(GenomicVariant variant){
		List<GenomicVariantAnnotation> annotations = new ArrayList<>();
		GenomicVariantAnnotation svLengthAnnot = new GenomicVariantAnnotation(variant,
				GenomicVariantAnnotation.ATTRIBUTE_SVLEN, variant.length());
		GenomicVariantAnnotation svEndAnnot = new GenomicVariantAnnotation(variant,
				GenomicVariantAnnotation.ATTRIBUTE_END, variant.length());
		GenomicVariantAnnotation svTypeAnnot = new GenomicVariantAnnotation(variant,
				GenomicVariantAnnotation.ATTRIBUTE_SVTYPE, variant.length());
		annotations.add(svLengthAnnot);
		annotations.add(svEndAnnot);
		annotations.add(svTypeAnnot);
		return annotations;
	}

	public VCFFileHeader createVCFHeader() {
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		header.addDefaultSample("C24");
		return header;
	}
	
	public List<VCFRecord> getRecords(List<GenomicVariant> variants, VCFFileHeader header){
		List<VCFRecord> records = new ArrayList<>();
		for(GenomicVariant variant:variants) {
			CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(variant, -1);
			List<GenomicVariantAnnotation> annotations = annotateStructuralVariant(variant);
			VCFRecord record = new VCFRecord(variant, VCFRecord.DEF_FORMAT_ARRAY_NONE,
					calledVariant, header);
			records.add(record);
		}
		return records;
	}
	
	public void printVCFFile(List<VCFRecord> records,VCFFileHeader header, String file) throws IOException {
		PrintStream pr = new PrintStream(file);
		VCFFileWriter writer = new VCFFileWriter();
		writer.printHeader(header, pr);
		writer.printVCFRecords(records, pr);
	}
	
	public void run(String algorithm, String refFile, String alnFile) throws IOException {
		List<GenomicVariant> variants = new ArrayList<>();
		setRefGenome(refFile);
		readAlignments(alnFile);
		if(MAX_CLIQUE_FINDER_ALGORITHM.equals(algorithm)) {
			LongReadVariantCallerAlgorithm caller = new MaxCliqueClusteringDetectionAlgorithm();
			caller.setSignatures(signatures);
			caller.setReferenceGenome(refGenome);
			variants = caller.callVariants();
		}
		for(GenomicVariant v:variants) {
 			System.out.println(v.getSequenceName() + " begin: " + v.getFirst() + " end: " + v.getLast() + " type: " + 
					GenomicVariantImpl.getVariantTypeName(v.getType())+ " length: " + v.length() + " QS: " + v.getVariantQS());
		}
		VCFFileHeader header = createVCFHeader();
		List<VCFRecord> records = getRecords(variants, header);
		String saveFile = "variants.vcf";
		printVCFFile(records, header, saveFile);
	}


	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		RunVariantCaller caller = new RunVariantCaller();
		caller.run(MAX_CLIQUE_FINDER_ALGORITHM, args[0], args[1]);
		/**Map<String, List<Signature>> signatures = caller.getSignatures();
		List<String> keys = new ArrayList<String>(signatures.keySet());
		int nSigns = 0;
		for(String chr:keys) {
			for(Signature sign:signatures.get(chr)) {
				System.out.println(sign.toString());
				nSigns += 1;
			}
		}
		System.out.println(nSigns + " were detected");
		List<List<Integer>> clusters = run.findVariantClusters();
		for(List<Integer> cluster:clusters) {
			System.out.println("Cluster composed of signatures: ");
			for(int c:cluster) {
				System.out.println(c);
			}
		}
		MaxCliqueClusteringDetectionAlgorithm alg = new MaxCliqueClusteringDetectionAlgorithm();
		alg.setSignatures(signatures);
		List<GenomicVariant> variants = alg.callVariants();
		System.out.println("End of test");
		for(GenomicVariant v:variants) {
 			System.out.println(v.getSequenceName() + " begin: " + v.getFirst() + " end: " + v.getLast() + " type: " + 
					GenomicVariantImpl.getVariantTypeName(v.getType())+ " length: " + v.length() + " QS: " + v.getVariantQS());
		}**/
	}
}



class Signature implements GenomicRegion{
	
	public static final String TYPENAME_UNDETERMINED = "UND";
	public static final String TYPENAME_INDEL = "INDEL";
	public static final String TYPENAME_DELETION = "DEL";
	public static final String TYPENAME_INSERTION = "INS";
	public static final String TYPENAME_DUPLICATION = "DUP";
	public static final String TYPENAME_INVERSION = "INV";
	public static final String TYPENAME_TANDEM = "TAN";
	public static final String TYPENAME_TRANSLOCATION_BREAKPOINT = "BRK";
	public static final String ALN_TYPENAME_INTRALIGNMENT = "INTRA";
	public static final String ALN_TYPENAME_INTERALIGNMENT = "INTER";

	private String signatureType;
	private String sequenceName;
	private String type;
	private int first;
	private int last;
	private int length;
	
	public Signature(String sn, int begin, int end, int length, String sType) {
		this.first = begin;
		this.last = end;
		this.sequenceName = sn;
		this.signatureType = sType;
		this.length = length;
	}
	
	public void setType(char cigOp) {
		if(cigOp == 'N' || cigOp == 'D') {
			type = TYPENAME_DELETION;
		}else if(cigOp == 'I') {
			type = TYPENAME_INSERTION;
		}else type = TYPENAME_UNDETERMINED;
	}
	
	public void setType(String type) {
		this.type = type;
	}
	
	public String getSequenceName() {
		return sequenceName;
	}

	public String getType() {
		return type;
	}
	
	public String getSignatureType() {
		return signatureType;
	}
	
	public int getFirst() {
		return first;
	}

	public int getLast() {
		return last;
	}

	@Override
	public boolean isNegativeStrand() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isPositiveStrand() {
		// TODO Auto-generated method stub
		return true;
	}

	@Override
	public int length() { 
		return length;
	}
	
	public String toString() {
		return "" + sequenceName + " " + first + " " + last + " " + type + " " + signatureType;
	}
}
