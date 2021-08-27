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
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.CliquesFinder;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class LongReadStructuralVariantDetector {
	
	public static final int DEF_MIN_MQ_UNIQUE_ALIGNMENT = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	public static final String MAX_CLIQUE_FINDER_ALGORITHM = "Clique";
	public static final String FILTER_SETTING_MISSING = ".";
	
	private Map<String, List<GenomicVariant>> signatures;
	private ReferenceGenome refGenome;
	private GenomicRegionSortedCollection<ReadAlignment> alignments;
	private int minMQ = DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	LongReadStructuralVariantDetector(){
		signatures = new LinkedHashMap<>();
		alignments = new GenomicRegionSortedCollection<ReadAlignment>();
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
			//List<ReadAlignment> tempAlns = new ArrayList<>();
			while(it.hasNext()) {
				aln = it.next();
				findSignatures(aln);
				//tempAlns.add(aln);
				//System.out.println(" Alignment processed at: " + aln.getSequenceName() + " pos: " + aln.getFirst());
			}
			//QualifiedSequenceList seqs = new QualifiedSequenceList();
			//for(String key:signatures.keySet()) seqs.add(new QualifiedSequence(key));
			//alignments = new GenomicRegionSortedCollection<ReadAlignment>(seqs);
			//alignments.addAll(tempAlns);
		}
	}
	public Map<String, List<GenomicVariant>> getSignatures() {
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
					GenomicVariant sig = new GenomicVariant(chr, begin, end, length, GenomicVariant.ALN_TYPENAME_INTRALIGNMENT);
					sig.setType(op);
					if(GenomicVariants.containsKey(chr)) {
						List<GenomicVariant> signs = GenomicVariants.get(chr);
						signs.add(sig);
						GenomicVariants.put(chr, signs);
					}else {
						List<GenomicVariant> signs = new ArrayList<>();
						signs.add(sig);
						GenomicVariants.put(chr, signs);
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
				GenomicVariant sig = new GenomicVariantImpl(seqName, begin, end, varType);
				sig.setType(varType);
				System.out.println(" Indel first: " + indel.getFirst() + " last: " + indel.getLast() + " chr: " + indel.getSequenceName() + " with length: " + indel.length());
				if(signatures.containsKey(seqName)) {
					List<GenomicVariant> signs = signatures.get(seqName);
					signs.add(sig);
					signatures.put(seqName, signs);
				}
				else {
					List<GenomicVariant> signs = new ArrayList<>();
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
				GenomicVariantAnnotation.ATTRIBUTE_END, variant.getLast());
		GenomicVariantAnnotation svTypeAnnot = new GenomicVariantAnnotation(variant,
				GenomicVariantAnnotation.ATTRIBUTE_SVTYPE, GenomicVariantImpl.getVariantTypeName(variant.getType()));
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
	public List<String> getVariantFilters(GenomicVariant variants){
		List<String> filters = new ArrayList<>();
		filters.add(FILTER_SETTING_MISSING);
		return filters;
	}
	public List<VCFRecord> getRecords(List<GenomicVariant> variants, VCFFileHeader header){
		List<VCFRecord> records = new ArrayList<>();
		for(GenomicVariant variant:variants) {
			CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(variant, -1);
			List<CalledGenomicVariant> calls = new ArrayList<>();
			calls.add(calledVariant);
			List<GenomicVariantAnnotation> infoFields = annotateStructuralVariant(variant);
			List<String> filters = getVariantFilters(variant);
			VCFRecord record = new VCFRecord(variant, filters,
					 infoFields, VCFRecord.DEF_FORMAT_ARRAY_NONE,  calls, header);
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
		setRefGenome(refFile);
		readAlignments(alnFile);
		GenomicRegionSortedCollection<GenomicVariant> variants = new GenomicRegionSortedCollection<>();
		if(MAX_CLIQUE_FINDER_ALGORITHM.equals(algorithm)) {
			LongReadVariantDetectorAlgorithm caller = new MaxCliqueClusteringDetectionAlgorithm();
			caller.setSignatures(signatures);
			caller.setReferenceGenome(refGenome);
			variants = caller.callVariants();
		}
		List<GenomicVariant> variantsList = variants.asList();
		for(GenomicVariant v:variantsList) {
 			System.out.println(v.getSequenceName() + " begin: " + v.getFirst() + " end: " + v.getLast() + " type: " + 
					GenomicVariantImpl.getVariantTypeName(v.getType())+ " length: " + v.length() + " QS: " + v.getVariantQS());
		}
		VCFFileHeader header = createVCFHeader();
		List<VCFRecord> records = getRecords(variantsList, header);
		String saveFile = "variants.vcf";
		printVCFFile(records, header, saveFile);
	}
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		LongReadStructuralVariantDetector caller = new LongReadStructuralVariantDetector();
		caller.run(MAX_CLIQUE_FINDER_ALGORITHM, args[0], args[1]);
	}
}
	
