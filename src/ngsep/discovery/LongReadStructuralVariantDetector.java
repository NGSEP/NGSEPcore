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
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.CliquesFinder;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.CalledLargeIndel;
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
	
	private GenomicRegionSortedCollection<GenomicVariant> signatures;
	private ReferenceGenome refGenome;
	private int minMQ = DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	private GenomicRegionSortedCollection<GenomicVariant> variants;
	
	LongReadStructuralVariantDetector(){
	}
	public void setRefGenome(String referenceFile) throws IOException {
		this.refGenome = new ReferenceGenome(referenceFile);
	}
	public void readAlignments(String alignmentFile) throws IOException{
		try(ReadAlignmentFileReader alignmentReader = new ReadAlignmentFileReader(alignmentFile, refGenome)){
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			alignmentReader.setFilterFlags(filterFlags);
			alignmentReader.setMinMQ(minMQ);
			Map<String, List<GenomicVariant>> candidates = new LinkedHashMap<>();
			Iterator<ReadAlignment> it = alignmentReader.iterator();
			ReadAlignment aln = it.next();
			while(it.hasNext()) {
				aln = it.next();
				findSignatures(aln, candidates);
			}
			long startTime = System.currentTimeMillis();
			System.out.println(" Passing to sorted ds");
			QualifiedSequenceList seqs = new QualifiedSequenceList();
			for(String k:candidates.keySet()) seqs.add(new QualifiedSequence(k));
			signatures = new GenomicRegionSortedCollection<>(seqs);
			for(List<GenomicVariant> s:candidates.values()) signatures.addAll(s);
			signatures.forceSort();
			System.out.println(" sorted!");
			long endTime = System.currentTimeMillis();
			System.out.println("Total execution time of ds: " + (endTime-startTime) + "ms"); 
		}
	}
	public GenomicRegionSortedCollection<GenomicVariant> getSignatures() {
		return signatures;
	}
	public ReferenceGenome getRefGenome() {
		return refGenome;
	}
	public GenomicVariantImpl createSignature(GenomicVariant indel) {
		int first = indel.getFirst();
		int last = indel.getLast();
		String seqName = indel.getSequenceName(); 
		byte type = first + 1 == last ? GenomicVariant.TYPE_LARGEINS : GenomicVariant.TYPE_LARGEDEL;
		GenomicVariantImpl signature = new GenomicVariantImpl(seqName, first, last, type);
		if(signature.getType() == GenomicVariant.TYPE_LARGEINS) signature.setLast(first + 1); 
		return signature;
	}
	public void findSignatures(ReadAlignment aln, Map<String, List<GenomicVariant>> signatures) {
		Map<Integer, GenomicVariant> calls = aln.getIndelCalls();
		if(calls != null) {
			for (Map.Entry<Integer, GenomicVariant> call : calls.entrySet()) {
				GenomicVariant indel = call.getValue();
				if(indel.length() < 10) continue;
				GenomicVariant signature = createSignature(indel);
				String seqName = signature.getSequenceName();
				System.out.println(" Indel first: " + indel.getFirst() + " last: " + indel.getLast() + " chr: "
				+ indel.getSequenceName() + " with length: " + indel.length());
				signatures.computeIfAbsent(seqName, (k) -> new ArrayList<>()).add(signature);
			}
		}
	}	
	public List<GenomicVariantAnnotation> annotateStructuralVariant(GenomicVariant variant){
		List<GenomicVariantAnnotation> annotations = new ArrayList<>();
		int length = variant.getType() == GenomicVariant.TYPE_LARGEDEL ? variant.length()*(-1) : variant.length();
		GenomicVariantAnnotation svLengthAnnot = new GenomicVariantAnnotation(variant,
				GenomicVariantAnnotation.ATTRIBUTE_SVLEN, length);
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
	public void makeGenotypeCalls(String alignmentFile, List<GenomicVariant> variantsList) throws IOException{
		try(ReadAlignmentFileReader alignmentReader = new ReadAlignmentFileReader(alignmentFile, refGenome)){
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			alignmentReader.setFilterFlags(filterFlags);
			alignmentReader.setMinMQ(minMQ);
			GenomicRegionComparator cmpClassInstance = new GenomicRegionComparator(variants.getSequenceNames());
			Iterator<ReadAlignment> it = alignmentReader.iterator();
			ReadAlignment aln = it.next();
			int i = 0;
			GenomicVariant var = variantsList.get(i);
			int totalAlnVarInt = 0;
			int refAlnCalls = 0;
			int altAlnCalls = 0;
			while(it.hasNext()) {
				int cmp = cmpClassInstance.compare(var, aln);
				if(cmp < -1) {
					aln = it.next();
				}else if(cmp > 1) {
					if(totalAlnVarInt > 0) {
						CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(var, -1);
					}
					if(i < variantsList.size() - 1) {
						i++;
						var = variantsList.get(i);
						totalAlnVarInt = 0;
						refAlnCalls = 0;
						altAlnCalls = 0;
					}
					else break;
				}
				else {
					int begin = var.getFirst();
					if(aln.getIndelCall(begin) == null) refAlnCalls++;
					else altAlnCalls++;
					totalAlnVarInt++;
				}
			}
		}
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
		if(MAX_CLIQUE_FINDER_ALGORITHM.equals(algorithm)) {
			LongReadVariantDetectorAlgorithm caller = new MaxCliqueClusteringDetectionAlgorithm();
			caller.setSignatures(signatures);
			variants = caller.callVariants();
		}
		List<GenomicVariant> variantsList = variants.asList();
		//for(GenomicVariant v:variantsList) {
 			//System.out.println(v.getSequenceName() + " begin: " + v.getFirst() + " end: " + v.getLast() + " type: " + 
				//	GenomicVariantImpl.getVariantTypeName(v.getType())+ " length: " + v.length() + " QS: " + v.getVariantQS());
		//}
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
	
