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
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
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
import ngsep.variants.CalledInversion;
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
	public static final String KEY_SEPARATOR = ",";
	
	private GenomicRegionSortedCollection<GenomicRegion> signatures;
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
			Map<String, Integer> candidates = new LinkedHashMap<>();
			Set<String> sequenceNames = new LinkedHashSet<>();
			Iterator<ReadAlignment> it = alignmentReader.iterator();
			ReadAlignment aln = it.next();
			while(it.hasNext()) {
				aln = it.next();
				findSignatures(aln, candidates, sequenceNames);
			}
			QualifiedSequenceList qualifiedSequenceNames = new QualifiedSequenceList();
			for(String seqName:sequenceNames) qualifiedSequenceNames.add(new QualifiedSequence(seqName));
			signatures = new GenomicRegionSortedCollection<>(qualifiedSequenceNames);
			for(Map.Entry<String, Integer> entry : candidates.entrySet()) {
				decodeKeyToSignature(entry.getKey(), entry.getValue());
			}
			signatures.forceSort();
		}
	}
	public GenomicRegionSortedCollection<GenomicRegion> getSignatures() {
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
	public void decodeKeyToSignature(String info, int count) {
		String [] signInfo = info.split(KEY_SEPARATOR);
		String seqName = signInfo[0];
		int first =  Integer.parseInt(signInfo[1]);
		int last =  Integer.parseInt(signInfo[2]);
		int length =  Integer.parseInt(signInfo[3]);
		byte type = first + 1 == last ? GenomicVariant.TYPE_LARGEINS : GenomicVariant.TYPE_LARGEDEL;
		if(type == GenomicVariant.TYPE_LARGEDEL || type == GenomicVariant.TYPE_LARGEINS) {
			CalledLargeIndel sign = new CalledLargeIndel(seqName, first, last, 
					type == GenomicVariant.TYPE_LARGEDEL, length, count);
			signatures.add(sign);
		}
		else if(type == GenomicVariant.TYPE_INVERSION) {
			CalledInversion sign = new CalledInversion(seqName, first, last, count);
			signatures.add(sign);
		}
	}
	public void computeIndelToSignatureKey(GenomicVariant indel, Map<String, Integer> signatures) {
		String signatureInfo = "";
		String seqName = indel.getSequenceName();
		String first = Integer.toString(indel.getFirst());
		String last = Integer.toString(indel.getLast());
		String length = Integer.toString(indel.length());
		signatureInfo = seqName + KEY_SEPARATOR + first + KEY_SEPARATOR + last +
				KEY_SEPARATOR + length;
		signatures.compute(signatureInfo, (k,v)-> (v == null) ? 1 : v + 1);
	}
	public void findSignatures(ReadAlignment aln, Map<String, Integer> signatures, Set<String> sequenceNames) {
		Map<Integer, GenomicVariant> calls = aln.getIndelCalls();
		if(calls != null) {
			for (Map.Entry<Integer, GenomicVariant> call : calls.entrySet()) {
				GenomicVariant indel = call.getValue();
				if(indel.length() < 10) continue;
				String seqName = indel.getSequenceName();
				computeIndelToSignatureKey(indel, signatures);
				sequenceNames.add(seqName);
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
	
