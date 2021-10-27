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
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class LongReadStructuralVariantDetector {
	
	public static final int DEF_MIN_MQ_UNIQUE_ALIGNMENT = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	public static final double HETEROCIGOZITY_TRESHOLD = 0.8;
	public static final String MAX_CLIQUE_FINDER_ALGORITHM = "Clique";
	public static final String FILTER_SETTING_MISSING = ".";
	public static final String KEY_SEPARATOR = ",";
	
	private GenomicRegionSortedCollection<GenomicRegion> signatures;
	private GenomicRegionSortedCollection<GenomicVariant> variants;
	private ReferenceGenome refGenome;
	private int indelTresholdSize = 30;
	private int minMQ = DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
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
				if(indel.length() < indelTresholdSize) continue;
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
	public VCFFileHeader createVCFHeader(String sampleId) {
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		//header.addDefaultSample("C24");
		header.addSample(new Sample(sampleId), true);
		return header;
	}
	public List<String> getVariantFilters(GenomicVariant variants){
		List<String> filters = new ArrayList<>();
		filters.add(FILTER_SETTING_MISSING);
		return filters;
	}
	public List<CalledGenomicVariant> makeGenotypeCalls(String alignmentFile, List<GenomicVariant> variantsList) throws IOException{
		List<CalledGenomicVariant> genotypeCalls = new ArrayList<>();
		try(ReadAlignmentFileReader alignmentReader = new ReadAlignmentFileReader(alignmentFile, refGenome)){
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			alignmentReader.setFilterFlags(filterFlags);
			alignmentReader.setMinMQ(minMQ);
			GenomicRegionComparator cmpClassInstance = new GenomicRegionComparator(variants.getSequenceNames());
			Iterator<ReadAlignment> it = alignmentReader.iterator();
			ReadAlignment aln = it.next();	
			int varIdx = 0;
			int genotype = -1;
			int currIntersectingVars = 0;
			int [] alnCalls = new int [variantsList.size()];
			int [] alnCoveredPerVar = new int [variantsList.size()];
			boolean lastVar = false;
			while(it.hasNext()) {
				GenomicVariant var = variantsList.get(varIdx);
				int cmp = cmpClassInstance.compare(var, aln);
				if(cmp > 1 || lastVar) {
					varIdx -= currIntersectingVars;
					currIntersectingVars = 0;
					if(it.hasNext()) aln = it.next();
					lastVar = false;
				}
				else if(cmp < -1) {
					genotype = (double) alnCalls[varIdx] / alnCoveredPerVar[varIdx] > HETEROCIGOZITY_TRESHOLD ? 2 : 1;
					//System.out.println("het% " + (double) alnCalls[varIdx] / alnCoveredPerVar[varIdx]);
					CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(var, genotype);
					genotypeCalls.add(calledVariant);
					//System.out.println("genotype " + genotype);
					if(varIdx + 1 != variantsList.size()) varIdx++;
					else lastVar = true;
				}
				else {
					alnCoveredPerVar[varIdx]++;
					if(isSVAlleleAlnCall(aln, var)) alnCalls[varIdx]++;
					currIntersectingVars++;
					if(varIdx + 1 != variantsList.size()) varIdx++;
					else lastVar = true;
				}
			}
			//add remaining calls
			for(;varIdx < variantsList.size(); varIdx++) {
				GenomicVariant var = variantsList.get(varIdx);	
				genotype = (double) alnCalls[varIdx] / alnCoveredPerVar[varIdx] > HETEROCIGOZITY_TRESHOLD ? 2 : 1;
				CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(var, genotype);
				genotypeCalls.add(calledVariant);
				varIdx++;
			}
		}
		return genotypeCalls;
	}
	public boolean isSVAlleleAlnCall(ReadAlignment aln, GenomicVariant var) {
		GenomicVariant firstIndelCall = aln.getIndelCall(var.getFirst());
		byte firstIndelCallType = 0;
		if(firstIndelCall != null) {
			firstIndelCallType = firstIndelCall.getFirst() + 1 == firstIndelCall.getLast() ? 
					GenomicVariant.TYPE_LARGEINS : GenomicVariant.TYPE_LARGEDEL;
		}
		if(firstIndelCall == null || 
				firstIndelCallType != var.getType()) return false;
		return true;
	}
	public List<VCFRecord> buildRecords(List<GenomicVariant> variants, List<CalledGenomicVariant> genotypeCalls, VCFFileHeader header){
		List<VCFRecord> records = new ArrayList<>();
		for(int i = 0; i < variants.size(); i++) {
			GenomicVariant variant = variants.get(i);
			CalledGenomicVariant calledVariant = genotypeCalls.get(i);
			List<CalledGenomicVariant> calls = new ArrayList<>();
			calls.add(calledVariant);
			List<GenomicVariantAnnotation> infoFields = annotateStructuralVariant(variant);
			List<String> filters = getVariantFilters(variant);
			VCFRecord record = new VCFRecord(variant, filters,
					 infoFields, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL,  calls, header);
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
	public void run(String algorithm, String refFile, String alnFile, String sampleId) throws IOException {
		setRefGenome(refFile);
		readAlignments(alnFile);
		if(MAX_CLIQUE_FINDER_ALGORITHM.equals(algorithm)) {
			LongReadVariantDetectorAlgorithm caller = new MaxCliqueClusteringDetectionAlgorithm(refGenome, signatures);
			variants = caller.callVariants();
		}
		List<GenomicVariant> variantsList = variants.asList();
		List<CalledGenomicVariant> genotypeCalls = makeGenotypeCalls(alnFile, variantsList);
		VCFFileHeader header = createVCFHeader(sampleId);
		List<VCFRecord> records = buildRecords(variantsList, genotypeCalls, header);
		String saveFile = sampleId + ".variants.vcf";
		printVCFFile(records, header, saveFile);
	}
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		LongReadStructuralVariantDetector caller = new LongReadStructuralVariantDetector();
		caller.run(MAX_CLIQUE_FINDER_ALGORITHM, args[0], args[1], args[2]);
	}
}
	
