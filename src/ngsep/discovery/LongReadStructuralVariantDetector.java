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
import ngsep.genome.GenomicRegionImpl;
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
	public static final int DUP_DETERMINING_MAX_DISTANCE = 1000;
	public static final int DEL_INTER_DETERMINING_MAX_DISTANCE = 15000;
	public static final String MAX_CLIQUE_FINDER_ALGORITHM = "Clique";
	public static final String FILTER_SETTING_MISSING = ".";
	public static final String KEY_SEPARATOR = ",";
	
	private GenomicRegionSortedCollection<GenomicVariant> signatures;
	private GenomicRegionSortedCollection<GenomicVariant> variants;
	private Map<String, List<GenomicRegion>> alignmentRegions;
	
	private ReferenceGenome refGenome;
	private int indelTresholdSize;
	private int minMQ = DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	LongReadStructuralVariantDetector(){
	}
	public void setRefGenome(String referenceFile) throws IOException {
		this.refGenome = new ReferenceGenome(referenceFile);
	}
	private void readAlignments(String alignmentFile,  int SVLength) throws IOException{
		try(ReadAlignmentFileReader alignmentReader = new ReadAlignmentFileReader(alignmentFile, refGenome)){
			setIndelTresholdSize(SVLength);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			alignmentReader.setFilterFlags(filterFlags);
			alignmentReader.setMinMQ(minMQ);
			QualifiedSequenceList qualifiedSequenceNames = refGenome.getSequencesList();
			signatures = new GenomicRegionSortedCollection<>(qualifiedSequenceNames);
			alignmentRegions = new HashMap<>();
			Iterator<ReadAlignment> it = alignmentReader.iterator();
			ReadAlignment aln = it.next();
			while(it.hasNext()) {
				aln = it.next();
				computeAlnAsRegion(aln);
				findIntraAlnSignatures(aln);
				System.out.println("Only intra: " + signatures.size());
			}
			for(List <GenomicRegion> alnRegions : alignmentRegions.values()) {
				if(alnRegions.size() > 1) findInterAlnSignatures(alnRegions);
			}
			System.out.println("All: " + signatures.size());
			signatures.forceSort();
		}
	}
	private void computeAlnAsRegion(ReadAlignment aln) {
		// TODO Auto-generated method stub
		GenomicRegion alnRegion = new GenomicRegionImpl(aln.getSequenceName(), aln.getFirst(), aln.getLast());
		List<GenomicRegion> alnRegions = alignmentRegions.computeIfAbsent(aln.getReadName(), v -> new ArrayList<>());
		alnRegions.add(alnRegion);
	}
	public GenomicRegionSortedCollection<GenomicVariant> getSignatures() {
		return signatures;
	}
	public ReferenceGenome getRefGenome() {
		return refGenome;
	}
	private GenomicVariant createIndelInterAlnSignature(String seqName, int first, int last, int length,
			byte type) {
		// TODO Auto-generated method stub
		GenomicVariantImpl signature = new GenomicVariantImpl(seqName, first, last, type);
		signature.setLength(length);
		return signature;
	}
	private GenomicVariantImpl createIndelIntraAlnSignature(GenomicVariant indel) {
		int first = indel.getFirst();
		int last = indel.getLast();
		int length = indel.length();
		String seqName = indel.getSequenceName(); 
		byte type = first + 1 == last ? GenomicVariant.TYPE_LARGEINS : GenomicVariant.TYPE_LARGEDEL;
		GenomicVariantImpl signature = new GenomicVariantImpl(seqName, first, last, type);
		signature.setLength(length);
		return signature;
	}
	private void findIntraAlnSignatures(ReadAlignment aln) {
		Map<Integer, GenomicVariant> calls = aln.getIndelCalls();
		if(calls != null) {
			for (Map.Entry<Integer, GenomicVariant> call : calls.entrySet()) {
				GenomicVariant indel = call.getValue();
				if(indel.length() < indelTresholdSize) continue;
				signatures.add(createIndelIntraAlnSignature(indel));
			}
		}
	}	
	private void findInterAlnSignatures(List<GenomicRegion> alnRegions) {
		// TODO Auto-generated method stub
		//System.out.println("Read aln list: ");
		int nAln = alnRegions.size();
		GenomicRegion firstAln = alnRegions.get(0);
		GenomicRegion lastAln = alnRegions.get(nAln - 1);
		String sequenceName = firstAln.getSequenceName();
		if(nAln == 2 ) {
			// && Math.abs(firstAln.getLast() - lastAln.getFirst())
			//Two aln that are very far away could be a DUP
			if(sequenceName.equals(lastAln.getSequenceName())) {
				int distance = Math.abs(firstAln.getLast() - lastAln.getFirst());
				if(distance <= DEL_INTER_DETERMINING_MAX_DISTANCE) {
					int first = firstAln.getLast() + 1;
					int last = lastAln.getFirst() - 1;
					int length = last - first + 1;
					byte type = GenomicVariant.TYPE_LARGEDEL;
					if(length >= indelTresholdSize) signatures.add(createIndelInterAlnSignature(sequenceName, first, last, length, type));
				}
			}
		}
		if(nAln == 3) {
			GenomicRegion midAln = alnRegions.get(1);
			if(sequenceName.equals(midAln.getSequenceName()) && sequenceName.equals(lastAln.getSequenceName())) {
				//A DEL in two aln could be an INS far away
				int firstDistance = Math.abs(midAln.getFirst() - firstAln.getLast());
				int secondDistance = Math.abs(midAln.getLast() - lastAln.getFirst());
				if(secondDistance > DUP_DETERMINING_MAX_DISTANCE && firstDistance > DUP_DETERMINING_MAX_DISTANCE) {
					int first = midAln.getFirst();
					int last = midAln.getLast();
					int length = last - first + 1;
					byte type = GenomicVariant.TYPE_INVERSION;
					if(length >= indelTresholdSize) signatures.add(createIndelInterAlnSignature(sequenceName, first, last, length, type));
				}
			}
		}
		/**
		for(int i = 0; i < alnRegions.size(); i++) {
			GenomicRegion current = alnRegions.get(i);
			System.out.println(current.getSequenceName() + "\t" + current.getFirst() + "\t" + current.getLast());
		}**/
	}
	private List<GenomicVariantAnnotation> annotateStructuralVariant(GenomicVariant variant){
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
	private VCFFileHeader createVCFHeader(String sampleId) throws IOException {
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		//header.addDefaultSample("C24");
		header.addSample(new Sample(sampleId), true);
		String end = "##INFO=<ID=END,Number=1,Type=Integer,Description=" + "\"End position of the structural variant\"" + ">";
		String svtype = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description="
				+ "\"Type of SV:DEL=Deletion, INS=Insertion, DUP=Duplication, INV=Inversion\"" +">";
		String svlen = "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=" + "\"Difference in length between REF and ALT alleles\"" + ">";
		header.loadHeaderLine(end);
		header.loadHeaderLine(svtype);
		header.loadHeaderLine(svlen);
		return header;
	}
	private List<String> getVariantFilters(GenomicVariant variants){
		List<String> filters = new ArrayList<>();
		filters.add(FILTER_SETTING_MISSING);
		return filters;
	}
	private List<CalledGenomicVariant> makeGenotypeCalls(String alignmentFile, List<GenomicVariant> variantsList) throws IOException{
		List<CalledGenomicVariant> genotypeCalls = new ArrayList<>();
		try(ReadAlignmentFileReader alignmentReader = new ReadAlignmentFileReader(alignmentFile, refGenome)){
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			alignmentReader.setFilterFlags(filterFlags);
			alignmentReader.setMinMQ(minMQ);
			GenomicRegionComparator cmpClassInstance = new GenomicRegionComparator(variants.getSequenceNames());
			Iterator<ReadAlignment> it = alignmentReader.iterator();
			ReadAlignment aln = it.next();	
			int varIdx = 0;
			int genotype;
			int currIntersectingVars = 0;
			int [] alnCalls = new int [variantsList.size()];
			int [] alnCoveredPerVar = new int [variantsList.size()];
			boolean lastVar = false;
			while(it.hasNext()) {
				GenomicVariant var = variantsList.get(varIdx);
				//System.out.println("var=" + var.getSequenceName() + " " + var.getFirst() 
				//+ " " + var.getLast() + " "  + var.length());
				//System.out.println("aln=" + aln.getSequenceName() + " " + aln.getFirst() 
				//+ " " + aln.getLast() + " "  + aln.length());
				int cmp;
				//if(aln.getLast() == aln.getFirst()) cmp = modifiedAlnComparison(var, aln, cmpClassInstance);
				//else cmp = cmpClassInstance.compare(var, aln);
				cmp = cmpClassInstance.compare(var, aln);
				//System.out.println("cmp=" + cmp);
				if(cmp > 1 || lastVar) {
					varIdx -= currIntersectingVars;
					currIntersectingVars = 0;
					if(it.hasNext()) aln = it.next();
					lastVar = false;
				}
				else if(cmp < -1) {
					genotype = decideGenotype(alnCalls[varIdx], alnCoveredPerVar[varIdx]);
					//System.out.println("calls% " + alnCalls[varIdx]);
					//System.out.println("all% " + alnCoveredPerVar[varIdx]);
					//System.out.println("het% " + (double) alnCalls[varIdx] / alnCoveredPerVar[varIdx]);
					CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(var, genotype);
					genotypeCalls.add(calledVariant);
					//System.out.println("genotype " + genotype);
					if(varIdx + 1 != variantsList.size()) varIdx++;
					else lastVar = true;
				}
				else {
					alnCoveredPerVar[varIdx]++;
					//Implement something for interaln calls
					if(hasSVAlleleAlnCall(aln, var)) alnCalls[varIdx]++;
					currIntersectingVars++;
					if(varIdx + 1 != variantsList.size()) varIdx++;
					else lastVar = true;
				}
			}
			//add remaining calls
			for(;varIdx < variantsList.size(); varIdx++) {
				GenomicVariant var = variantsList.get(varIdx);	
				genotype = decideGenotype(alnCalls[varIdx], alnCoveredPerVar[varIdx]);
				CalledGenomicVariant calledVariant = new CalledGenomicVariantImpl(var, genotype);
				genotypeCalls.add(calledVariant);
				varIdx++;
			}
		}
		return genotypeCalls;
	}
	private boolean hasSVAlleleAlnCall(ReadAlignment aln, GenomicVariant var) {
		int varFirst = var.getFirst();
		int varLength = var.length();
		int baseQueueMidWindowSize = (int) (varLength * 0.10);
		int baseQueueBegin = varFirst - baseQueueMidWindowSize;
		int baseQueueLimit = varFirst + baseQueueMidWindowSize;
		int lengthDiffLimit = (int) (varLength * 0.20);
		for(int i = baseQueueBegin; i <= baseQueueLimit; i++) {
			GenomicVariant candidateIndelCall = aln.getIndelCall(i);
			byte firstIndelCallType;
			if(candidateIndelCall != null) {
				firstIndelCallType = candidateIndelCall.getFirst() + 1 == candidateIndelCall.getLast() ? 
						GenomicVariant.TYPE_LARGEINS : GenomicVariant.TYPE_LARGEDEL;
				//System.out.println("call=" + candidateIndelCall.getSequenceName() + " " + candidateIndelCall.getFirst() 
					//	+ " " + candidateIndelCall.getLast() + " "  + candidateIndelCall.length());
				if(firstIndelCallType == var.getType()) {
					int callLength = candidateIndelCall.length();
					int lengthDiff = Math.abs(callLength - varLength);
					if(lengthDiff <= lengthDiffLimit) return true;
				}
			}
		}
		/**
		if(aln.getIndelCalls() != null) {
			for(Map.Entry<Integer, GenomicVariant> entry : aln.getIndelCalls().entrySet()) {
				GenomicVariant tv = entry.getValue();
				System.out.println("&sign=" + tv.getSequenceName() + " " + tv.getFirst() 
				+ " " + tv.getLast() + " "  + tv.length());
			}
		}**/
		return false;
	}
	private void filterRecordsByGenotype(List<VCFRecord> records) {
		// TODO Auto-generated method stub
		List<VCFRecord> toRemove = new ArrayList<>();
		for(int i = 0; i < records.size(); i++) {
			VCFRecord record = records.get(i);
			CalledGenomicVariant call = record.getCalls().get(0);
			if(call.isUndecided()) toRemove.add(record); 
		}
		records.removeAll(toRemove);
	}
	private int decideGenotype(int calls, int coverage) {
		double filterCallTreshold = 0.025;
		double percentage = (double) calls / coverage;
		return percentage > HETEROCIGOZITY_TRESHOLD ? 2 : (percentage > filterCallTreshold ? 1 : -1);
	}
	private List<VCFRecord> buildRecords(List<GenomicVariant> variants, List<CalledGenomicVariant> genotypeCalls, VCFFileHeader header){
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
	private void printVCFFile(List<VCFRecord> records,VCFFileHeader header, String file) throws IOException {
		PrintStream pr = new PrintStream(file);
		VCFFileWriter writer = new VCFFileWriter();
		writer.printHeader(header, pr);
		writer.printVCFRecords(records, pr);
	}
	public void run(String algorithm, String refFile, String alnFile, String sampleId, int SVLength) throws IOException {
		setRefGenome(refFile);
		readAlignments(alnFile, SVLength);
		if(MAX_CLIQUE_FINDER_ALGORITHM.equals(algorithm)) {
			LongReadVariantDetectorAlgorithm caller = new MaxCliqueClusteringDetectionAlgorithm(refGenome, signatures, SVLength);
			variants = caller.callVariants();
		}
		List<GenomicVariant> variantsList = variants.asList();
		List<CalledGenomicVariant> genotypeCalls = makeGenotypeCalls(alnFile, variantsList);
		VCFFileHeader header = createVCFHeader(sampleId);
		List<VCFRecord> records = buildRecords(variantsList, genotypeCalls, header);
		filterRecordsByGenotype(records);
		String saveFile = sampleId + ".variants.vcf";
		printVCFFile(records, header, saveFile);
	}
	private void setIndelTresholdSize(int SVLength) {
		this.indelTresholdSize = SVLength;
	}
	/**
	public int modifiedAlnComparison(GenomicVariant var, ReadAlignment aln, GenomicRegionComparator cmpClassInstance) {
		String seqName = aln.getSequenceName();
		int first =  aln.getFirst();
		int length = aln.length();
		int last =  aln.getFirst() + length;
		ReadAlignment modAln = new ReadAlignment(seqName, first, last, length, aln.getFlags());
		return cmpClassInstance.compare(var, modAln);
	}*/
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		LongReadStructuralVariantDetector caller = new LongReadStructuralVariantDetector();
		caller.run(MAX_CLIQUE_FINDER_ALGORITHM, args[0], args[1], args[2], Integer.parseInt(args[3]));
	}
}
	
