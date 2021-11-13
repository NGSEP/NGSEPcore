package ngsep.benchmark;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.io.FastaFileReader;
import ngsep.variants.AlleleCompatibilityGenomicVariantComparator;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class VCFDirectSVComparison {
	/**
	 * Finds the intersecting structural variants between two different VCF files
	 * @author Nicolas Gaitan
	 */
	private static final double PERCENTAGE_INTERSECTION_TRESHOLD = 0.3;
	private static final String CREATE_NOT_CALLED_VARIANTS_LENGTH_DISTRIBUTION_COMMAND = "Distribution";
	
	private ReferenceGenome genome;
	private List<GenomicVariant> referenceVariants;
	private List<GenomicVariant> testVariants;
	private List<GenomicVariant> notIntersectingReferenceVariants;
	private List<GenomicVariant> notIntersectingTestVariants;
	private List<GenomicVariant> intersectingTestVariants;
	private List<GenomicVariant> intersectingRefVariants;
	private List<Integer> truePositivesPerVar;
	private List<Integer> falsePositivesPerVar;
	private int lengthToDefineSVEvent = 30;
	private int treshold = 5000;
	private int testVariantCount;
	private int refVariantCount;
	private int truePositiveBps;
	private int falsePositiveBps;
	private int truePositiveBpsTreshold;
	private int falsePositiveBpsTreshold;
	
	public VCFDirectSVComparison() {
		referenceVariants = new ArrayList<>();
		testVariants = new ArrayList<>();
		truePositivesPerVar = new ArrayList<>();
		falsePositivesPerVar = new ArrayList<>();
		notIntersectingReferenceVariants = new ArrayList<>();
		notIntersectingTestVariants = new ArrayList<>();
		intersectingRefVariants = new ArrayList<>();
		intersectingTestVariants = new ArrayList<>();
	}
	public int findIntersectingVariants(String vcf1, String vcf2) throws IOException{
		try(VCFFileReader input1 = new VCFFileReader(vcf1);
				VCFFileReader input2 = new VCFFileReader(vcf2);){
			AlleleCompatibilityGenomicVariantComparator comparator = new AlleleCompatibilityGenomicVariantComparator(genome);
			input1.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			input2.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			Iterator<VCFRecord> it1 = input1.iterator();
			Iterator<VCFRecord> it2 = input2.iterator();
			int numIntersectingVariants = 0;
			VCFRecord testRecords = it1.next();
			VCFRecord refRecords = it2.next();
			GenomicVariantImpl testVariant = (GenomicVariantImpl) testRecords.getVariant();
			GenomicVariantImpl refVariant = (GenomicVariantImpl) refRecords.getVariant();
			int firstTest = -1;
			int firstRef = -1;
			int lastTest = -1;
			int lastRef = -1;
			testVariantCount = 1;
			while(testRecords != null && refRecords != null) {
				testVariant = (GenomicVariantImpl) testRecords.getVariant();
				refVariant = (GenomicVariantImpl) refRecords.getVariant();
				firstTest = testVariant.getFirst();
				firstRef = refVariant.getFirst();
				lastTest = testVariant.getLast();
				lastRef = refVariant.getLast();
				if(testVariant.getType() == GenomicVariant.TYPE_INDEL) {
					if(lastTest-firstTest < 2) {
						testVariant.setType(GenomicVariant.TYPE_LARGEINS);
					}else {
						testVariant.setType(GenomicVariant.TYPE_LARGEDEL);
					}
				}
				if(lastRef-firstRef < 2) {
					lastRef = firstRef + refVariant.length() - 1;
					refVariant.setLast(lastRef);
				}
				if(lastTest-firstTest < 2) {
					lastTest = firstTest + testVariant.length() - 1;
					testVariant.setLast(lastTest);
				}
				int cmp = comparator.compare(testVariant, refVariant);
				//System.out.println("Comparison: "+cmp + " for: " + 
					//	"firstTest: "+ firstTest + " lastTest: " + lastTest +
						//" firstRef: " + firstRef + " lastRef: " + lastRef);
				if(cmp > 1) {
					if(it2.hasNext()) {
						refRecords = it2.next();
						processRefVariantForOtherCalculations(refVariant);	
					}
					else break;
				}
				else if(cmp < -1) {
					if(it1.hasNext()) {
						testRecords = it1.next();
						testVariantCount++;
						falsePositiveBps += Math.abs(testVariant.length());
						if(testVariant.length() < treshold) falsePositiveBpsTreshold+= Math.abs(testVariant.length());
						testVariants.add(testVariant);
					}
					else break;
				}
				else {
					System.out.println(" ------------------------------------------------------------------------- ");
					if(testVariant.getType() == refVariant.getType()) {
						VariantIntersectionMetricsCalculator calc = 
								new VariantIntersectionMetricsCalculator(testVariant, refVariant);
						boolean isIntersection = isIntersection(calc.estimateIntersectionPercentageOverReference());
						updateStatisticsForBaseAnalysis(calc.estimateTruePositiveNucleotides(),
								calc.estimateFalsePositiveNucleotides(), calc.estimateTruePositiveNucleotides(treshold), 
								calc.estimateFalsePositiveNucleotides(treshold));
						if(isIntersection) {
							if(isSVByLength(refVariant.length())) {
								numIntersectingVariants++;
								intersectingRefVariants.add(refVariant);
								intersectingTestVariants.add(testVariant);
								
							}
						}
					}else {
						falsePositiveBps += Math.abs(testVariant.length());
						if(testVariant.length() < treshold) falsePositiveBpsTreshold+= Math.abs(testVariant.length());
						
					}
					if(it1.hasNext()) {
						testRecords = it1.next();
						testVariantCount++;
						testVariants.add(testVariant);
					}
					else break;
					if(it2.hasNext()) {
						refRecords = it2.next();
						processRefVariantForOtherCalculations(refVariant);
					}
					else break;
				}
			}
			processRefVariantForOtherCalculations(refVariant);
			estimatePositiveBps();
			testVariants.add(testVariant);
			return numIntersectingVariants;
		}
	}
	public void printReport(int ints) throws IOException {
		try(PrintWriter writer = new PrintWriter("statistics_report.txt")){
			writer.println(" ------------------------------------------------------------------------- ");
			writer.println("#Number of test variants: " + testVariantCount);
			writer.println("#Number of reference variants: " + refVariantCount);
			writer.println("#Intersections: " + ints);
			writer.println("#Recall (sensitivity) over intersections: " +  ((double) ints / refVariantCount)*100);
			double precision = (double) ints / (testVariantCount);
			writer.println("#Precision over intersections: " + precision*100);
			writer.println("#FDR over intersections: " + ((1 - precision)*100));
			writer.println("#Total True positive bases: " + truePositiveBps);
			writer.println("#Total False positive bases: " + falsePositiveBps);
			int positiveBps = estimatePositiveBps();
			writer.println("#Total positive bases: " + positiveBps);
			int positivesUnderTreshold = estimatePositiveBps(treshold);
			writer.println("#Total positive bases under " + treshold + " : " + positivesUnderTreshold);
			writer.println("#Recall (sensitivity) over bases: " +  ((double) truePositiveBps / positiveBps)*100);
			writer.println("#Recall (sensitivity) under " + treshold + " : " +
					((double) truePositiveBpsTreshold / positivesUnderTreshold)*100);
			double perBasePrecision = (double) truePositiveBps / (truePositiveBps + falsePositiveBps);
			writer.println("#Precision over bases: " + perBasePrecision*100);
			writer.println("#FDR over bases: " + ((1 - perBasePrecision)*100));
			writer.println(" ------------------------------------------------------------------------- ");
		}
	}
	public static void main(String[] args) throws IOException {
		VCFDirectSVComparison v = new VCFDirectSVComparison();
		v.setGenome(args[0]);
		String testFile = args[1];
		String refFile = args[2];
		int ints = v.findIntersectingVariants(testFile, refFile);
		if(args.length > 3 && CREATE_NOT_CALLED_VARIANTS_LENGTH_DISTRIBUTION_COMMAND.equals(args[3])) {
			 SortedMap<Integer, Integer> notCalled = v.getNotCalledVariantsLengthDistribution();
			 System.out.println("------------------------------------------------------------------------------------------ ");
			 int val = 0;
			 for(Map.Entry<Integer, Integer> entry : notCalled.entrySet()) {
				 if(entry.getKey() < 2000) {
					 System.out.println(entry.getKey() + "\t" + entry.getValue());
				 }else {
					 val += entry.getValue();
				 }
			 }System.out.println("more" + "\t" + val);
		}
		else v.printReport(ints);
		//List<GenomicVariant> tests = v.getTestVariants();
		//List<GenomicVariant> refs = v.getReferenceVariants();
		//for(int i = 0; i < tests.size(); i++) {
			//System.out.println(" test: " + tests.get(i).getAlleles()[0] + " refs: " + refs.get(i).getAlleles()[0]);
		//}
	}
	public List<GenomicVariant> getReferenceVariants() {
		return referenceVariants;
	}
	public List<GenomicVariant> getTestVariants() {
		return testVariants;
	}
	public int getTotalTruePositiveBps() {
		return truePositiveBps;
	}
	public int getTotalFalsePositiveBps() {
		return falsePositiveBps;
	}
	public void setGenome(String file) throws IOException{
		genome = new ReferenceGenome(file);
	}
	public boolean isSVByLength(int length) {
		return Math.abs(length) >= lengthToDefineSVEvent;
	}
	public boolean isIndel(byte type) {
		return (type == GenomicVariant.TYPE_LARGEDEL || 
				type == GenomicVariant.TYPE_LARGEINS ||
				type == GenomicVariant.TYPE_INDEL);
	}
	public int getMinRefLength() {
		int min = Integer.MAX_VALUE;
		for (GenomicVariant v:referenceVariants) if(Math.abs(v.length()) < min) min = Math.abs(v.length());
		return min;
	}
	public int estimatePositiveBps() {
		int positiveBps = 0;
		for(GenomicVariant ref:referenceVariants) {
			if(isSVByLength(ref.length()) && isIndel(ref.getType())) positiveBps += Math.abs(ref.length());
		}
		return positiveBps;
	}
	
	public int estimatePositiveBps(int treshold) {
		int positivesUnderTreshold = 0;
		for(GenomicVariant ref:referenceVariants) {
			int length = Math.abs(ref.length());
			if(isSVByLength(length) && isIndel(ref.getType()) && length <= treshold) positivesUnderTreshold += length;
		}
		return positivesUnderTreshold;
	}
	
	public void processRefVariantForOtherCalculations(GenomicVariant refVariant) {
		if(isSVByLength(refVariant.length()) && isIndel(refVariant.getType())) {
			referenceVariants.add(refVariant);
			refVariantCount += 1;
		}
	}
	public boolean isIntersection(double percentage) {
		return percentage >= PERCENTAGE_INTERSECTION_TRESHOLD;
	}
	public void updateStatisticsForBaseAnalysis(int truePositiveNucleotides, int falsePositiveNucleotides,
			int truePositiveNucleotidesTreshold, int falsePositiveNucleotidesTreshold) {
		truePositivesPerVar.add(truePositiveNucleotides);
		falsePositivesPerVar.add(falsePositiveNucleotides);	
		truePositiveBps += truePositiveNucleotides;
		falsePositiveBps += falsePositiveNucleotides;
		truePositiveBpsTreshold += truePositiveNucleotidesTreshold;
		falsePositiveBpsTreshold += falsePositiveNucleotidesTreshold;
		System.out.println("Current cumulative true positives: " + truePositiveBps);
		System.out.println("Current cumulative false positives: " + falsePositiveBps);
		System.out.println("Current cumulative true positives under th: " + truePositiveBpsTreshold);
		System.out.println("Current cumulative false positives under th: " + falsePositiveBpsTreshold);
	}
	public List<GenomicVariant> getNotIntersectingReferenceVariants() {
		List<GenomicVariant> notIntersecting = referenceVariants;
		notIntersecting.removeAll(intersectingRefVariants);
		return notIntersecting;
	}
	public SortedMap<Integer, Integer> getNotCalledVariantsLengthDistribution() {
		SortedMap<Integer, Integer> distribution = new TreeMap<>();
		List<GenomicVariant> notIntersecting = getNotIntersectingReferenceVariants();
		for(GenomicVariant var : notIntersecting) {
			int length = Math.abs(var.length());
			distribution.compute(length, (k,v)-> (v == null) ? 1 : v + 1);
		}
		for(int i = 0; i < distribution.lastKey(); i++) {
			if(!distribution.containsKey(i)) distribution.put(i, 0); 
		}
		return distribution;
	}
	public SortedMap<Integer, Integer> getTestVariantsLengthDistribution() {
		SortedMap<Integer, Integer> distribution = new TreeMap<>();
		List<GenomicVariant> testVar = getNotIntersectingReferenceVariants();
		for(GenomicVariant var : testVar) {
			int length = Math.abs(var.length());
			distribution.compute(length, (k,v)-> (v == null) ? 1 : v + 1);
		}
		return distribution;
	}
	
	class VariantIntersectionMetricsCalculator{
		int firstTest;
		int lastTest;
		int firstRef;
		int lastRef;
		int higherLast;
		int lowerLast;
		int higherFirst;
		int lowerFirst;
		int refLength;
		int testLength;
		int span;
		int complement;
		int intersection;
		public VariantIntersectionMetricsCalculator(GenomicVariant test, GenomicVariant ref) {
			firstTest = test.getFirst();
			lastTest = test.getLast();
			firstRef = ref.getFirst();
			lastRef = ref.getLast();
			refLength = lastRef - firstRef + 1;
			testLength = lastTest - firstTest + 1;
			System.out.println("Test -> " + " sequence: " + test.getSequenceName() + " type: " 
					+ GenomicVariantImpl.getVariantTypeName(test.getType()) + " first: " + firstTest 
					+ " last: " + lastTest + " length: " + testLength
					+ " \nRef -> " + " sequence: " + ref.getSequenceName() + " type: " 
					+ GenomicVariantImpl.getVariantTypeName(ref.getType()) + " first: " 
					+ firstRef +  " last: " + lastRef + " length: " + refLength);
		}
		private void estimateLimits() {
			if(lastRef > lastTest) {
				higherLast = lastRef;
				lowerLast = lastTest;
			}
			else {
				higherLast = lastTest;
				lowerLast = lastRef;
			}
			if(firstRef > firstTest) {
				higherFirst = firstRef;
				lowerFirst = firstTest;
			}
			else {
				higherFirst = firstTest;
				lowerFirst = firstRef;
			}
			System.out.println("Higher last: " + higherLast +
					" Higher First: " + higherFirst +
					" Lower last: " + lowerLast +
					" Lower first: " + lowerFirst);
		}
		private void estimateSpan() {
			span =  higherLast - lowerFirst;
			System.out.println("Span: " + span);
		}
		private void estimateComplement() {
			complement = (higherFirst - lowerFirst) + (higherLast - lowerLast);
			System.out.println("Complement: " + complement);
		}
		
		private void estimateIntersection() {
			intersection = span -  complement;
			System.out.println("Intersection: " + intersection);
		}
		public double estimateIntersectionPercentageOverReference() {
			estimateLimits();
			estimateSpan();
			if(span == 0) return 0;
			estimateComplement();
			estimateIntersection();
			System.out.println("Percentage: " + (double) intersection / refLength);
			return (double) intersection / refLength;
		}
		public int estimateTruePositiveNucleotides() {
			System.out.println( "True Positive Nucleotides: " + intersection);
			return intersection;
		}
		public int estimateFalsePositiveNucleotides() {
			System.out.println( "False Positive Nucleotides: " + ( testLength > intersection ? testLength - intersection : 0));
			return testLength > intersection ? testLength - intersection : 0;
		}
		public int estimateTruePositiveNucleotides(int treshold) {
			System.out.println("True Positive Nucleotides under th: " + (testLength <= treshold ? intersection : 0));
			return testLength <= treshold ? intersection : 0;
		}
		public int estimateFalsePositiveNucleotides(int treshold) {
			System.out.println("False Positive Nucleotides under th: " + (testLength > intersection && testLength <= treshold ? testLength - intersection : 0));
			return testLength > intersection && testLength <= treshold ? testLength - intersection : 0;
		}
	}
}

	
