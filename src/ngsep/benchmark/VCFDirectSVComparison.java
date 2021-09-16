package ngsep.benchmark;

import java.io.IOException;
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
	
	private ReferenceGenome genome;
	private List<GenomicVariant> referenceVariants;
	private List<GenomicVariant> testVariants;
	private List<Integer> truePositivesPerVar;
	private List<Integer> falsePositivesPerVar;
	private int lengthToDefineSVEvent = 30;
	private int testVariantCount;
	private int refVariantCount;
	private int truePositiveBps;
	private int falsePositiveBps;
	private int positiveBps;
	
	public VCFDirectSVComparison() {
		referenceVariants = new ArrayList<>();
		testVariants = new ArrayList<>();
		truePositivesPerVar = new ArrayList<>();
		falsePositivesPerVar = new ArrayList<>();
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
	
	public int addIfIndel(int addAmount, byte type) {
		addAmount = Math.abs(addAmount);
		return (type == GenomicVariant.TYPE_LARGEDEL || 
				type == GenomicVariant.TYPE_LARGEINS) ? addAmount : 0;
	}
	public int getMinRefLength() {
		int min = Integer.MAX_VALUE;
		for (GenomicVariant v:referenceVariants) if(Math.abs(v.length()) < min) min = Math.abs(v.length());
		return min;
	}
	public int findIntersectingVariants(String vcf1, String vcf2) throws IOException{
		try(VCFFileReader input1 = new VCFFileReader(vcf1);
				VCFFileReader input2 = new VCFFileReader(vcf2);){
			AlleleCompatibilityGenomicVariantComparator comparator = new AlleleCompatibilityGenomicVariantComparator(genome);
			input1.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			input2.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			Iterator<VCFRecord> it1 = input1.iterator();
			Iterator<VCFRecord> it2 = input2.iterator();
			int intersectingVariants = 0;
			VCFRecord testRecords = it1.next();
			VCFRecord refRecords = it2.next();
			GenomicVariantImpl testVariant = (GenomicVariantImpl) testRecords.getVariant();
			GenomicVariantImpl refVariant = (GenomicVariantImpl) refRecords.getVariant();
			int firstTest = -1;
			int firstRef = -1;
			int lastTest = -1;
			int lastRef = -1;
			testVariantCount = 1;
			if(Math.abs(refRecords.length()) >= lengthToDefineSVEvent) {
				refVariantCount += addIfIndel( 1, refVariant.getType());
				positiveBps += addIfIndel(refVariant.length(), refVariant.getType());
			}
			while(testRecords != null && refRecords != null) {
				testVariant = (GenomicVariantImpl) testRecords.getVariant();
				refVariant = (GenomicVariantImpl) refRecords.getVariant();
				firstTest = testVariant.getFirst();
				firstRef = refVariant.getFirst();
				lastTest = testVariant.getLast();
				lastRef = refVariant.getLast();
				if(lastRef-firstRef < 2) {
					lastRef = firstRef + refVariant.length()- 1;
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
						referenceVariants.add(refVariant);
						if(isSVByLength(refRecords.length())) {
							positiveBps += addIfIndel(refVariant.length(), refVariant.getType());
							refVariantCount += addIfIndel(1, refVariant.getType());
						}	
					}
					else break;
				}
				else if(cmp < -1) {
					if(it1.hasNext()) {
						testRecords = it1.next();
						testVariantCount++;
					}
					else break;
				}
				else {
					if(isIntersection(firstTest, firstRef,
							lastTest, lastRef)) {
						if(testVariant.getType() == refVariant.getType() &&
								isSVByLength(refRecords.length())) {
							intersectingVariants++;
						}
						testVariants.add(testVariant);
						referenceVariants.add(refVariant);
					}
					if(it1.hasNext()) {
						testRecords = it1.next();
						testVariantCount++;
					}
					else break;
					if(it2.hasNext()) {
						refRecords = it2.next();
						referenceVariants.add(refVariant);
						if(isSVByLength(refRecords.length())) {
							positiveBps += addIfIndel(refVariant.length(), refVariant.getType());
							refVariantCount += addIfIndel(1, refVariant.getType());
						}
					}
					else break;
				}
			}
			return intersectingVariants;
		}
	}
	public int estimateTruePositiveNucleotides(int intersection) {
		return intersection;
	}
	public int estimateFalsePositiveNucleotides(int testLength, int intersection) {
		return testLength - intersection;
	}
	public double estimateIntersectionPercentageOverReference(int first1, int first2, int last1, int last2) {
		int higherLast = 0;
		int lowerLast = 0;
		int higherFirst = 0;
		int lowerFirst = 0;
		int refLength = last2 - first2 + 1;
		int testLength = last1 - first1 + 1;
		if(last2 > last1) {
			higherLast = last2;
			lowerLast = last1;
		}
		else {
			higherLast = last1;
			lowerLast = last2;
		}
		if(first2 > first1) {
			higherFirst = first2;
			lowerFirst = first1;
		}
		else {
			higherFirst = first1;
			lowerFirst = first2;
		}
		int span = higherLast - lowerFirst;
		if(span == 0) return 0;
		int complement = (higherFirst - lowerFirst) + (higherLast - lowerLast);
		int intersection = span - complement;
		int truePositiveNucleotides = estimateTruePositiveNucleotides(intersection);
		int falsePositiveNucleotides = estimateFalsePositiveNucleotides(testLength, intersection);
		truePositivesPerVar.add(truePositiveNucleotides);
		falsePositivesPerVar.add(falsePositiveNucleotides);
		truePositiveBps += truePositiveNucleotides;
		falsePositiveBps += falsePositiveNucleotides;
		double percentage = (double) intersection / refLength;
		//System.out.println(" Intersection: " + intersection + " union: " + refLength + " percentage: " + percentage);
		return percentage;
	}
	public boolean isIntersection(int first1, int first2, int last1, int last2) {
		double percentage = estimateIntersectionPercentageOverReference(first1,first2,last1,last2);
		return percentage >= PERCENTAGE_INTERSECTION_TRESHOLD;
	}
	public void printReport(String testFile, String refFile) throws IOException {
		int ints = findIntersectingVariants(testFile, refFile);
		System.out.println("Number of test variants: " + testVariantCount);
		System.out.println("Number of reference variants: " + refVariantCount);
		System.out.println("Intersections: " + ints);
		System.out.println("Recall (sensitivity) over intersections: " +  ((double) ints / refVariantCount)*100);
		double precision = (double) ints / (testVariantCount);
		System.out.println("Precision over intersections: " + precision*100);
		System.out.println("FDR over intersections: " + ((1 - precision)*100));
		System.out.println("Total True positive bases: " + truePositiveBps);
		System.out.println("Total False positive bases: " + falsePositiveBps);
		double perBasePrecision = (double) truePositiveBps / (truePositiveBps + falsePositiveBps);
		System.out.println("Recall (sensitivity) over bases: " +  ((double) truePositiveBps / positiveBps)*100);
		System.out.println("Precision over bases: " + perBasePrecision*100);
		System.out.println("FDR over bases: " + ((1 - perBasePrecision)*100));
	}
	public static void main(String[] args) throws IOException {
		VCFDirectSVComparison v = new VCFDirectSVComparison();
		v.setGenome(args[0]);
		String testFile = args[1];
		String refFile = args[2];
		v.printReport(testFile, refFile);
		//List<GenomicVariant> tests = v.getTestVariants();
		//List<GenomicVariant> refs = v.getReferenceVariants();
		//for(int i = 0; i < tests.size(); i++) {
			//System.out.println(" test: " + tests.get(i).getAlleles()[0] + " refs: " + refs.get(i).getAlleles()[0]);
		//}
	}
}
	
