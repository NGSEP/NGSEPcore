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
	
	public VCFDirectSVComparison() {
		referenceVariants = new ArrayList<>();
		testVariants = new ArrayList<>();
	}
	
	public List<GenomicVariant> getReferenceVariants() {
		return referenceVariants;
	}

	public List<GenomicVariant> getTestVariants() {
		return testVariants;
	}
	
	public void setGenome(String file) throws IOException{
		genome = new ReferenceGenome(file);
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
			int ref = 0;
			while(testRecords != null && refRecords != null) {
				/*while(refRecords.getVariant().isSNV() && it2.hasNext()) { 
					refRecords = it2.next();
				}
				/**if(testVariant.getFirst() - firstTest < -1000) {
					System.out.println("Cambio de cromosoma: test");
					while(refVariant.getFirst() >= firstRef) {
						firstRef = refVariant.getFirst();
						if(it2.hasNext())refRecords = it2.next();
						else break;
						refVariant = refRecords.getVariant();
					}
				}else if(refVariant.getFirst() - firstRef < -1000) {
					System.out.println("Cambio de cromosoma: referencia");
					while(testVariant.getFirst() >= firstTest) {
						firstTest = testVariant.getFirst();
						if(it1.hasNext()) testRecords = it1.next();
						else break;
						testVariant = testRecords.getVariant();
					}
				}**/
				testVariant = (GenomicVariantImpl) testRecords.getVariant();
				refVariant = (GenomicVariantImpl) refRecords.getVariant();
				firstTest = testVariant.getFirst();
				firstRef = refVariant.getFirst();
				lastTest = testVariant.getLast();
				lastRef = refVariant.getLast();
				//int cmp = compare(firstTest, firstRef, lastTest, lastRef);
				if(lastRef-firstRef < 2) {
					lastRef = firstRef + refVariant.length()- 1;
					refVariant.setLast(lastRef);
				}
				if(lastTest-firstTest < 2) {
					lastTest = firstTest + testVariant.length() - 1;
					testVariant.setLast(lastTest);
				}
				int cmp = comparator.compare(testVariant, refVariant);
				//System.out.println(testRecords.getVariant().length() + " " + (lastTest-firstTest + 1) );
				//System.out.println("Comparison: "+cmp + " for: " + 
					//	"firstTest: "+ firstTest + " lastTest: " + lastTest +
						//" firstRef: " + firstRef + " lastRef: " + lastRef);
				if(cmp > 1) {
					if(it2.hasNext()) {
						refRecords = it2.next();
						//ref++;
					}
					else break;
				}
				else if(cmp < -1) {
					if(it1.hasNext()) testRecords = it1.next();
					else break;
				}
				else {
					if(intersectionTreshold(firstTest, firstRef,
							lastTest, lastRef)) {
						intersectingVariants++;
						testVariants.add(testVariant);
						referenceVariants.add(refVariant);
						//ref++;
					}
					if(it1.hasNext()) testRecords = it1.next();
					else break;
					if(it2.hasNext())refRecords = it2.next();
					else break;
				}
			}//System.out.println(" La cantidad de referencia: " + ref);
			return intersectingVariants;
		}
	}
	
	/**public int compare(int first1, int first2, int last1, int last2) {
		if(last1 < first2) return -1;
		else if(last2 < first1) return 1;
		else return 0;
	}**/
	
	public boolean intersectionTreshold(int first1, int first2, int last1, int last2) {
		int higherLast = 0;
		int lowerLast = 0;
		int higherFirst = 0;
		int lowerFirst = 0;
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
		if(span == 0) return false;
		int complement = (higherFirst - lowerFirst) + (higherLast - lowerLast);
		int intersection = span - complement;
		double percentage = (double) intersection / span;
		//System.out.println(" Intersection: " + intersection + " union: " + span + " percentage: " + percentage);
		return percentage >= PERCENTAGE_INTERSECTION_TRESHOLD;
	}
	
	public static void main(String[] args) throws IOException {
		VCFDirectSVComparison v = new VCFDirectSVComparison();
		v.setGenome(args[0]);
		String testFile = args[1];
		String refFile = args[2];
		int ints = v.findIntersectingVariants(testFile, refFile);
		System.out.println(ints);
		//List<GenomicVariant> tests = v.getTestVariants();
		//List<GenomicVariant> refs = v.getReferenceVariants();
		//for(int i = 0; i < tests.size(); i++) {
			//System.out.println(" test: " + tests.get(i).getAlleles()[0] + " refs: " + refs.get(i).getAlleles()[0]);
		//}
	}
}
	
