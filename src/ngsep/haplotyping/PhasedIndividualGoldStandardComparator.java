package ngsep.haplotyping;

import java.io.IOException;
import java.util.Iterator;

import ngsep.genome.GenomicRegionComparator;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.SimpleSequenceListLoader;
import ngsep.variants.CalledGenomicVariant;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class PhasedIndividualGoldStandardComparator {

	private QualifiedSequenceList sequenceNames;
	
	public static void main(String[] args) throws Exception {
		PhasedIndividualGoldStandardComparator instance = new PhasedIndividualGoldStandardComparator();
		instance.loadSequenceNames(args[0]);
		instance.compareVCFs(args[1],args[2]);

	}

	private void loadSequenceNames(String filename) throws IOException {
		SimpleSequenceListLoader loader = new SimpleSequenceListLoader();
		sequenceNames = loader.loadSequences(filename);
	}

	private void compareVCFs(String goldStandard, String testFile) throws IOException {
		int switchErrors=0;
		int totalVariants = 0;
		GenomicRegionComparator comparator = new GenomicRegionComparator(sequenceNames);
		try (VCFFileReader r1 = new VCFFileReader(goldStandard);
				VCFFileReader r2 = new VCFFileReader(testFile);	) {
			Iterator<VCFRecord> it1 = r1.iterator();
			Iterator<VCFRecord> it2 = r2.iterator();
			Boolean sameMolecule =null; 
			String seqName = null;
			VCFRecord recordGS;
			if(it1.hasNext()) recordGS = it1.next();
			else recordGS = null;
			VCFRecord recordTest;
			if(it2.hasNext()) recordTest = it2.next();
			else recordTest = null;
			while(recordGS!=null && recordTest!=null) {
				int cmp = comparator.compare(recordGS, recordTest);
				if(cmp < 0) {
					if(it1.hasNext()) recordGS = it1.next();
					else recordGS = null;
				} else if (cmp>0) {
					if(it2.hasNext()) recordTest = it2.next();
					else recordTest = null;
				} else {
					if(!recordGS.getSequenceName().equals(seqName)) {
						sameMolecule = null;
						seqName = recordGS.getSequenceName();
					}
					CalledGenomicVariant callGS = recordGS.getCalls().get(0);
					CalledGenomicVariant callTest = recordTest.getCalls().get(0);
					if(callGS.isBiallelic() && callGS.isHeterozygous() && callTest.isBiallelic() && callTest.isHeterozygous()) {
						String [] phasedAllelesGS = callGS.getPhasedAlleles();
						String [] phasedAllelesTest = callTest.getPhasedAlleles();
						if(phasedAllelesGS.length == 2 && phasedAllelesTest.length==2) {
							totalVariants++;
							boolean coincide = phasedAllelesGS[0].equals(phasedAllelesTest[0]);
							if(sameMolecule==null) {
								sameMolecule = coincide;
							} else if (sameMolecule!=coincide) {
								System.out.println("Switch error at "+seqName+":"+callGS.getFirst());
								switchErrors++;
								sameMolecule = coincide;
							}
						}
					}
					
					if(it1.hasNext()) recordGS = it1.next();
					else recordGS = null;
					if(it2.hasNext()) recordTest = it2.next();
					else recordTest = null;
				}
			}
		}
		System.out.println("SE: "+switchErrors+" variants: "+totalVariants+" rate: "+(double)switchErrors/totalVariants);
		
	}

}
