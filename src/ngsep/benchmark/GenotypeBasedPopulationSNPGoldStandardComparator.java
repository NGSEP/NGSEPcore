package ngsep.benchmark;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ngsep.genome.GenomicRegionComparator;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.SimpleSequenceListLoader;
import ngsep.variants.CalledGenomicVariant;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFRecord;

public class GenotypeBasedPopulationSNPGoldStandardComparator {

	private QualifiedSequenceList sequenceNames;
	private int [][] variantCounts = new int [10][3];
	private GoldStandardComparisonCounts genotypingCounts = new GoldStandardComparisonCounts();
	public static void main(String[] args) throws Exception {
		GenotypeBasedPopulationSNPGoldStandardComparator instance = new GenotypeBasedPopulationSNPGoldStandardComparator();
		SimpleSequenceListLoader loader = new SimpleSequenceListLoader();
		instance.sequenceNames = loader.loadSequences(args[0]);
		instance.compareVCFs(args[1],args[2]);
		instance.printResults(System.out);
	}

	public void compareVCFs(String gsVCF, String testVCF) throws IOException {
		genotypingCounts.setCountNonGSAsFP(true);
		GenomicRegionComparator comparator = new GenomicRegionComparator(sequenceNames);
		try (VCFFileReader r1 = new VCFFileReader(gsVCF);
			 VCFFileReader r2 = new VCFFileReader(testVCF)) {
			Iterator<VCFRecord> it1 = r1.iterator();
			Iterator<VCFRecord> it2 = r2.iterator();
			VCFRecord next1 = it1.next();
			VCFRecord next2 = it2.next();
			while (next1!=null && next2!=null) {
				int cmp = comparator.compare(next1, next2);
				if (cmp == 0) {
					compareRecords(next1,next2);
					if(it1.hasNext()) next1 = it1.next();
					else next1=null;
					if(it2.hasNext()) next2 = it2.next();
					else next2=null;
				} else if (cmp<0) {
					registerDifferences(next1,1);
					if(it1.hasNext()) next1 = it1.next();
					else next1=null;
				} else {
					registerDifferences(next2,2);
					if(it2.hasNext()) next2 = it2.next();
					else next2=null;
				}
			}
		}
		
	}

	private void compareRecords(VCFRecord recordGoldStandard, VCFRecord recordTest) {
		int lastRowCounts = GoldStandardComparisonCounts.NUM_ROWS_COUNTS-1;
		Map<String,CalledGenomicVariant> callsMapGS = createCallsMap(recordGoldStandard);
		Map<String,CalledGenomicVariant> callsMapTest = createCallsMap(recordTest);
		int rowV = Math.min(variantCounts.length-1, recordTest.getVariant().getVariantQS()/10);
		for(int i=0;i<=rowV;i++) {
			(variantCounts[i][0])++;
		}
		if(!recordGoldStandard.getVariant().getReference().equals(recordTest.getVariant().getReference())) {
			System.err.println("Inconsistent refrence for SNP at "+recordGoldStandard.getSequenceName()+":"+recordGoldStandard.getFirst()+" ref GS "+recordGoldStandard.getVariant().getReference()+" ref test: "+recordTest.getVariant().getReference());
		}
		for(String sampleId:callsMapGS.keySet()) {
			CalledGenomicVariant callGS = callsMapGS.get(sampleId);
			CalledGenomicVariant callTest = callsMapTest.get(sampleId);
			
			if(callGS==null || callTest==null) continue;
			int genotypeGS = VCFGoldStandardComparator.getGenotypeNumber(callGS);
			int genotypeTest = VCFGoldStandardComparator.getGenotypeNumber(callTest);
			short qualTest = VCFGoldStandardComparator.loadGenotypeQuality(callTest);
			if(genotypeGS==CalledGenomicVariant.GENOTYPE_UNDECIDED) continue;
			int k=9+genotypeGS;
			if(genotypeTest!=CalledGenomicVariant.GENOTYPE_UNDECIDED)  k = 3*genotypeGS+genotypeTest;
			int row = qualTest>0?Math.min(qualTest/10, lastRowCounts):-1; 
			genotypingCounts.update(0,row,k);
			genotypingCounts.update(row+1,lastRowCounts,9+genotypeGS);
		}
	}

	private Map<String, CalledGenomicVariant> createCallsMap(VCFRecord record) {
		List<CalledGenomicVariant> calls = record.getCalls();
		Map<String, CalledGenomicVariant> callsMap = new HashMap<String, CalledGenomicVariant>(calls.size());
		for(CalledGenomicVariant call:calls) {
			callsMap.put(call.getSampleId(), call);
		}
		return callsMap;
	}

	private void registerDifferences(VCFRecord record, int j) {
		List<CalledGenomicVariant> calls = record.getCalls();
		
		int rowV = Math.min(variantCounts.length-1, record.getVariant().getVariantQS()/10);
		for(int i=0;i<=rowV;i++) {
			(variantCounts[i][j])++;
		}
		for(CalledGenomicVariant call:calls) {
			if(call.isUndecided()) continue;
			int genotypeNumber = VCFGoldStandardComparator.getGenotypeNumber(call);
			int columnIndex = 9;
			if(j==2) columnIndex = 12;
			columnIndex+=genotypeNumber;
			int row = GoldStandardComparisonCounts.NUM_ROWS_COUNTS-1;
			
			if(j==2) {
				short qualTest = VCFGoldStandardComparator.loadGenotypeQuality(call);
				row = qualTest>0?Math.min(qualTest/10, row):-1;
			}
			genotypingCounts.update(0,row,columnIndex);
		}
	}
	
	private void printResults(PrintStream out) {
		out.println("VARIANT COUNTS");
		for(int i=0;i<variantCounts.length;i++) {
			out.print(""+(10*i));
			for(int j=0;j<variantCounts[i].length;j++) {
				out.print("\t"+variantCounts[i][j]);
			}
			out.println();
		}
		out.println();
		out.println("GENOTYPING COUNTS");
		genotypingCounts.print(out);
		
	}

}
