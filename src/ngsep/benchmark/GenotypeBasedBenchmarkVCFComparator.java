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

public class GenotypeBasedBenchmarkVCFComparator {

	private QualifiedSequenceList sequenceNames;
	private int [][] genotypingCounts = new int [10][4];
	private int [][] falseVariantCounts = new int [10][4];
	public static void main(String[] args) throws Exception {
		GenotypeBasedBenchmarkVCFComparator instance = new GenotypeBasedBenchmarkVCFComparator();
		SimpleSequenceListLoader loader = new SimpleSequenceListLoader();
		instance.sequenceNames = loader.loadSequences(args[0]);
		instance.compareVCFs(args[1],args[2]);
		instance.printResults(System.out);
	}

	public void compareVCFs(String vcfFilename1, String vcfFilename2) throws IOException {
		GenomicRegionComparator comparator = new GenomicRegionComparator(sequenceNames);
		try (VCFFileReader r1 = new VCFFileReader(vcfFilename1);
			 VCFFileReader r2 = new VCFFileReader(vcfFilename2)) {
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
					registerDifferences(next1,0);
					if(it1.hasNext()) next1 = it1.next();
					else next1=null;
				} else {
					registerDifferences(next2,1);
					if(it2.hasNext()) next2 = it2.next();
					else next2=null;
				}
			}
		}
		
	}

	private void compareRecords(VCFRecord record1, VCFRecord record2) {
		Map<String,CalledGenomicVariant> callsMap1 = createCallsMap(record1);
		Map<String,CalledGenomicVariant> callsMap2 = createCallsMap(record2);
		for(String sampleId:callsMap1.keySet()) {
			CalledGenomicVariant call1 = callsMap1.get(sampleId);
			CalledGenomicVariant call2 = callsMap2.get(sampleId);
			//System.out.println("Variant 1 "+record1.getSequenceName()+":"+record1.getFirst()+" Variant 2 "+record2.getSequenceName()+":"+record2.getFirst()+" call1 "+call1+" call2 "+call2);
			if(call1==null || call2==null) continue;
			if(call1.isUndecided() || call2.isUndecided()) continue;
			int j = 0;
			if(call1.isHomozygous() && call2.isHeterozygous()) j=1;
			if(call1.isHeterozygous() && call2.isHomozygous()) j=2;
			String [] a1 = call1.getCalledAlleles();
			String [] a2 = call2.getCalledAlleles();
			if(!a1[0].equals(a2[0])) j=3;
			int row = Math.min(call1.getGenotypeQuality(),call2.getGenotypeQuality())/10;
			if(row>=genotypingCounts.length) row = genotypingCounts.length-1;
			for(int i=0;i<=row;i++) {
				(genotypingCounts[i][j])++;
			}
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
		for(CalledGenomicVariant call:calls) {
			if(call.isUndecided()) continue;
			if(call.isHomozygousReference()) continue;
			int columnIndex = 2*j;
			if(!call.isHeterozygous()) columnIndex++;
			int row = call.getGenotypeQuality()/10;
			if(row>=falseVariantCounts.length) row = falseVariantCounts.length-1;
			for(int i=0;i<=row;i++) {
				(falseVariantCounts[i][columnIndex])++;
			}
		}
	}
	
	private void printResults(PrintStream out) {
		for(int i=0;i<genotypingCounts.length;i++) {
			out.print(""+(10*i));
			for(int j=0;j<genotypingCounts[i].length;j++) {
				out.print("\t"+genotypingCounts[i][j]);
			}
			for(int j=0;j<falseVariantCounts[i].length;j++) {
				out.print("\t"+falseVariantCounts[i][j]);
			}
			out.println();
		}
		
	}

}
