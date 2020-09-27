package ngsep.simulation;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.discovery.TillingPoolsIndividualGenotyper;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class TillingIndividualVCF2PoolVCF {

	public static void main(String[] args) throws Exception {
		TillingIndividualVCF2PoolVCF instance = new TillingIndividualVCF2PoolVCF();
		String individualVCF = args[0];
		String poolsFile = args[1];
		Map<String,List<Integer>> poolsData = TillingPoolsIndividualGenotyper.loadPoolsFile(poolsFile);
		instance.processIndividualsVCF (individualVCF,poolsData,System.out);
	}

	public void processIndividualsVCF(String individualVCF, Map<String, List<Integer>> poolsData, PrintStream out) throws IOException {
		System.err.println("Loaded "+poolsData.size()+" samples");
		Map<Integer, Sample> poolSamples = createPoolSamples(poolsData);
		System.err.println("Loaded "+poolSamples.size()+" pools");
		List<Sample> poolSamplesList = new ArrayList<Sample>(poolSamples.size());
		poolSamplesList.addAll(poolSamples.values());
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		VCFFileWriter writer = new VCFFileWriter();
		header.addSamples(poolSamplesList);
		writer.printHeader(header, out);
		try (VCFFileReader reader = new VCFFileReader(individualVCF)) {
			Iterator<VCFRecord> it = reader.iterator();
			while(it.hasNext()) {
				VCFRecord indvRecord = it.next();
				GenomicVariant variant = indvRecord.getVariant();
				if(!variant.isBiallelic()) continue;
				List<CalledGenomicVariant> indvCalls = indvRecord.getCalls();
				Map<Integer,Integer> poolAltCounts = calculateAlternativeCounts(poolsData, indvCalls);
				List<CalledGenomicVariant> poolCalls = new ArrayList<CalledGenomicVariant>(poolSamplesList.size());
				int i=0;
				for(int poolId:poolSamples.keySet()) {
					Sample sample = poolSamplesList.get(i);
					Integer count = poolAltCounts.get(poolId);
					CalledGenomicVariantImpl call;
					short [] allelesCN = new short[2];
					if(count==null) {
						call = new CalledGenomicVariantImpl(variant, CalledGenomicVariant.GENOTYPE_HOMOREF);
						allelesCN[0] = sample.getNormalPloidy();
						allelesCN[1] = 0;
					} else {
						call = new CalledGenomicVariantImpl(variant, CalledGenomicVariant.GENOTYPE_HETERO);
						allelesCN[0] = (short) (sample.getNormalPloidy()-count);
						allelesCN[1] = (short)count.intValue();
					}
					call.setGenotypeQuality((short) 100);
					call.setSampleId(sample.getId());
					call.setAllelesCopyNumber(allelesCN);
					poolCalls.add(call);
					i++;
				}
				VCFRecord outRecord = new VCFRecord(variant,VCFRecord.DEF_FORMAT_ARRAY_QUAL_CN, poolCalls,header);
				writer.printVCFRecord(outRecord, out);
			}
		}	
	}

	private Map<Integer, Sample> createPoolSamples(Map<String, List<Integer>> poolsData) {
		Map<Integer,Sample> poolSamples = new TreeMap<Integer,Sample>();
		for(String indvSampleId:poolsData.keySet()) {
			List<Integer> poolsSample = poolsData.get(indvSampleId);
			for(int i:poolsSample) {
				Sample s = poolSamples.get(i);
				if(s==null) {
					s = new Sample(""+i);
					poolSamples.put(i, s);
				} else {
					s.setNormalPloidy((short) (s.getNormalPloidy()+2));
				}
			}
		}
		return poolSamples;
	}

	private Map<Integer, Integer> calculateAlternativeCounts(Map<String, List<Integer>> poolsData, List<CalledGenomicVariant> indvCalls) {
		Map<Integer, Integer> poolCounts = new HashMap<Integer, Integer>();
		for(CalledGenomicVariant call:indvCalls) {
			if(call.isUndecided() || call.isHomozygousReference()) continue;
			int altAlleles = call.isHeterozygous()?1:2;
			String sampleId = call.getSampleId();
			List<Integer> pools = poolsData.get(sampleId);
			for(int poolId:pools) {
				poolCounts.compute(poolId, (k,v)->v!=null?v+altAlleles:altAlleles);
			}
		}
		return poolCounts;
	}

}
