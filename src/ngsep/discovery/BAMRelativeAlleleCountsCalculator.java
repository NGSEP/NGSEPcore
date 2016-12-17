package ngsep.discovery;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.math.Distribution;




public class BAMRelativeAlleleCountsCalculator implements PileupListener {

	private Distribution distProp = new Distribution(0, 0.5, 0.01);
	private Distribution distNumAlleles = new Distribution(0, 10, 1);
	private GenomicRegionSortedCollection<GenomicRegion> repeats;
	private GenomicRegionSortedCollection<GenomicRegion> selectedRegions;
	
	private int minCoverage = 10;
	private int maxCoverage = 1000;
	private int minBaseQualityScore = 20;
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		BAMRelativeAlleleCountsCalculator instance = new BAMRelativeAlleleCountsCalculator();
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		generator.addListener(instance);
		String repeatsFile = null;
		String selectedRegionsFile = null;
		int i=0;
		while(i<args.length && args[i].charAt(0)=='-') {
			if("-m".equalsIgnoreCase(args[i])) {
				i++;
				instance.minCoverage = Integer.parseInt(args[i]);
			} else if("-M".equalsIgnoreCase(args[i])) {
				i++;
				instance.maxCoverage = Integer.parseInt(args[i]);
			} else if("-q".equalsIgnoreCase(args[i])) {
				i++;
				instance.minBaseQualityScore = Integer.parseInt(args[i]);
			} else if("-r".equalsIgnoreCase(args[i])) {
				i++;
				repeatsFile = args[i];
			}  else if("-f".equalsIgnoreCase(args[i])) {
				i++;
				selectedRegionsFile = args[i];
			} else if ("-s".equalsIgnoreCase(args[i])) {
				generator.setProcessSecondaryAlignments(true);
			}
			i++;
		}
		String filename = args[i++];
		generator.setMaxAlnsPerStartPos(instance.maxCoverage);
		if(repeatsFile!=null) {
			SimpleGenomicRegionFileHandler grfh = new SimpleGenomicRegionFileHandler();
			List<GenomicRegion> repeats = grfh.loadRegions(repeatsFile);
			instance.repeats = new GenomicRegionSortedCollection<GenomicRegion>();
			instance.repeats.addAll(repeats);
		}
		if(selectedRegionsFile!=null) {
			SimpleGenomicRegionFileHandler grfh = new SimpleGenomicRegionFileHandler();
			List<GenomicRegion> regionsList = grfh.loadRegions(selectedRegionsFile);
			instance.selectedRegions = new GenomicRegionSortedCollection<GenomicRegion>();
			instance.selectedRegions.addAll(regionsList);
		}
		generator.processFile(filename);
		System.out.println("Distribution of allele proportions");
		instance.distProp.printDistribution(System.out);
		System.out.println("Distribution of number of alleles");
		instance.distNumAlleles.printDistribution(System.out);
	}

	@Override
	public void onPileup(PileupRecord pileup) {
		if(repeats!=null) {
			GenomicRegionSortedCollection<GenomicRegion> spanningRepeats = repeats.findSpanningRegions(pileup.getSequenceName(), pileup.getPosition());
			if(spanningRepeats.size()>0) return;
		}
		if(selectedRegions!=null) {
			GenomicRegionSortedCollection<GenomicRegion> spanningRegions = selectedRegions.findSpanningRegions(pileup.getSequenceName(), pileup.getPosition());
			if(spanningRegions.size()==0) return;
		}
		List<String> calls = pileup.getAlleleCalls(1);
		if(calls.size()<2*minCoverage) return;
		Map<String, Integer> alleleCounts = new TreeMap<String, Integer>();
		for(int i=0;i<calls.size();i+=2) {
			String call = calls.get(i);
			char qsC = calls.get(i+1).charAt(0);
			int qs = ((byte)qsC)-33;
			if(qs>=minBaseQualityScore) {
				Integer count = alleleCounts.get(call);
				if(count==null) {
					alleleCounts.put(call, 1);
				} else {
					alleleCounts.put(call, count+1);
				}
			}
			
		}
		distNumAlleles.processDatapoint(alleleCounts.size());
		String alleleMax = getMaximum (alleleCounts,null);
		
		int countMax = 0;
		int countSecondMax = 0;
		if(alleleMax!=null) {
			countMax = alleleCounts.get(alleleMax);
			String alleleSecondMax = getMaximum (alleleCounts,alleleMax);
			if(alleleSecondMax!=null)  countSecondMax = alleleCounts.get(alleleSecondMax);
		}
		if(countMax>0) {
			double prop = (double)countSecondMax/(countMax+countSecondMax);
			distProp.processDatapoint(prop);
		}
	}

	private String getMaximum(Map<String, Integer> alleleCounts,String alleleIgnore) {
		String max = null;
		int maxCount = -1;
		for(String allele:alleleCounts.keySet()) {
			if(allele.equals(alleleIgnore)) continue;
			int count = alleleCounts.get(allele);
			if(max == null || maxCount<count) {
				max = allele;
				maxCount = count;
			}
			
		}
		return max;
	}

	@Override
	public void onSequenceStart(String sequenceName) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void onSequenceEnd(String sequenceName) {
		// TODO Auto-generated method stub
		
	}
	

}
