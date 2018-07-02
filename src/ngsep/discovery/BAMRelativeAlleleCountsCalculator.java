package ngsep.discovery;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.math.Distribution;
import ngsep.sequences.QualifiedSequence;




public class BAMRelativeAlleleCountsCalculator implements PileupListener {

	private Logger log = Logger.getLogger(BAMRelativeAlleleCountsCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private long coveredGenomeSize = 0;
	private AlignmentsPileupGenerator generator;
	private GenomicRegionSortedCollection<GenomicRegion> repeats;
	private GenomicRegionSortedCollection<GenomicRegion> selectedRegions;
	
	
	private int minRD = 10;
	private int maxRD = 1000;
	private int minBaseQualityScore = 20;
	private String repeatsFile = null;
	private String selectedRegionsFile = null;
	private boolean secondaryAlns = false;
	
	private Distribution distProp = new Distribution(0, 0.5, 0.01);
	private List<String> sequenceNamesList = new ArrayList<>();
	private List<Distribution> distPropPerSequenceList = new ArrayList<>();
	private Distribution currentSequencePropDist=null;
	private Distribution distNumAlleles = new Distribution(1, 10, 1);
	
	private static DecimalFormat format = new DecimalFormat("##0.0#");
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		BAMRelativeAlleleCountsCalculator instance = new BAMRelativeAlleleCountsCalculator();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		
		String filename = args[i++];
		instance.runProcess(filename);
		instance.printResults(System.out);
	}
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
	
	public int getMinRD() {
		return minRD;
	}

	public void setMinRD(int minRD) {
		this.minRD = minRD;
	}
	
	public void setMinRD(Integer minRD) {
		this.setMinRD(minRD.intValue());
	}

	public int getMaxRD() {
		return maxRD;
	}

	public void setMaxRD(int maxRD) {
		this.maxRD = maxRD;
	}

	public void setMaxRD(Integer maxRD) {
		this.setMaxRD(maxRD.intValue());
	}
	
	public int getMinBaseQualityScore() {
		return minBaseQualityScore;
	}

	public void setMinBaseQualityScore(int minBaseQualityScore) {
		this.minBaseQualityScore = minBaseQualityScore;
	}
	
	public void setMinBaseQualityScore(Integer minBaseQualityScore) {
		this.setMinBaseQualityScore(minBaseQualityScore.intValue());
	}

	public String getRepeatsFile() {
		return repeatsFile;
	}

	public void setRepeatsFile(String repeatsFile) {
		this.repeatsFile = repeatsFile;
	}

	public String getSelectedRegionsFile() {
		return selectedRegionsFile;
	}
	
	public void setSelectedRegionsFile(String selectedRegionsFile) {
		this.selectedRegionsFile = selectedRegionsFile;
	}

	public boolean isSecondaryAlns() {
		return secondaryAlns;
	}

	public void setSecondaryAlns(boolean secondaryAlns) {
		this.secondaryAlns = secondaryAlns;
	}

	public void setSecondaryAlns(Boolean secondaryAlns) {
		this.setSecondaryAlns(secondaryAlns.booleanValue());
	}
	
	public void runProcess (String filename) throws IOException {
		coveredGenomeSize = 0;
		generator = new AlignmentsPileupGenerator();
		generator.addListener(this);
		generator.setLog(log);
		generator.setProcessSecondaryAlignments(secondaryAlns);
		generator.setMaxAlnsPerStartPos(maxRD);
		if(repeatsFile!=null) {
			repeats = new GenomicRegionSortedCollection<GenomicRegion>();
			SimpleGenomicRegionFileHandler grfh = new SimpleGenomicRegionFileHandler();
			repeats.addAll(grfh.loadRegions(repeatsFile));
		}
		if(selectedRegionsFile!=null) {
			SimpleGenomicRegionFileHandler grfh = new SimpleGenomicRegionFileHandler();
			selectedRegions = new GenomicRegionSortedCollection<GenomicRegion>();
			selectedRegions.addAll(grfh.loadRegions(selectedRegionsFile));
		}
		generator.processFile(filename);
		
	}
	public void printResults (PrintStream out) {
		out.println("Distribution of allele proportions");
		distProp.printDistribution(out);
		out.println("Distribution of number of alleles");
		distNumAlleles.printDistributionInt(out);
		
		if(sequenceNamesList.size()==0) return;
		out.println("Distribution of allele proportions per sequence");
		out.print("Proportion");
		for(int i=0;i<sequenceNamesList.size();i++) {
			out.print("\t"+sequenceNamesList.get(i));
		}
		out.println();
		int nBins = distProp.getDistribution().length;
		double min= 0;
		for(int i=0;i<nBins;i++) {
			out.print(""+format.format(min));
			for(Distribution d:distPropPerSequenceList) {
				out.print("\t"+format.format(d.getDistribution()[i]));
			}
			out.println();
			min+=0.01;
		}
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
		List<PileupAlleleCall> calls = pileup.getAlleleCalls(1);
		if(calls.size()<minRD) return;
		Map<String, Integer> alleleCounts = new TreeMap<String, Integer>();
		for(PileupAlleleCall call:calls) {
			String sequence = call.getSequence().toString();
			char qsC = call.getQualityScores().charAt(0);
			int qs = ((byte)qsC)-33;
			if(qs>=minBaseQualityScore) {
				Integer count = alleleCounts.get(sequence);
				if(count==null) {
					alleleCounts.put(sequence, 1);
				} else {
					alleleCounts.put(sequence, count+1);
				}
			}
			
		}
		if(alleleCounts.size()==0) return;
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
			if(currentSequencePropDist!=null) currentSequencePropDist.processDatapoint(prop);
		}
		coveredGenomeSize++;
		if(progressNotifier!=null && coveredGenomeSize%10000==0) {
			int progress = (int)(coveredGenomeSize/10000);
			generator.setKeepRunning(progressNotifier.keepRunning(progress));
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
	public void onSequenceStart(QualifiedSequence sequence) {
		//log.info("Starting sequence: "+sequenceName);
		if(sequence.getLength()>100000) {
			sequenceNamesList.add(sequence.getName());
			currentSequencePropDist = new Distribution(0, 0.5, 0.01);
			distPropPerSequenceList.add(currentSequencePropDist);
		} else {
			currentSequencePropDist = null;
		}
		
	}

	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {
		//log.info("Finished sequence: "+sequenceName);
		
	}
	

}
