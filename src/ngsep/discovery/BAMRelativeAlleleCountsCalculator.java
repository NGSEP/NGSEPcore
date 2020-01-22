package ngsep.discovery;
import java.io.ByteArrayOutputStream;
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
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.math.Distribution;
import ngsep.sequences.QualifiedSequence;




public class BAMRelativeAlleleCountsCalculator implements PileupListener {

	// Constants for default values
	public static final int DEF_MIN_RD = 10;
	public static final int DEF_MAX_RD = 1000;
	public static final int DEF_MIN_BASE_QUALITY_SCORE = 20;
	
	
	// Logging and progress
	private Logger log = Logger.getLogger(BAMRelativeAlleleCountsCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	private long coveredGenomeSize = 0;
	
	// Parameters
	private String inputFile = null;
	private String outputFile;
	private int minRD = DEF_MIN_RD;
	private int maxRD = DEF_MAX_RD;
	private int minBaseQualityScore = DEF_MIN_BASE_QUALITY_SCORE;
	private GenomicRegionSortedCollection<GenomicRegion> regionsToFilter = null;
	private GenomicRegionSortedCollection<GenomicRegion> regionsToSelect = null;
	private boolean secondaryAlns = false;
	
	// Model attributes
	private AlignmentsPileupGenerator generator;
	
	
	private Distribution distProp = new Distribution(0, 0.5, 0.01);
	private List<String> sequenceNamesList = new ArrayList<>();
	private List<Distribution> distPropPerSequenceList = new ArrayList<>();
	private Distribution currentSequencePropDist=null;
	private Distribution distNumAlleles = new Distribution(1, 10, 1);
	
	private static DecimalFormat format = new DecimalFormat("##0.0#");
	
	// Get and set methods
	
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
	
	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public int getMinRD() {
		return minRD;
	}
	public void setMinRD(int minRD) {
		this.minRD = minRD;
	}
	public void setMinRD(String value) {
		this.setMinRD((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	public int getMaxRD() {
		return maxRD;
	}
	public void setMaxRD(int maxRD) {
		this.maxRD = maxRD;
	}
	public void setMaxRD(String value) {
		this.setMaxRD((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getMinBaseQualityScore() {
		return minBaseQualityScore;
	}
	public void setMinBaseQualityScore(int minBaseQualityScore) {
		this.minBaseQualityScore = minBaseQualityScore;
	}
	public void setMinBaseQualityScore(String value) {
		this.setMinBaseQualityScore((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public GenomicRegionSortedCollection<GenomicRegion> getRegionsToFilter() {
		return regionsToFilter;
	}
	public void setRegionsToFilter(GenomicRegionSortedCollection<GenomicRegion> regionsToFilter) {
		this.regionsToFilter = regionsToFilter;
	}
	public void setRegionsToFilter(String regionsFile) throws IOException {
		if(regionsFile==null || regionsFile.length()==0) {
			this.regionsToFilter = null;
			return;
		}
		SimpleGenomicRegionFileHandler regionFileHandler = new SimpleGenomicRegionFileHandler();
		List<GenomicRegion> regions = regionFileHandler.loadRegions(regionsFile);
		this.regionsToFilter = new GenomicRegionSortedCollection<GenomicRegion>(regions);
	}
	
	public GenomicRegionSortedCollection<GenomicRegion> getRegionsToSelect() {
		return regionsToSelect;
	}
	public void setRegionsToSelect(GenomicRegionSortedCollection<GenomicRegion> regionsToSelect) {
		this.regionsToSelect = regionsToSelect;
	}
	public void setRegionsToSelect(String regionsFile) throws IOException {
		if(regionsFile==null || regionsFile.length()==0) {
			this.regionsToSelect = null;
			return;
		}
		SimpleGenomicRegionFileHandler regionFileHandler = new SimpleGenomicRegionFileHandler();
		List<GenomicRegion> regions = regionFileHandler.loadRegions(regionsFile);
		this.regionsToSelect = new GenomicRegionSortedCollection<GenomicRegion>(regions);
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
	
	public static void main(String[] args) throws Exception {
		BAMRelativeAlleleCountsCalculator instance = new BAMRelativeAlleleCountsCalculator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run () throws IOException {
		logParameters();
		if (inputFile == null) throw new IOException("The alignments input file is a required parameter");
		runProcess(inputFile);
		if (outputFile==null ) printResults(System.out);
		else {
			try (PrintStream out = new PrintStream(outputFile)) {
				printResults(out);
			}
		}
		log.info("Process finished");
	}
	
	private void logParameters () {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		if(inputFile != null) out.println("Input file: "+inputFile);
		if(outputFile != null) out.println("Output file: "+outputFile);
		else out.println("Write to standard output");
		out.println("Minimum read depth: "+getMinRD());
		out.println("Maximum read depth: "+getMaxRD());
		out.println("Minimum base quality score: "+getMinBaseQualityScore());
		if (regionsToFilter!=null) out.println("Loaded "+regionsToFilter.size()+" regions to filter");
		if (regionsToSelect!=null) out.println("Loaded "+regionsToSelect.size()+" regions to select");
		if (secondaryAlns) out.println("Include secondary alignments");
		log.info(""+os.toString());
	}
	
	public void runProcess (String filename) throws IOException {
		coveredGenomeSize = 0;
		generator = new AlignmentsPileupGenerator();
		generator.addListener(this);
		generator.setLog(log);
		generator.setProcessSecondaryAlignments(secondaryAlns);
		generator.setMaxAlnsPerStartPos(maxRD);
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
		if(regionsToFilter!=null) {
			GenomicRegionSortedCollection<GenomicRegion> spanningFilterRegions = regionsToFilter.findSpanningRegions(pileup.getSequenceName(), pileup.getPosition());
			if(spanningFilterRegions.size()>0) return;
		}
		if(regionsToSelect!=null) {
			GenomicRegionSortedCollection<GenomicRegion> spanningRegions = regionsToSelect.findSpanningRegions(pileup.getSequenceName(), pileup.getPosition());
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
