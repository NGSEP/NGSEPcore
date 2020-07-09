package ngsep.vcf;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.transcriptome.Gene;
import ngsep.transcriptome.Transcript;
import ngsep.transcriptome.Transcriptome;
import ngsep.transcriptome.VariantFunctionalAnnotation;
import ngsep.transcriptome.VariantFunctionalAnnotationType;
import ngsep.transcriptome.io.GFF3TranscriptomeHandler;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.Sample;
import ngsep.variants.io.SimpleSamplesFileHandler;


public class VCFAlleleSharingStatisticsCalculator {
	
	// Constants for default values
	public static final int DEF_WINDOW_LENGTH = 100000;
	public static final int DEF_STEP_LENGTH = 10000;
	
	// Logging and progress
	private Logger log = Logger.getLogger(VCFAlleleSharingStatisticsCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile;
	private List<Sample> samples = new ArrayList<Sample>();
	private String outputFile;
	private String populationsGroup1;
	private String populationsGroup2;
	private String transcriptomeFile;
	private Transcriptome transcriptome; 
	private boolean includeIntrons = false;
	private int windowLength = DEF_WINDOW_LENGTH;
	private int stepLength = DEF_STEP_LENGTH;
	
	// Model attributes
	private Set<String> samplesG1;
	private Set<String> samplesG2;
	
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
	
	public List<Sample> getSamples() {
		return samples;
	}
	public void setSamples(List<Sample> samples) {
		this.samples = samples;
	}
	public void setSamples(String samplesFile) throws IOException {
		SimpleSamplesFileHandler handler = new SimpleSamplesFileHandler();
		samples = handler.loadSamples(samplesFile);
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public String getPopulationsGroup1() {
		return populationsGroup1;
	}
	public void setPopulationsGroup1(String populationsGroup1) {
		this.populationsGroup1 = populationsGroup1;
	}
	
	public String getPopulationsGroup2() {
		return populationsGroup2;
	}
	public void setPopulationsGroup2(String populationsGroup2) {
		this.populationsGroup2 = populationsGroup2;
	}
	
	public String getTranscriptomeFile() {
		return transcriptomeFile;
	}
	public void setTranscriptomeFile(String transcriptomeFile) {
		this.transcriptomeFile = transcriptomeFile;
	}
	public Transcriptome getTranscriptome() {
		return transcriptome;
	}

	public void setTranscriptome(Transcriptome transcriptome) {
		this.transcriptome = transcriptome;
	}
	
	public int getWindowLength() {
		return windowLength;
	}
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}
	public void setWindowLength(String value) {
		setWindowLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getStepLength() {
		return stepLength;
	}
	public void setStepLength(int stepLength) {
		this.stepLength = stepLength;
	}
	public void setStepLength(String value) {
		setStepLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	public boolean isIncludeIntrons() {
		return includeIntrons;
	}
	public void setIncludeIntrons(boolean includeIntrons) {
		this.includeIntrons = includeIntrons;
	}
	public void setIncludeIntrons(Boolean includeIntrons) {
		this.setIncludeIntrons(includeIntrons.booleanValue());
	}

	public static void main(String[] args) throws Exception {
		VCFAlleleSharingStatisticsCalculator instance = new VCFAlleleSharingStatisticsCalculator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
		
		
	}

	public void run() throws IOException {
		logParameters();
		samplesG1 = Sample.getSampleIds(samples,populationsGroup1.split(","));
		log.info("Samples group 1: "+samplesG1);
		samplesG2 = Sample.getSampleIds(samples,populationsGroup2.split(","));
		log.info("Samples group 2: "+samplesG1);
		if(transcriptomeFile!=null) {
			GFF3TranscriptomeHandler transcriptomeHandler = new GFF3TranscriptomeHandler();
			transcriptome = transcriptomeHandler.loadMap(transcriptomeFile);
			log.info("Loaded transcriptome");
		}
		Map<String,List<Double>> sharingStats;
		if(inputFile!=null) sharingStats = calculateSharingStatistics(inputFile);
		else sharingStats = calculateSharingStatistics(System.in);
		if(sharingStats==null) {
			log.info("Process cancelled");
			return;
		}
		if(outputFile==null) printSharingStats (sharingStats,System.out);
		else {
			try (PrintStream out = new PrintStream(outputFile)) {
				printSharingStats(sharingStats,out);
			}
		}
		log.info("Process finished");
	}
	
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		if(inputFile != null) out.println("Input VCF file: "+inputFile);
		else out.println("Read variants from standard input");
		if(outputFile != null) out.println("Output file: "+outputFile);
		else out.println("Write to standard output");
		out.println("Populations group 1: "+populationsGroup1);
		out.println("Populations group 2: "+populationsGroup2);
		if(transcriptomeFile!=null) {
			out.println("Transcriptome file to calculate statistics by gene: "+transcriptomeFile);
			if (includeIntrons) out.println("Include introns in the calculation");
		} else {
			out.println("Window length (bp): "+getWindowLength());
			out.println("Length of step to next window (bp): "+getStepLength());
		}
		log.info(""+os.toString());
		
	}
	public Map<String, List<Double>> calculateSharingStatistics(String vcfFile) throws IOException {
		try (VCFFileReader in = new VCFFileReader(vcfFile)){
			return calculateSharingStatistics(in);
		}
	}
	public Map<String, List<Double>> calculateSharingStatistics(InputStream is) throws IOException {
		try (VCFFileReader in = new VCFFileReader(is)) {
			return calculateSharingStatistics(in);
		}
	}
	public Map<String, List<Double>> calculateSharingStatistics(VCFFileReader in) throws IOException  {
		in.setLog(log);
		if(transcriptome!=null) return calculateSharingStatisticsByGene(in);
		else return calculateSharingStatisticsByWindow(in);
	}

	public Map<String, List<Double>> calculateSharingStatisticsByWindow(VCFFileReader in) throws IOException {
		Map<String, List<Double>> answer = new LinkedHashMap<String, List<Double>>();
		Iterator<VCFRecord> it = in.iterator();
		List<String> sampleIds = in.getSampleIds();
		String seqName = "";
		int curWindowStart = 0;
		Queue<VCFRecord> variantsInProgress = new LinkedList<VCFRecord>();
		int n=0;
		
		while(it.hasNext()) {
			VCFRecord record = it.next();
			GenomicVariant var = record.getVariant();
			boolean seqNameEq = seqName.equals(record.getVariant().getSequenceName()); 
			if(!seqNameEq || curWindowStart + windowLength <= var.getFirst()) {
				if(!seqNameEq) {
					while(variantsInProgress.size()>0) {
						answer.put(seqName+"\t"+curWindowStart,processWindow(variantsInProgress,sampleIds,curWindowStart));
						curWindowStart+=stepLength;
					}
					seqName = record.getVariant().getSequenceName();
					curWindowStart = 1;
					
				} else {
					while(curWindowStart + windowLength <= var.getFirst() ) {
						answer.put(seqName+"\t"+curWindowStart,processWindow(variantsInProgress,sampleIds,curWindowStart));
						curWindowStart+=stepLength;
					}
				}
			}
			variantsInProgress.add(record);
			n++;
			if (progressNotifier!=null && n%1000==0) {
				int progress = n/1000;
				if (!progressNotifier.keepRunning(progress)) {
					return null;
				}
			}
		}
		while(variantsInProgress.size()>0) {
			answer.put(seqName+"\t"+curWindowStart,processWindow(variantsInProgress,sampleIds,curWindowStart));
			curWindowStart+=stepLength;
		}
		return answer;
	}

	private List<Double> processWindow(Queue<VCFRecord> variantsInProgress, List<String> sampleIds, int windowStart) {
		//Select window vars
		List<VCFRecord> windowVars = new ArrayList<VCFRecord>();
		VCFRecord firstNextWindow = null;
		for(VCFRecord record:variantsInProgress) {
			int pos = record.getVariant().getFirst();
			if(windowStart<=pos && pos < windowStart + windowLength ) {
				windowVars.add(record);
			}
			if(firstNextWindow==null && windowStart + stepLength < pos) {
				firstNextWindow = record;
			}
		}
		List<Double> answer = calculateSharingStatisticsRegion(windowVars, sampleIds, windowLength);
		//Remove variants that will not be used in the next window
		if(firstNextWindow==null) variantsInProgress.clear();
		else {
			VCFRecord r = variantsInProgress.peek();
			while(r!=firstNextWindow) {
				variantsInProgress.poll();
				r = variantsInProgress.peek();
			}
		}
		return answer;
	}

	
	public Map<String, List<Double>> calculateSharingStatisticsByGene(VCFFileReader in) throws IOException {
		Map<String, List<Double>> answer = new TreeMap<String, List<Double>>();
		Iterator<VCFRecord> it = in.iterator();
		List<String> sampleIds = in.getSampleIds();
		String currentTranscriptId = "";
		String currentGeneId = "";
		int currentGeneLength = 0;
		List<VCFRecord> currentGeneVars = new ArrayList<VCFRecord>();
		int n = 0;
		while(it.hasNext()) {
			VCFRecord record = it.next();
			Transcript transcript = getTranscript(record);
			Gene g = null;
			if(transcript!=null) {
				g = transcript.getGene();
			}
			if(g==null || !currentGeneId.equals(g.getId())) {
				if(currentGeneId.length()>0) {
					List<Double> stats = calculateSharingStatisticsRegion(currentGeneVars, sampleIds, currentGeneLength);
					answer.put(currentTranscriptId+"\t"+currentGeneId,stats);
				}
				if(g!=null) {
					currentTranscriptId = transcript.getId();
					currentGeneId = g.getId();
					currentGeneLength = calculateLength (g);
					currentGeneVars.clear();
				}
				else {
					currentTranscriptId = "";
					currentGeneId = "";
					currentGeneLength = 0;
				}
			}
			if(currentGeneId.length()>0 && includeRecordByAnnotation(record)) {
				currentGeneVars.add(record);
			}
			n++;
			if (progressNotifier!=null && n%1000==0) {
				int progress = n/1000;
				if (!progressNotifier.keepRunning(progress)) {
					return null;
				}
			}
		}
		if(currentGeneId.length()>0) {
			List<Double> stats = calculateSharingStatisticsRegion(currentGeneVars, sampleIds, currentGeneLength);
			answer.put(currentTranscriptId+"\t"+currentGeneId,stats);
		}
		return answer;
	}
	private boolean includeRecordByAnnotation(VCFRecord record) {
		VariantFunctionalAnnotation ann = record.getNGSEPFunctionalAnnotation();
		if(ann==null) return false;
		if(includeIntrons && VariantFunctionalAnnotationType.ANNOTATION_INTRON.equals(ann.getTypeName())) return true;
		return ann.isCoding();
	}

	private Transcript getTranscript(VCFRecord record) {
		GenomicVariantAnnotation ann = record.getInfoField(GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ID);
		if(ann !=null) {
			String transcriptId = (String)ann.getValue();
			Transcript t = transcriptome.getTranscript(transcriptId);
			return t;
		}
		return null;
	}

	private int calculateLength(Gene g) {
		List<Transcript> transcripts = transcriptome.getTranscriptsByGene(g.getId());
		int maxLength = 0;
		for(Transcript t:transcripts) {
			int l = t.length();
			if(includeIntrons) l = t.getLast()-t.getFirst()+1;
			if(l>maxLength) maxLength = l;
		}
		return maxLength;
	}
	private List<Double> calculateSharingStatisticsRegion(List<VCFRecord> regionVars, List<String> sampleIds, int regionLength) {
		int nSamples = sampleIds.size();
		int [][] sampleDifferences = new int [nSamples][nSamples];
		resetDifferences(sampleDifferences);
		int [] segregatingSites = new int [3];
		Arrays.fill(segregatingSites, 0);
		int nVars = 0;
		for(VCFRecord record: regionVars) {
			updateCounts(record,sampleIds,sampleDifferences,segregatingSites);
			nVars++;
		}
		return calculateStatistics(sampleIds, sampleDifferences, nVars, segregatingSites, regionLength);
	}

	public void resetDifferences(int[][] sampleDifferences) {
		for(int i=0;i<sampleDifferences.length;i++) {
			for(int j=0;j<sampleDifferences[0].length;j++) {
				sampleDifferences[i][j] = 0;
			}
		}
	}
	private void updateCounts(VCFRecord record, List<String> sampleIds, int[][] sampleDifferences, int [] segregatingSites ) {
		List<CalledGenomicVariant> calls = record.getCalls();
		boolean variableG1 = false;
		boolean variableG2 = false;
		boolean variableG1G2 = false;
		for(int i=0;i<calls.size();i++) {
			CalledGenomicVariant v1 = calls.get(i);
			String [] alleles1 = v1.getCalledAlleles();
			String sampleId1 = sampleIds.get(i);
			boolean s1G1 = samplesG1.contains(sampleId1);
			boolean s1G2 = samplesG2.contains(sampleId1);
			for(int j=i+1;j<calls.size();j++) {
				CalledGenomicVariant v2 = calls.get(j);
				String [] alleles2 = v2.getCalledAlleles();
				String sampleId2 = sampleIds.get(j);
				boolean s2G1 = samplesG1.contains(sampleId2);
				boolean s2G2 = samplesG2.contains(sampleId2);
				int diffs = countDifferences (alleles1,alleles2);
				sampleDifferences[i][j]+=diffs;
				sampleDifferences[j][i]+=diffs;
				if(diffs>0) {
					if(s1G1 && s2G1) variableG1 = true;
					if(s1G2 && s2G2) variableG2 = true;
					if((s1G1 && s2G2) || (s1G2 && s2G1)) variableG1G2 = true;
				}
			}
		}
		if(variableG1) segregatingSites[0]++;
		if(variableG2) segregatingSites[1]++;
		if(variableG1G2) segregatingSites[2]++;
		
	}
	private int countDifferences(String[] alleles1, String[] alleles2) {
		if(alleles1.length == 0 || alleles2.length == 0 ) return 0;
		//Assumes that samples are  homozygous
		if(alleles1.length!=alleles2.length) return 0;
		if(alleles1.length==1 && alleles2.length==1 && !alleles1[0].equals(alleles2[0])) return 1;
		return 0;
	}
	private List<Double> calculateStatistics(List<String> sampleIds,int[][] sampleDifferences, int nVars, int [] segregatingSites, int regionLength) {
		double avgDiffsG1=0;
		double avgDiffsG2=0;
		double avgDiffsBetween=0;
		double avgDiffsG1G2=0;
		double avgDiffsTotal=0;
		int harmonicValue = 0;
		int harmonicG1 = 0;
		int harmonicG2 = 0;
		int numPairsG1 = 0;
		int numPairsG2 = 0;
		int numPairsBetween = 0;
		
		for(int i=0;i<sampleIds.size();i++) {
			String sampleId1 = sampleIds.get(i);
			boolean s1G1 = samplesG1.contains(sampleId1);
			boolean s1G2 = samplesG2.contains(sampleId1);
			boolean s1G = s1G1 || s1G2;
			for(int j=i+1;j<sampleIds.size();j++) {
				String sampleId2 = sampleIds.get(j);
				boolean s2G1 = samplesG1.contains(sampleId2);
				boolean s2G2 = samplesG2.contains(sampleId2);
				boolean s2G = s2G1 || s2G2;
				int diffsPair = sampleDifferences[i][j];
				avgDiffsTotal+=diffsPair;
				if(s1G && s2G) {
					avgDiffsG1G2+=diffsPair;
					if(s1G1 && s2G1) {
						avgDiffsG1 +=diffsPair;
						numPairsG1++;
					}
					if(s1G2 && s2G2) {
						avgDiffsG2 +=diffsPair;
						numPairsG2++;
					}
					if((s1G1 && s2G2) || (s1G2 && s2G1)  ) {
						avgDiffsBetween += diffsPair;
						numPairsBetween++;
					}
				}
			}
			if(i==samplesG1.size()-1) harmonicG1 = harmonicValue;
			if(i==samplesG2.size()-1) harmonicG2 = harmonicValue;
			harmonicValue+=(1/(i+1));
			
		}
		double avgDiffsWithin = avgDiffsG1+avgDiffsG2;
		List<Double> answer = new ArrayList<Double>();
		double regionLengthKBP = (double)regionLength/1000.0;
		answer.add(regionLengthKBP);
		answer.add((double) nVars);
		for(int i=0;i<segregatingSites.length;i++) {
			answer.add((double) segregatingSites[i]);
		}
		avgDiffsG1/=numPairsG1;
		avgDiffsG1/=regionLengthKBP;
		answer.add(avgDiffsG1);
		//answer.add((double) numPairsG1);
		//double dG1 = avgDiffsG1-((double)segregatingSites[0])/harmonicG1;
		//answer.add(dG1);
		avgDiffsG2/=numPairsG2;
		avgDiffsG2/=regionLengthKBP;
		answer.add(avgDiffsG2);
		//answer.add((double) numPairsG2);
		//double dG2 = avgDiffsG2-((double)segregatingSites[1])/harmonicG2;
		//answer.add(dG2);
		avgDiffsBetween/=numPairsBetween;
		avgDiffsBetween/=regionLengthKBP;
		answer.add(avgDiffsBetween);
		//answer.add((double) numPairsBetween);
		avgDiffsWithin/=(numPairsG1+numPairsG2);
		avgDiffsWithin/=regionLengthKBP;
		answer.add(avgDiffsWithin);
		int n12 = samplesG1.size()+samplesG2.size();
		int numPairsG1G2 = n12*(n12-1)/2;
		avgDiffsG1G2/=numPairsG1G2;
		avgDiffsG1G2/=regionLengthKBP;
		answer.add(avgDiffsG1G2);
		int numPairsTotal = sampleIds.size()*(sampleIds.size()-1)/2;
		avgDiffsTotal/=numPairsTotal;
		avgDiffsTotal/=regionLengthKBP;
		answer.add(avgDiffsTotal);
		double fstG1G2 = 0;
		if(avgDiffsBetween>0) fstG1G2 = (avgDiffsBetween-avgDiffsWithin)/avgDiffsBetween;
		answer.add(fstG1G2);
		double dG1 = regionLengthKBP*avgDiffsG1-((double)segregatingSites[0])/harmonicG1;
		answer.add(dG1);
		double dG2 = regionLengthKBP*avgDiffsG2-((double)segregatingSites[1])/harmonicG2;
		answer.add(dG2);
		
		return answer;
	}
	public void printSharingStats(Map<String, List<Double>> sharingStats, PrintStream out) {
		out.println("#Samples group 1: "+samplesG1);
		out.println("#Samples group 2: "+samplesG2);
		DecimalFormat fmt = new DecimalFormat("###0.00");
		for(String key:sharingStats.keySet()) {
			out.print(key);
			List<Double> stats = sharingStats.get(key);
			for(double d:stats) {
				out.print("\t"+fmt.format(d));
			}
			out.println();
		}
		
	}
}
