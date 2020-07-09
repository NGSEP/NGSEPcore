/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.vcf;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.math.NumberArrays;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;
import ngsep.variants.io.SimpleSamplesFileHandler;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class VCFWindowIntrogressionAnalysis {
	
	// Constants for default values
	public static final double DEF_MIN_PCT_GENOTYPED = 80;
	public static final double DEF_MIN_DIFF_AF = 0.6;
	public static final double DEF_MAX_MAF_WITHIN = 0.4;
	public static final int DEF_WINDOW_SIZE = 50;
	public static final int DEF_OVERLAP = 0;
	public static final int DEF_MATCH_SCORE = 1;
	public static final int DEF_MISMATCH_SCORE = -1;
	public static final int DEF_MIN_SCORE = 30;
	public static final String UNASSIGNED_GROUP_ID="U";
	
	// Logging and progress
	private Logger log = Logger.getLogger(VCFWindowIntrogressionAnalysis.class.getName());
	private ProgressNotifier progressNotifier=null;

	// Parameters
	private String inputFile;
	private String populationsFile;
	private String outputPrefix;
	private double minPCTGenotyped = DEF_MIN_PCT_GENOTYPED;
	private double minDiffAF = DEF_MIN_DIFF_AF;
	private double maxMAFWithin = DEF_MAX_MAF_WITHIN;
	private int windowSize = DEF_WINDOW_SIZE;
	private int overlap = DEF_OVERLAP;
	private int matchScore = DEF_MATCH_SCORE;
	private int mismatchScore = DEF_MISMATCH_SCORE;
	private int minScore = DEF_MIN_SCORE;
	private boolean printVCF = false;
	private boolean printUnassigned = false;
	
	// Model attributes
	private Queue<WindowScores> currentWindows = new LinkedList<WindowScores>();
	private int [][] assignmentStats;
	private Introgression [] currentIntrogressions;
	private List<Introgression> sequenceIntrogressions;
	
	// Get and set methds
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
	
	public String getPopulationsFile() {
		return populationsFile;
	}
	public void setPopulationsFile(String populationsFile) {
		this.populationsFile = populationsFile;
	}
	
	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	
	public double getMinPCTGenotyped() {
		return minPCTGenotyped;
	}
	public void setMinPCTGenotyped(double minPctGenotyped) {
		this.minPCTGenotyped = minPctGenotyped;
	}
	public void setMinPCTGenotyped(String value) {
		setMinPCTGenotyped((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	public double getMinDiffAF() {
		return minDiffAF;
	}
	public void setMinDiffAF(double minDiffAF) {
		this.minDiffAF = minDiffAF;
	}
	public void setMinDiffAF(String value) {
		this.setMinDiffAF((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	public double getMaxMAFWithin() {
		return maxMAFWithin;
	}
	public void setMaxMAFWithin(double maxMAFWithinPops) {
		this.maxMAFWithin = maxMAFWithinPops;
	}
	public void setMaxMAFWithin(String value) {
		this.setMaxMAFWithin((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	public int getWindowSize() {
		return windowSize;
	}
	public void setWindowSize(int windowSize) {
		this.windowSize = windowSize;
	}
	public void setWindowSize(String value) {
		this.setWindowSize((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public int getOverlap() {
		return overlap;
	}
	public void setOverlap(int overlap) {
		this.overlap = overlap;
	}
	public void setOverlap(String value) {
		this.setOverlap((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public int getMatchScore() {
		return matchScore;
	}
	public void setMatchScore(int matchScore) {
		this.matchScore = matchScore;
	}
	public void setMatchScore(String value) {
		this.setMatchScore((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public int getMismatchScore() {
		return mismatchScore;
	}
	public void setMismatchScore(int mismatchScore) {
		this.mismatchScore = mismatchScore;
	}
	public void setMismatchScore(String value) {
		this.setMismatchScore((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public int getMinScore() {
		return minScore;
	}
	public void setMinScore(int minScore) {
		this.minScore = minScore;
	}
	public void setMinScore(String value) {
		this.setMinScore((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public boolean isPrintVCF() {
		return printVCF;
	}
	public void setPrintVCF(boolean printVCF) {
		this.printVCF = printVCF;
	}
	public void setPrintVCF(Boolean printVCF) {
		this.setPrintVCF(printVCF.booleanValue());
	}
	
	public boolean isPrintUnassigned() {
		return printUnassigned;
	}
	public void setPrintUnassigned(boolean printUnassigned) {
		this.printUnassigned = printUnassigned;
	}
	public void setPrintUnassigned(Boolean printUnassigned) {
		this.setPrintUnassigned(printUnassigned.booleanValue());
	}
	
	public static void main(String[] args) throws Exception {
		VCFWindowIntrogressionAnalysis instance = new VCFWindowIntrogressionAnalysis();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
		
	}
	
	public void run () throws IOException {
		logParameters();
		if(outputPrefix == null) throw new IOException("The prefix of the output files is a required parameter");
		if(populationsFile==null) throw new IOException("The fie with the description of the populations is a required parameter");
		if(inputFile == null) {
			runIntrogressions(System.in);
		} else {
			runIntrogressions(inputFile);
		}
		log.info("Process finished");
	}
	
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		if(inputFile != null) out.println("VCF input file: "+inputFile);
		else out.println("Read variants from standard input");
		out.println("Populations information file: "+getPopulationsFile());
		out.println("Prefix for output files: "+outputPrefix);
		out.println("Minimim percentage of variants genotyped: "+getMinPCTGenotyped());
		out.println("Minimum difference in allele frequencies to call a variant segregating: "+getMinDiffAF());
		out.println("Maximum MAF within a group to call a variant segregating: "+getMaxMAFWithin());
		out.println("Number of variants per window: "+getWindowSize());
		out.println("Overlap between windows: "+getOverlap());
		out.println("Score for a genotype match between a sample haplotype and a population major allele: "+getMatchScore());
		out.println("Score for a genotype mismatch between a sample haplotype and a population major allele: "+getMismatchScore());
		out.println("Minimum score to assign a sample haplotype to a population: "+getMinScore());
		if (printVCF) out.println("Output a VCF file with the variants showing segregation between at least two populations");
		if (printUnassigned) out.println("Report introgression events for unassigned haplotypes");
		log.info(""+os.toString());
	}
	
	public void runIntrogressions(InputStream is) throws IOException{
		try (VCFFileReader reader = new VCFFileReader(is)){
			runIntrogressions(reader);
		}
	}
	public void runIntrogressions(String vcfFile) throws IOException {
		
		try (VCFFileReader reader = new VCFFileReader(vcfFile)){
			runIntrogressions(reader);
		}
	}
	public void runIntrogressions(VCFFileReader reader) throws IOException {
		PrintStream outVCF = null;
		VCFFileWriter writer = new VCFFileWriter();
		try (PrintStream outAssignments = new PrintStream(outputPrefix+"_assignments.txt");
			 PrintStream outIntrogressions = new PrintStream(outputPrefix+"_introgressions.txt");
			 PrintStream outStatistics = new PrintStream(outputPrefix+"_assignmentStats.txt");) {
			
			
			VCFFileHeader header = reader.getHeader();
			if(printVCF) {
				outVCF = new PrintStream(outputPrefix+"_segregating.vcf");
				writer.printHeader(header, outVCF);
				reader.setLoadMode(VCFFileReader.LOAD_MODE_QUALITY);
			} else {
				reader.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			}
			List<Sample> samples = header.getSamples();
			printHeaderAssignments(samples,outAssignments);
			int nSamples = samples.size();
			currentIntrogressions = new Introgression[nSamples];
			Arrays.fill(currentIntrogressions, null);
			sequenceIntrogressions = new ArrayList<Introgression>();
			String currentSeqName = null;
			SimpleSamplesFileHandler handler = new SimpleSamplesFileHandler();
			Set<String> groupIds = handler.fillGroups (populationsFile, samples);
			List<String> groupIdsList = new ArrayList<String>(groupIds);
			Map<String,Integer> groupIdsReverseIdxMap = new HashMap<>();
			for(int i=0;i<groupIdsList.size();i++) {
				groupIdsReverseIdxMap.put(groupIdsList.get(i), i);
			}
			printHeaderIntrogressions(groupIdsList,outIntrogressions);
			int nG = groupIds.size();
			assignmentStats = new int [nSamples][nG+3];
			NumberArrays.initializeIntMatrix(assignmentStats);
			
			int [][] numDiscriminative = new int[nG][nG];
			NumberArrays.initializeIntMatrix(numDiscriminative);
			int [] numHeterozygous = new int [nG];
			Arrays.fill(numHeterozygous, 0);
			int [] numUndecided = new int [nG];
			Arrays.fill(numUndecided, 0);
			int n=0;
			int numVarsLastWindow=0;
			Iterator<VCFRecord> it = reader.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				boolean seqChange = !record.getSequenceName().equals(currentSeqName);
				if(seqChange) {
					processWindows(samples, groupIdsList, true, outAssignments);
					processSequenceIntrogressions(outIntrogressions,groupIdsReverseIdxMap);
					currentSeqName = record.getSequenceName();
					log.info("Starting sequence: "+currentSeqName);
					numVarsLastWindow = 0;
					currentWindows.add(new WindowScores(currentSeqName, nSamples, nG));
				}
				GenomicVariant var = record.getVariant();
				if(!var.isBiallelic()) continue;
				List<CalledGenomicVariant> varCalls = record.getCalls();
				
				double [] refAlleleFreqs = calculateFrequencies(varCalls, groupIdsList, samples);
				boolean selectVariant = false;
				//System.out.print(""+var.getFirst());
				for(int i = 0; i<refAlleleFreqs.length;i++) {
					//System.out.print(" "+refAlleleFreqs[i]);
					if(refAlleleFreqs[i]<0) continue;
					for(int j=i+1;j<refAlleleFreqs.length;j++) {
						if(refAlleleFreqs[j]<0) continue;
						double diff = Math.abs(refAlleleFreqs[i]-refAlleleFreqs[j]);
						if(diff>minDiffAF) {
							selectVariant = true;
							(numDiscriminative[i][j])++;
							(numDiscriminative[j][i])++;
						}
					}
				}
				//System.out.println();
				if(!selectVariant) continue;
				
				List<CalledGenomicVariant> groupCalls = calculateGroupCalls(refAlleleFreqs,var);
				for(int i=0;i<groupCalls.size();i++) {
					CalledGenomicVariant cv = groupCalls.get(i);
					if(cv.isUndecided()) numUndecided[i]++;
					else if (cv.isHeterozygous()) numHeterozygous[i]++;
				}
				for(WindowScores window:currentWindows) {
					window.addVariant(var, varCalls, groupCalls, matchScore, mismatchScore);
				}
				numVarsLastWindow++;
				processWindows(samples, groupIdsList, false, outAssignments);
				if(numVarsLastWindow+overlap==windowSize) {
					log.finest("Starting window after "+numVarsLastWindow+" biallelic discriminative variants in the last window. Last variant pos: "+var.getFirst()+". Number of windows: "+currentWindows.size());
					currentWindows.add(new WindowScores(var.getSequenceName(), nSamples, nG));
					numVarsLastWindow = 0;
				}
				if(outVCF!=null) writer.printVCFRecord(record, outVCF);
				n++;
				if (progressNotifier!=null && n%1000==0) {
					int progress = n/1000;
					if (!progressNotifier.keepRunning(progress)) {
						return;
					}
				}
			}
			processWindows(samples, groupIdsList, true, outAssignments);
			processSequenceIntrogressions(outIntrogressions, groupIdsReverseIdxMap);
			printReport(samples,groupIdsList,numDiscriminative,numHeterozygous,numUndecided,outStatistics);
		} finally {
			if(outVCF!=null) outVCF.close();
		}
	}
	private void printHeaderAssignments(List<Sample> samples, PrintStream outAssignments) {
		outAssignments.print("Chr\tFirst\tLast");
		for(Sample s:samples) {
			outAssignments.print("\t"+s.getId());
		}
		outAssignments.println();
	}
	private void printHeaderIntrogressions(List<String> groupIdsList, PrintStream outIntrogressions) {
		outIntrogressions.print("Chr\tFirst\tLast\tSample\tSample group\tHaplotype group\tTotal variants\tGenotyped variants");
		for(String groupId:groupIdsList) {
			outIntrogressions.print("\tScore "+groupId);
		}
		outIntrogressions.print("\tScore self group");
		outIntrogressions.println();
		
		
	}
	private void processWindows(List<Sample> samples, List<String> groupIdsList, boolean seqEnd ,PrintStream outAssignments) {
		while(currentWindows.size()>0) {
			WindowScores window = currentWindows.peek();
			if(window.getNumVariants()==windowSize) {
				processWindow(window,samples, groupIdsList, outAssignments);
			} else if (!seqEnd) {
				return;
			}
			currentWindows.remove();
		}
		if(seqEnd) {
			for(int i=0;i<currentIntrogressions.length;i++) {
				if(currentIntrogressions[i]!=null) {
					sequenceIntrogressions.add(currentIntrogressions[i]);
					currentIntrogressions[i] = null;
				}
			}
		}
	}
	
	private double[] calculateFrequencies(List<CalledGenomicVariant> varCalls,List<String> groupIds, List<Sample> samples) {
		double [] countsRefAllele = new double[groupIds.size()];
		Arrays.fill(countsRefAllele, 0);
		double [] countsAltAllele = new double[groupIds.size()];
		Arrays.fill(countsAltAllele, 0);
		double [] countsUndecided = new double[groupIds.size()];
		Arrays.fill(countsUndecided, 0);
		for(int i=0;i<varCalls.size();i++) {
			CalledGenomicVariant call = varCalls.get(i);
			String groupId = samples.get(i).getGroup();
			if(groupId==null) continue;
			int j = groupIds.indexOf(groupId);
			if(j<0) continue;
			if(call.isUndecided()) {
				(countsUndecided[j])+=2;
			} else if (call.isHeterozygous()) {
				(countsRefAllele[j])++;
				(countsAltAllele[j])++;
			} else if (call.isHomozygousReference()) {
				(countsRefAllele[j])+=2;
			} else {
				(countsAltAllele[j])+=2;
			}
		}
		double [] answer = new double[groupIds.size()];
		for(int i=0;i<answer.length;i++) {
			double totalCount = countsRefAllele[i]+countsAltAllele[i]+countsUndecided[i];
			double genotypedCount = countsRefAllele[i]+countsAltAllele[i];
			double pctGenoyped = 100.0*genotypedCount/totalCount;
			if(pctGenoyped<minPCTGenotyped) {
				answer[i]=-1;
			} else {
				answer[i] = countsRefAllele[i]/genotypedCount;
			}
			
		}
		return answer;
	}
	private List<CalledGenomicVariant> calculateGroupCalls(double[] refAlleleFreqs, GenomicVariant var) {
		List<CalledGenomicVariant> answer = new ArrayList<CalledGenomicVariant>();
		for(int i=0;i<refAlleleFreqs.length;i++) {
			double refAF = refAlleleFreqs[i];
			double maf = refAF;
			if(maf>0.5) maf = 1-maf;
			if(refAF<0) {
				answer.add(new CalledGenomicVariantImpl(var, new byte[0]));
				continue;
			} 
			byte [] genotype;
			if(maf>maxMAFWithin) {
				genotype = new byte[2];
				genotype[0]=0;
				genotype[1]=1;
				
			} else {
				genotype = new byte[1];
				genotype[0]=0;
				if(refAF<0.5) genotype[0]=1;
			}
			CalledGenomicVariantImpl call = new CalledGenomicVariantImpl(var, genotype);
			call.setGenotypeQuality((short)40);
			answer.add(call);
		}
		return answer;
	}
	private void processWindow(WindowScores window, List<Sample> samples, List<String> groupIds, PrintStream outMap) {
		outMap.print(window.getSequenceName()+"\t"+window.getFirst()+"\t"+window.getLast());
		int [] genotyped = window.getCountsGenotyped();
		int [] [] scores = window.getScores();
		int nVars = window.getNumVariants();
		for(int i=0;i<scores.length;i++) {
			double g = genotyped[i];
			if(g*100.0/nVars<minPCTGenotyped) {
				outMap.print("\tM");
				(assignmentStats[i][groupIds.size()])++;
				processCurrentIntrogressions(i, samples.get(i), window, null);
				continue;
			}
			int jMax = NumberArrays.getIndexMaximum(scores[i], -1);
			int bestScore = scores[i][jMax];
			String bestGroup = groupIds.get(jMax);
			if(bestScore<minScore) {
				outMap.print("\t"+UNASSIGNED_GROUP_ID);
				(assignmentStats[i][groupIds.size()+1])++;
				processCurrentIntrogressions(i, samples.get(i), window, UNASSIGNED_GROUP_ID);
				continue;
			}
			int jSecondMax = NumberArrays.getIndexMaximum(scores[i], jMax);
			int secondBest = scores[i][jSecondMax];
			if(bestScore-secondBest<10) {
				outMap.print("\t"+bestGroup+":"+bestScore+"/"+groupIds.get(jSecondMax)+":"+secondBest);
				(assignmentStats[i][groupIds.size()+2])++;
			} else {
				outMap.print("\t"+bestGroup);
				(assignmentStats[i][jMax])++;
				processCurrentIntrogressions(i,samples.get(i),window, bestGroup);
			}
		}
		outMap.println();
	}
	private void processCurrentIntrogressions(int sampleIdx, Sample sample, WindowScores window, String bestGroup) {
		if(sample.getGroup()==null) return;
		Introgression introgression = currentIntrogressions[sampleIdx];
		if(introgression!=null) {
			if(!introgression.getForeignGroup().equals(bestGroup)) {
				sequenceIntrogressions.add(introgression);
				currentIntrogressions[sampleIdx] = introgression = null;
			} else {
				introgression.addWindow(window,sampleIdx);
			}
		}
		
		if(introgression == null) {
			if(bestGroup!=null && !sample.getGroup().equals(bestGroup)) {
				currentIntrogressions[sampleIdx] = new Introgression(sample, window, sampleIdx, bestGroup);
			}
		} 
	}
	private void processSequenceIntrogressions(PrintStream outIntrogressions, Map<String,Integer> groupIdsReverseIdxMap) {
		Collections.sort(sequenceIntrogressions,GenomicRegionPositionComparator.getInstance());
		for(Introgression i:sequenceIntrogressions) {
			if(!printUnassigned && UNASSIGNED_GROUP_ID.equals(i.getForeignGroup())) continue;
			int sampleGroupIdx = groupIdsReverseIdxMap.get(i.getSample().getGroup());
			printIntrogression(i, sampleGroupIdx, outIntrogressions);
		}
		sequenceIntrogressions.clear();
	}
	private void printIntrogression(Introgression i, int sampleGroupIdx, PrintStream out) {
		Sample s = i.getSample();
		out.print(""+i.getSequenceName()+"\t"+i.getFirst()+"\t"+i.getLast());
		out.print("\t"+s.getId()+"\t"+s.getGroup()+"\t"+i.getForeignGroup()+"\t"+i.getNumVariants());
		out.print("\t"+i.getNumGenotyped());
		int [] scores = i.getScoresGroups();
		for(int j=0;j<scores.length;j++) {
			out.print("\t"+scores[j]);
		}
		out.print("\t"+scores[sampleGroupIdx]);
		out.println();
	}
	private void printReport(List<Sample> samples, List<String> groupIds, int[][] numDiscriminative, int[] numHeterozygous, int[] numUndecided, PrintStream out) {
		
		out.println("Variants segregating between each pair of groups");
		out.print("Group");
		for(int i=0;i<groupIds.size();i++) out.print("\t"+groupIds.get(i));
		out.println();
		for(int i=0;i<groupIds.size();i++) {
			out.print(groupIds.get(i));
			for(int j=0;j<groupIds.size();j++) {
				 out.print("\t"+numDiscriminative[i][j]);
			}
			out.println();
		}
		out.println();
		
		out.println("Variants with not enough individuals genotyped or with high MAF within each group");
		out.print("Group");
		for(int i=0;i<groupIds.size();i++) out.print("\t"+groupIds.get(i));
		out.println();
		out.print("Not genotyped");
		for(int i=0;i<groupIds.size();i++) out.print("\t"+numUndecided[i]);
		out.println();
		out.print("Heterozygous");
		for(int i=0;i<groupIds.size();i++) out.print("\t"+numHeterozygous[i]);
		out.println();
		out.println();
		out.println("Regions assigned to each group for each sample");
		out.print("Sample\tSample group");
		for(int i=0;i<groupIds.size();i++) out.print("\t"+groupIds.get(i));
		out.println("\tNon genotyped\tUnassigned\tMore than one group");
		for(int i=0;i<assignmentStats.length;i++) {
			Sample s = samples.get(i); 
			out.print(s.getId()+"\t");
			String sG = s.getGroup();
			out.print((sG!=null)?sG:"None");
			for(int j=0;j<assignmentStats[i].length;j++) {
				out.print("\t"+assignmentStats[i][j]);
			}
			out.println();
		}
	}
}
class WindowScores implements GenomicRegion {
	private int [] countsGenotyped;
	private int [][] scores;
	private int numVariants = 0;
	private String sequenceName;
	private int first=-1;
	private int last=-1;
	public WindowScores (String sequenceName, int numSamples, int numGroups) {
		this.sequenceName = sequenceName;
		countsGenotyped = new int [numSamples];
		Arrays.fill(countsGenotyped, 0);
		this.scores = new int [numSamples][numGroups];
		NumberArrays.initializeIntMatrix(this.scores);
	}
	public void addVariant(GenomicVariant var, List<CalledGenomicVariant> callsSamples, List<CalledGenomicVariant> callsGroups, int matchScore, int mismatchScore) {
		for(int i=0;i<countsGenotyped.length;i++) {
			if(!callsSamples.get(i).isUndecided()) countsGenotyped[i]++;
		}
		for(int i=0;i<scores.length;i++) {
			for(int j=0;j<scores[i].length;j++) {
				scores[i][j] += getScore(callsSamples.get(i),callsGroups.get(j),matchScore,mismatchScore);
			}
		}
		if(first == -1 || first > var.getFirst() ) first = var.getFirst();
		if(last == -1 || last < var.getLast() ) last = var.getLast();
		numVariants++;
	}
	private int getScore(CalledGenomicVariant call1, CalledGenomicVariant call2, int matchScore, int mismatchScore) {
		if(call1.isUndecided() || call2.isUndecided()) return 0;
		byte [] a1 = call1.getIndexesCalledAlleles();
		byte [] a2 = call2.getIndexesCalledAlleles();
		if(a1.length>1 || a2.length>1) return 0;
		if(a1[0]==a2[0]) return matchScore;
		return mismatchScore;
	}
	public int getNumVariants() {
		return numVariants;
	}
	
	public int[] getCountsGenotyped() {
		return countsGenotyped;
	}
	public int[][] getScores() {
		return scores;
	}
	
	public int getNumGroups() {
		return scores[0].length;
	}
	@Override
	public String getSequenceName() {
		return sequenceName;
	}
	@Override
	public int getFirst() {
		return first;
	}
	@Override
	public int getLast() {
		return last;
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return true;
	}
	@Override
	public boolean isNegativeStrand() {
		return false;
	}
}
class Introgression extends GenomicRegionImpl {
	private Sample sample;
	private String foreignGroup;
	private int numVariants = 0;
	private int numGenotyped = 0;
	private int [] scoresGroups;
	public Introgression(Sample sample, WindowScores w, int sampleIdx, String foreignGroup) {
		super(w.getSequenceName(), w.getFirst(), w.getLast());
		this.sample = sample;
		this.foreignGroup = foreignGroup;
		scoresGroups = new int [w.getNumGroups()];
		Arrays.fill(scoresGroups, 0);
		addWindow(w, sampleIdx);
		this.numVariants = w.getNumVariants();
		this.numGenotyped = w.getCountsGenotyped()[sampleIdx];
		
	}
	public void addWindow (WindowScores w, int sampleIdx) {
		if(getLast()<w.getLast()) setLast(w.getLast());
		numVariants+=w.getNumVariants();
		numGenotyped+=w.getCountsGenotyped()[sampleIdx];
		int [] scores = w.getScores()[sampleIdx];
		for(int j=0;j<scores.length;j++) scoresGroups[j]+=scores[j];
	}
	public Sample getSample() {
		return sample;
	}
	public String getForeignGroup() {
		return foreignGroup;
	}
	public void setForeignGroup(String foreignGroup) {
		this.foreignGroup = foreignGroup;
	}
	public int getNumVariants() {
		return numVariants;
	}
	public int getNumGenotyped() {
		return numGenotyped;
	}
	public int[] getScoresGroups() {
		return scoresGroups;
	}
}
