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
package ngsep.variants.imputation;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.hmm.RecombinationHMM;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;


public class GenotypeImputer {
	
	// Constants for default values
	public static final int DEF_NUM_HAPLOTYPE_CLUSTERS = 8;
	public static final int DEF_WINDOW_SIZE = 5000;
	public static final int DEF_OVERLAP = 50;
	public static final double DEF_AVG_CM_PER_KBP = 0.001;
	
	// Logging and progress
	private Logger log = Logger.getLogger(GenotypeImputer.class.getName());
	private ProgressNotifier progressNotifier=null;
	private int progress = 0;
	
	// Parameters
	private String inputFile = null;
	private String outputPrefix = null;
	private List<String> parentIds = new ArrayList<String>() ;
	private int numHaplotypeClusters = DEF_NUM_HAPLOTYPE_CLUSTERS;
	private int windowSize = DEF_WINDOW_SIZE;
	private int overlap = DEF_OVERLAP;
	private double avgCMPerKbp = DEF_AVG_CM_PER_KBP;
	private boolean skipTransitionsTraining = false;
	private boolean inbredParents = false;
	private boolean inbredSamples = false;
	
	// Model attributes
	private PrintStream outAssignments;
	 
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
	
	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	
	public List<String> getParentIds() {
		return parentIds;
	}
	public void setParentIds(List<String> parentIds) {
		this.parentIds = parentIds;
	}
	public void setParentIds(String parentIdsStr) {
		this.setParentIds(Arrays.asList(parentIdsStr.split(",")));
	}
	
	public int getNumHaplotypeClusters() {
		return numHaplotypeClusters;
	}
	public void setNumHaplotypeClusters(int numHaplotypeClusters) {
		this.numHaplotypeClusters = numHaplotypeClusters;
	}
	public void setNumHaplotypeClusters(String value) {
		setNumHaplotypeClusters((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getWindowSize() {
		return windowSize;
	}
	public void setWindowSize(int windowSize) {
		this.windowSize = windowSize;
	}
	public void setWindowSize(String value) {
		setWindowSize((int)OptionValuesDecoder.decode(value, Integer.class));
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
	
	public double getAvgCMPerKbp() {
		return avgCMPerKbp;
	}
	public void setAvgCMPerKbp(double avgCMPerKbp) {
		this.avgCMPerKbp = avgCMPerKbp;
	}
	public void setAvgCMPerKbp(String value) {
		setAvgCMPerKbp((double)OptionValuesDecoder.decode(value, Double.class));
	}

	public boolean isSkipTransitionsTraining() {
		return skipTransitionsTraining;
	}
	public void setSkipTransitionsTraining(boolean skipTransitionsTraining) {
		this.skipTransitionsTraining = skipTransitionsTraining;
	}
	public void setSkipTransitionsTraining(Boolean skipTransitionsTraining) {
		this.setSkipTransitionsTraining(skipTransitionsTraining.booleanValue());
	}
	
	public boolean isInbredParents() {
		return inbredParents;
	}
	public void setInbredParents(boolean inbredParents) {
		this.inbredParents = inbredParents;
	}
	public void setInbredParents(Boolean inbredParents) {
		this.setInbredParents(inbredParents.booleanValue());
	}

	public boolean isInbredSamples() {
		return inbredSamples;
	}
	public void setInbredSamples(boolean inbredSamples) {
		this.inbredSamples = inbredSamples;
	}
	public void setInbredSamples(Boolean inbredSamples) {
		this.setInbredSamples(inbredSamples.booleanValue());
	}

	public PrintStream getOutAssignments() {
		return outAssignments;
	}
	public void setOutAssignments(PrintStream outAssignments) {
		this.outAssignments = outAssignments;
	}

	public static void main(String[] args) throws Exception {
		GenotypeImputer instance = new GenotypeImputer();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
		
	}

	public void run () throws IOException {
		logParameters();
		if(outputPrefix==null) throw new IOException("The prefix for the output files is a required parameter");
		outAssignments = null;
		try (PrintStream outGenotypes = new PrintStream(outputPrefix+"_imputed.vcf")){
			if(parentIds!=null && parentIds.size()>1) outAssignments = new PrintStream(outputPrefix+"_assignments.txt");
			if(inputFile!=null) impute(inputFile, outGenotypes);
			else {
				try (VCFFileReader reader = new VCFFileReader(System.in)) {
					impute(reader, outGenotypes);
				}
			}
		} finally {
			if(outAssignments!=null) outAssignments.close();
		}
		log.info("Process finished");
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		if(inputFile != null) out.println("Input file: "+inputFile);
		else out.println("Read variants from standard input");
		out.println("Prefix for output files: "+outputPrefix);
		if(parentIds!=null) out.println("Parents: "+parentIds);
		out.println("Number of haplotype clusters: "+getNumHaplotypeClusters());
		out.println("Window size: "+getWindowSize());
		out.println("Overlap: "+getOverlap());
		out.println("Average centimorgans per Kbp: "+avgCMPerKbp);
		if(skipTransitionsTraining) out.println("Transitions will not be modified during the HMM training");
		if(inbredParents) out.println("Parents of the population are assumed to be inbred");
		if(inbredSamples) out.println("Samples of the population are assumed to be inbred. All imputed genotype calls will be homozygous");
		log.info(""+os.toString());
	}
	public void impute(String filename, PrintStream outGenotypes) throws IOException {
		try (VCFFileReader reader = new VCFFileReader(filename)) {
			impute(reader, outGenotypes);
		}
	}
	public void impute(VCFFileReader reader, PrintStream outGenotypes) throws IOException {
		List<VCFRecord> records = new ArrayList<VCFRecord>();
		List<VCFRecord> lastRecords = new ArrayList<VCFRecord>();
		VCFFileWriter writer = new VCFFileWriter();
		
		if(log!=null) reader.setLog(log);
		//Ids of the samples
		VCFFileHeader header = reader.getHeader();
		writer.printHeader(header,outGenotypes);
		List<Sample> samples = header.getSamples();
		String lastSeqName = null;
		Iterator<VCFRecord> it = reader.iterator();
		while(it.hasNext()) {
			VCFRecord record = it.next();
			GenomicVariant var = record.getVariant();
			if(!(var instanceof SNV)) continue;
			boolean sequenceChange = !var.getSequenceName().equals(lastSeqName); 
			if(sequenceChange || records.size() == windowSize) {
				if(lastSeqName!=null) {
					processRecords(records,samples,lastRecords, writer, outGenotypes,sequenceChange);
					for(VCFRecord r:records) writer.printVCFRecord(r, outGenotypes);
					progress++;
					if(progressNotifier!=null && !progressNotifier.keepRunning(progress)) return;
				}
				lastSeqName = var.getSequenceName();
			}
			records.add(record);
		}
		if(lastSeqName!=null) {
			processRecords(records, samples, lastRecords, writer, outGenotypes, true);
		}
	}
	
	private void processRecords(List<VCFRecord> currentRecords, List<Sample> samples, List<VCFRecord> lastRecords, VCFFileWriter writer, PrintStream outGenotypes, boolean sequenceChange) {
		List<VCFRecord> recordsImpute = calculateRecordsImpute (currentRecords,lastRecords);
		Map<String, List<CalledSNV>> genotypes = convertToCalledGenotypes(samples, recordsImpute);
		imputeGenotypes(genotypes);
		int printStart = 0;
		if(lastRecords.size()>0) {
			printStart+=overlap;
		}
		int printEnd = recordsImpute.size();
		if(!sequenceChange) {
			printEnd-=overlap;
		}
		
		//Print records
		for(int i = printStart;i<printEnd;i++) {
			VCFRecord r = recordsImpute.get(i);
			writer.printVCFRecord(r, outGenotypes);
		}
		//Update last records
		lastRecords.clear();
		if(!sequenceChange) {
			lastRecords.addAll(currentRecords);
		}
		//Clean records for next window
		currentRecords.clear();
	}

	private List<VCFRecord> calculateRecordsImpute(List<VCFRecord> currentRecords, List<VCFRecord> lastRecords) {
		List<VCFRecord> answer = new ArrayList<>(windowSize+2*overlap);
		if(lastRecords.size()>0) {
			for(int i=lastRecords.size()-2*overlap;i<lastRecords.size();i++) {
				answer.add(lastRecords.get(i));
			}
		}
		answer.addAll(currentRecords);
		return answer;
	}

	private Map<String, List<CalledSNV>> convertToCalledGenotypes(List<Sample> samples, List<VCFRecord> recordsImpute) {
		Map<String, List<CalledSNV>> genotypes = new TreeMap<String, List<CalledSNV>>();
		for(Sample sample:samples) genotypes.put(sample.getId(), new ArrayList<CalledSNV>());
		for(VCFRecord record:recordsImpute) {
			GenomicVariant var = record.getVariant();
			SNV snv = (SNV) var;
			List<CalledGenomicVariant> genotypeCalls = record.getCalls();
			for(int i=0;i<genotypeCalls.size();i++) {
				String sampleId = samples.get(i).getId();
				CalledGenomicVariant genotypeCall = genotypeCalls.get(i);
				CalledSNV csnv;
				if(genotypeCall instanceof CalledSNV) csnv = (CalledSNV)genotypeCall;
				else csnv = new CalledSNV(snv, CalledSNV.GENOTYPE_UNDECIDED);
				genotypes.get(sampleId).add(csnv);
			}
		}
		return genotypes;
	}

	public void imputeGenotypes(Map<String,List<CalledSNV>> genotypes) {
		progress = 0;
		if(inbredSamples) imputeGenotypesHMMInbreds(genotypes);
		else  imputeGenotypesHMMDiploid(genotypes);
	}

	public void imputeGenotypesHMMInbreds(Map<String, List<CalledSNV>> genotypes) {
		GenotypeImputationHMM  hmm = GenotypeImputationHMM.createHMM(genotypes, parentIds, numHaplotypeClusters, inbredParents);
		hmm.setLog(log);
		hmm.setAvgCMPerKbp(avgCMPerKbp);
		hmm.setSkipTransitionsTraining(skipTransitionsTraining);
		hmm.setTrainingData(makeTrainingDataWithHomozygous(genotypes));
		if(progressNotifier!=null) {
			progress++;
			if(!progressNotifier.keepRunning(progress)) return;
		}
		
		int [][][] outClusters = new int [hmm.getStartsBaumWelch()][genotypes.size()][hmm.getSteps()];
		hmm.imputeGenotypes(genotypes,outClusters);
		List<CalledSNV> snvs = genotypes.values().iterator().next();
		//TODO: Conciliate more than one run of the Baum-Welch
		if(outAssignments!=null) printClusters(genotypes.keySet(),snvs,outClusters[0],hmm);
	}
	
	public void imputeGenotypesHMMDiploid(Map<String, List<CalledSNV>> genotypes) {
		DiploidGenotypeImputationHMM  hmm = DiploidGenotypeImputationHMM.createHMM(genotypes, parentIds, numHaplotypeClusters, inbredParents);
		hmm.setAvgCMPerKbp(avgCMPerKbp);
		hmm.setSkipTransitionsTraining(skipTransitionsTraining);
		hmm.setLog(log);
		hmm.setTrainingData(makeTrainingDataWithHomozygous(genotypes));
		
		if(progressNotifier!=null) {
			progress++;
			if(!progressNotifier.keepRunning(progress)) return;
		}
		int [][][] outClusters = new int [hmm.getStartsBaumWelch()][genotypes.size()][hmm.getSteps()];
		hmm.imputeGenotypes(genotypes,outClusters);
		List<CalledSNV> snvs = genotypes.values().iterator().next();
		//TODO: Conciliate more than one run of the Baum-Welch
		if(outAssignments!=null) printClusters(genotypes.keySet(),snvs,outClusters[0],hmm);
	}
	
	private void printClusters(Set<String>sampleIds, List<CalledSNV> snvs, int [][] outClusters,RecombinationHMM hmm) {
		int m = hmm.getSteps();
		outAssignments.print("Chr\tPos");
		for(String sampleId:sampleIds) outAssignments.print("\t"+sampleId);
		outAssignments.println();
		for (int j=0;j<m;j++) {
			CalledSNV snv = snvs.get(j);
			outAssignments.print(snv.getSequenceName()+"\t"+snv.getPosition());
			
			for(int i=0;i<sampleIds.size();i++) {
				int assignment = outClusters[i][j];
				String stateId = hmm.getState(assignment).getId();
				if(stateId==null) outAssignments.print("\t"+assignment);
				else outAssignments.print("\t"+stateId);
			}
			outAssignments.println();
		}
	}
	
	private List<List<? extends Object>> makeTrainingDataWithHomozygous(Map<String, List<CalledSNV>> genotypes) {
		List<List<? extends Object>> haplotypes = new ArrayList<>();
		for (String sampleId :genotypes.keySet()) {
			//System.out.println("Next id training data: "+sampleId);
			List<? extends Object> hapSample = HaplotypeClustersHMM.makeHaplotypeWithHomozygous(genotypes.get(sampleId));
			haplotypes.add(hapSample);
		}
		return haplotypes;
	}
	
}
