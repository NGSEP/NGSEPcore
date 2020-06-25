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

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.main.io.ParseUtils;
import ngsep.math.Distribution;
import ngsep.transcriptome.VariantFunctionalAnnotation;
import ngsep.transcriptome.VariantFunctionalAnnotationType;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.DiversityStatistics;
import ngsep.variants.GenomicVariant;
import ngsep.variants.SNV;

/**
 * This class calculates summary statistics from a VCF file
 * @author Jorge Duitama
 */
public class VCFSummaryStatisticsCalculator {
	
	// Constants for default values
	public static final int DEF_MIN_SAMPLES_GENOTYPED = 20;
	
	// Logging and progress
	private Logger log = Logger.getLogger(VCFSummaryStatisticsCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	//Parameters
	private String inputFile = null;
	private String outputFile = null;
	private int minSamplesGenotyped = DEF_MIN_SAMPLES_GENOTYPED;
	
	// Model attributes
	private static final String [] VARIANT_CATEGORIES= {"Biallelic SNVs","Biallelic Indels","Biallelic STRs","Other biallelic","Multiallelic SNVs","Multiallelic Indels","Multiallelic STRs","Other Multiallelic"};
	
	private List<String> sampleIds;
	//Counts for the summary section
	private VariantsBasicCounts [] summaryCounts = new VariantsBasicCounts[VARIANT_CATEGORIES.length];
	private int [] totalGenotypeCalls = new int [VARIANT_CATEGORIES.length];
	
	//MAF distribution per category 
	private Distribution [] mafDistribution = new Distribution[VARIANT_CATEGORIES.length];
	//MAF distribution per annotation. By now it is not distributed by categories
	private Map<String, Distribution> mafDistAnnBiallelicSNVs = new HashMap<>();
	private Map<String, Distribution> mafDistAnnBiallelicNonSNVs = new HashMap<>();
	
	//Distributions of genotyped accessions
	private Distribution [] genotypedAccessionsDistribution = new Distribution[VARIANT_CATEGORIES.length];
	//Counts for the data per sample
	private VariantsBasicCounts [][] countsPerSample;
	
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
	public int getMinSamplesGenotyped() {
		return minSamplesGenotyped;
	}
	public void setMinSamplesGenotyped(int minSamplesGenotyped) {
		this.minSamplesGenotyped = minSamplesGenotyped;
	}
	public void setMinSamplesGenotyped(String value) {
		this.setMinSamplesGenotyped((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		VCFSummaryStatisticsCalculator instance = new VCFSummaryStatisticsCalculator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	public void run() throws IOException {
		log.info("Minimum number of samples genotyped for calculation of population statistics: "+getMinSamplesGenotyped());
		if(inputFile==null) {
			log.info("Reading from standard input");
			if(outputFile == null) runStatistics(System.in, System.out);
			else {
				try (PrintStream out = new PrintStream(outputFile)) {
					runStatistics(System.in, out);
				}
			}
		} else {
			log.info("Reading from file: "+inputFile);
			if(outputFile == null) runStatistics(inputFile,System.out);
			else {
				try (PrintStream out = new PrintStream(outputFile)) {
					runStatistics(inputFile, out);
				}
			}
		}
		log.info("Process finished");
	}
	/**
	 * Calculates summary statistics
	 * @param filename Input VCF file
	 * @param out Stream where the output VCF will be written
	 * @throws IOException If the input file can not be read
	 */
	public void runStatistics(String filename, PrintStream out) throws IOException {
		try (VCFFileReader in = new VCFFileReader(filename)) { 
			runStatistics(in, out);
		}
	}
	public void runStatistics(InputStream fis, PrintStream out) throws IOException {
		try (VCFFileReader in = new VCFFileReader(fis)) {
			runStatistics(in, out);
		}	
	}
	public void runStatistics(VCFFileReader in, PrintStream out) throws IOException {
		if(log!=null)in.setLog(log);
		in.setLoadMode(VCFFileReader.LOAD_MODE_COPY_NUMBER);
		List<String> sampleIds = in.getHeader().getSampleIds();
		initStatistics(sampleIds);
		Iterator<VCFRecord> it = in.iterator();
		int n=0;
		while(it.hasNext()) {
			VCFRecord record = it.next();
			processRecord(record);
			n++;
			if (progressNotifier!=null && n%1000==0) {
				int progress = n/1000;
				if (!progressNotifier.keepRunning(progress)) {
					out.flush();
					return;
				}
			}
		}
		printStatistics(out);
	}

	private void initStatistics(List<String> sampleIds) {
		this.sampleIds = sampleIds;
		countsPerSample = new VariantsBasicCounts[VARIANT_CATEGORIES.length][sampleIds.size()];
		for(int i=0;i<VARIANT_CATEGORIES.length;i++) {
			summaryCounts[i] = new VariantsBasicCounts();
			totalGenotypeCalls[i] = 0;
			mafDistribution[i] = new Distribution(0, 0.5, 0.01);
			
			genotypedAccessionsDistribution[i] = new Distribution(0, sampleIds.size(), 1);
			for(int j=0;j<sampleIds.size();j++) {
				countsPerSample[i][j] = new VariantsBasicCounts();
			}
		}
		mafDistAnnBiallelicSNVs.put(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS, new Distribution(0, 0.5, 0.01));
		mafDistAnnBiallelicSNVs.put(VariantFunctionalAnnotationType.ANNOTATION_MISSENSE, new Distribution(0, 0.5, 0.01));
		mafDistAnnBiallelicSNVs.put(VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST, new Distribution(0, 0.5, 0.01));
		mafDistAnnBiallelicSNVs.put(VariantFunctionalAnnotationType.ANNOTATION_NONSENSE, new Distribution(0, 0.5, 0.01));
		mafDistAnnBiallelicSNVs.put(VariantFunctionalAnnotationType.ANNOTATION_START_LOST, new Distribution(0, 0.5, 0.01));
		mafDistAnnBiallelicNonSNVs.put(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS, new Distribution(0, 0.5, 0.01));
		mafDistAnnBiallelicNonSNVs.put(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS, new Distribution(0, 0.5, 0.01));
		mafDistAnnBiallelicNonSNVs.put(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL, new Distribution(0, 0.5, 0.01));
		mafDistAnnBiallelicNonSNVs.put(VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT, new Distribution(0, 0.5, 0.01));
	}
	public void processRecord(VCFRecord record) {
		GenomicVariant var = record.getVariant();
		//Gather variant characteristics to assign the category
		boolean isSNV = var.isSNV();
		boolean isBiallielic = var.isBiallelic();
		boolean isBiallelicSNV = isSNV && isBiallielic;
		boolean isIndel = var.getType()==GenomicVariant.TYPE_INDEL;
		boolean isSTR = var.getType()==GenomicVariant.TYPE_STR;
		if(isSNV && isIndel) {
			log.warning("Inconsistent Indel type in SNV at "+var.getSequenceName()+": "+var.getFirst()+". Ignoring record.");
			return;
		}
		if(isSNV && isSTR) {
			log.warning("Inconsistent STR type in SNV at "+var.getSequenceName()+": "+var.getFirst()+". Ignoring record.");
			return;
		}
		int idxVarType = calculateVariantCategory(isSNV, isIndel, isSTR, isBiallielic);
		List<CalledGenomicVariant> varCalls = record.getCalls();
		if(varCalls.size()!=countsPerSample[0].length) {
			throw new IllegalArgumentException("Inconsistent number of calls for variant "+record.getVariant().getSequenceName()+":"+record.getVariant().getFirst());
		}
		
		boolean isTransition = isBiallelicSNV && (var instanceof SNV) && ((SNV)var).isTransition();
		VariantFunctionalAnnotation annotation = record.getNGSEPFunctionalAnnotation();
		
		int populationStatus = 0;
		if(varCalls.size()==0) {
			//Variant without population information
			genotypedAccessionsDistribution[idxVarType].processDatapoint(0);
			summaryCounts[idxVarType].processGenotypeCall(VariantsBasicCounts.GENOTYPE_STATUS_HOMOALT, isTransition, annotation, populationStatus);
			return;
		}
		
		//Update diversity related statistics
		DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(varCalls, false);
		genotypedAccessionsDistribution[idxVarType].processDatapoint(divStats.getNumSamplesGenotyped());
		totalGenotypeCalls[idxVarType] += divStats.getNumSamplesGenotyped();
		double maf = divStats.getMaf();
		byte mafIdx = divStats.getMafIndex();
		byte wtIdx = divStats.getWtIndex();
		int [] alleleCounts = divStats.getAlleleCounts();
		if(divStats.getNumSamplesGenotyped()>=minSamplesGenotyped) {
			populationStatus = VariantsBasicCounts.POPULATION_STATUS_GENPOP;
			if(isBiallielic) {
				//if(call.getFirst()==181922)System.err.println("Call: "+i+" MAF: "+maf+". MAF idx: "+mafIdx+" wtIdx: "+wtIdx+" AC: "+alleleCounts[0]+" - "+ alleleCounts[1]+" status: "+popStatusSample);
				//Update MAF statistics
				mafDistribution[idxVarType].processDatapoint(maf);
				if(annotation!=null) {
					Distribution d;
					if(isSNV) {
						d = mafDistAnnBiallelicSNVs.get(annotation.getTypeName());
						if(d==null) {
							d = new Distribution(0, 0.5, 0.01);
							mafDistAnnBiallelicSNVs.put(annotation.getTypeName(), d);
						}
					} else {
						d = mafDistAnnBiallelicNonSNVs.get(annotation.getTypeName());
						if(d==null) {
							d = new Distribution(0, 0.5, 0.01);
							mafDistAnnBiallelicNonSNVs.put(annotation.getTypeName(), d);
						}
					}
					d.processDatapoint(maf);
				}
			}
		}
				
		summaryCounts[idxVarType].processGenotypeCall(VariantsBasicCounts.GENOTYPE_STATUS_HOMOALT, isTransition, annotation, populationStatus);
		//Update counts per sample
		for(int i=0;i<varCalls.size();i++) {
			CalledGenomicVariant call = varCalls.get(i);
			int genotypeStatus = calculateGenotypeStatus (call);
			int popStatusSample = populationStatus;
			//if(call.getFirst()==181922)System.err.println("Call: "+i+" wtIdx: "+wtIdx+" AC: "+alleleCounts[0]+" - "+ alleleCounts[1]+" status: "+popStatusSample);
			if(popStatusSample== VariantsBasicCounts.POPULATION_STATUS_GENPOP && wtIdx>=0) {
				byte [] calledAlleles = call.getIndexesCalledAlleles();
				short [] allelesCN = call.getAllelesCopyNumber();
				
				for(int j=0;j<calledAlleles.length;j++) {
					int callIdx = calledAlleles[j]; 
					if(callIdx != wtIdx) {
						popStatusSample = VariantsBasicCounts.POPULATION_STATUS_RARE;
						if(mafIdx == callIdx && allelesCN[callIdx] == alleleCounts[mafIdx]) {
							popStatusSample = VariantsBasicCounts.POPULATION_STATUS_UNIQUE;
						}
						break;
					}
				}
				//if(call.getFirst()==181922)System.err.println("Call: "+i+" MAF: "+maf+". MAF idx: "+mafIdx+" wtIdx: "+wtIdx+" AC: "+alleleCounts[0]+" - "+ alleleCounts[1]+" called: "+calledAlleles[0]+" hetero: "+ (calledAlleles.length>1)+ " ploidy a1: "+allelesCN[0]+" status: "+popStatusSample);
				//if(isBiallelicSNV && calledAlleles.length>0)System.out.println("MAF: "+maf+". MAF idx: "+mafIdx+" AC: "+alleleCounts[0]+" - "+ alleleCounts[1]+" called: "+calledAlleles[0]+" hetero: "+ (calledAlleles.length>1)+ " ploidy a1: "+calledAllelePloidies[0]+" status: "+populationStatus);
			}
			countsPerSample[idxVarType][i].processGenotypeCall(genotypeStatus, isTransition, annotation, popStatusSample);	
		}
		
	}
	private int calculateVariantCategory(boolean isSNV, boolean isIndel, boolean isSTR, boolean isBiallielic) {
		int idxVarType = 3;
		if(isSNV) {
			idxVarType = 0;
		} else if(isIndel) {
			idxVarType = 1;
		} else if (isSTR) {
			idxVarType = 2;
		}
		if(!isBiallielic) idxVarType+=4;
		return idxVarType;
	}

	private int calculateGenotypeStatus(CalledGenomicVariant call) {
		if(call.isUndecided()) return VariantsBasicCounts.GENOTYPE_STATUS_UNDECIDED;
		if(call.isHomozygousReference()) return VariantsBasicCounts.GENOTYPE_STATUS_HOMOREF;
		if(call.isHeterozygous()) return VariantsBasicCounts.GENOTYPE_STATUS_HETEROZYGOUS;
		return VariantsBasicCounts.GENOTYPE_STATUS_HOMOALT;
	}
	
	private void printStatistics(PrintStream out) {
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
		String [] annPrintOrder = new String [VariantFunctionalAnnotationType.getNumberSupportedTypes()];
		boolean [] printInSNVSummary = new boolean [VariantFunctionalAnnotationType.getNumberSupportedTypes()];
		boolean [] printInIndelSummary = new boolean [VariantFunctionalAnnotationType.getNumberSupportedTypes()];
		initPrintArraysdata(annPrintOrder,printInSNVSummary,printInIndelSummary);
		printGeneralSummary(out);
		printSummaryPerVariantType(out, fmt, annPrintOrder, printInSNVSummary, printInIndelSummary);
		printMAFDistributions(out, fmt, annPrintOrder);
		printSamplesGenotypedDistribution(out);
		printBiallelicSNVCountsPerSample(out, fmt);
		printBiallelicIndelCountsPerSample(out, fmt);
		printBiallelicSTRCountsPerSample(out, fmt);
		genericPrintCountsPerSample(out,"OTHER BIALLELIC VARIANTS COUNTS PER SAMPLE",3);
		genericPrintCountsPerSample(out,"MULTIALLELIC SNP COUNTS PER SAMPLE",4);
		genericPrintCountsPerSample(out,"MULTIALLELIC INDEL COUNTS PER SAMPLE",5);
		genericPrintCountsPerSample(out,"MULTIALLELIC STR COUNTS PER SAMPLE",6);
		genericPrintCountsPerSample(out,"OTHER MULTIALLELIC VARIANTS COUNTS PER SAMPLE",7);
		
	}
	public void printGeneralSummary(PrintStream out) {
		int l = VARIANT_CATEGORIES.length;
		out.println("GENERAL SUMMARY");
		out.print("Count");
		for(int i=0;i<l;i++) out.print("\t"+VARIANT_CATEGORIES[i]);
		out.println();
		out.print("Variants");
		for(int i=0;i<l;i++) out.print("\t"+summaryCounts[i].getGenotyped());
		out.println();
		out.print("Genotype calls");
		for(int i=0;i<l;i++) out.print("\t"+totalGenotypeCalls[i]);
		out.println();
		out.print("Coding variants");
		for(int i=0;i<l;i++) out.print("\t"+summaryCounts[i].getCodingTotalCount());
		out.println();
		out.print("Variants with at least "+minSamplesGenotyped+" samples genotyped");
		for(int i=0;i<l;i++) out.print("\t"+summaryCounts[i].getGenotypedPopCounts());
		out.println();
		out.println();
	}
	public void printSummaryPerVariantType(PrintStream out, DecimalFormat fmt, String[] annPrintOrder, boolean[] printInSNVSummary, boolean[] printInIndelSummary) {
		boolean [] printAll = new boolean[printInSNVSummary.length];
		Arrays.fill(printAll, true);
		printSummaryAnnotations("SUMMARY "+VARIANT_CATEGORIES[0], 0, out, annPrintOrder, printInSNVSummary);
		out.println("Transitions:\t"+summaryCounts[0].getTransitions()+"\tTr/Tv ratio:\t"+fmt.format(summaryCounts[0].getTrTvRatio()));
		out.println("Synonymous transitions:\t"+summaryCounts[0].getTransitionCount(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS)+"\tSynonymous Tr/Tv ratio:\t"+fmt.format(summaryCounts[0].getTrTvRatioSynonymous()));
		out.println("Non Synonymous transitions:\t"+summaryCounts[0].getNonSynonymousTransitionCount()+"\tNon synonymous Tr/Tv ratio:\t"+fmt.format(summaryCounts[0].getTrTvRatioNonSynonymous()));
		out.println();
		printSummaryAnnotations("SUMMARY "+VARIANT_CATEGORIES[1], 1, out, annPrintOrder, printInIndelSummary);
		out.println();
		printSummaryAnnotations("SUMMARY "+VARIANT_CATEGORIES[2], 2, out, annPrintOrder, printAll);
		out.println();
		printSummaryAnnotations("SUMMARY "+VARIANT_CATEGORIES[3], 3, out, annPrintOrder, printAll);
		out.println();
		printSummaryAnnotations("SUMMARY "+VARIANT_CATEGORIES[4], 4, out, annPrintOrder, printInSNVSummary);
		out.println();
		printSummaryAnnotations("SUMMARY "+VARIANT_CATEGORIES[5], 5, out, annPrintOrder, printInIndelSummary);
		out.println();
		printSummaryAnnotations("SUMMARY "+VARIANT_CATEGORIES[6], 6, out, annPrintOrder, printAll);
		out.println();
		printSummaryAnnotations("SUMMARY "+VARIANT_CATEGORIES[7], 7, out, annPrintOrder, printAll);
		out.println();
		
	}
	private void printSummaryAnnotations(String title, int idxSummaryCounts, PrintStream out, String[] annPrintOrder, boolean[] print) {
		out.println(title);
		out.println("Total:\t"+summaryCounts[idxSummaryCounts].getGenotyped());
		out.println("Coding:\t"+summaryCounts[idxSummaryCounts].getCodingTotalCount());
		for(int i=0;i<annPrintOrder.length;i++) {
			if(print[i]) out.println(getPrintName(annPrintOrder[i])+"\t"+summaryCounts[idxSummaryCounts].getTotalCount(annPrintOrder[i]));
		}
	}
	private String getPrintName(String annType) {
		if(VariantFunctionalAnnotationType.ANNOTATION_CODING.equals(annType)) return "Other coding";
		return annType;
	}
	public void printMAFDistributions(PrintStream out, DecimalFormat fmt, String[] annPrintOrder) {
		out.println("MAF DISTRIBUTIONS BIALLELIC SNVs WITH AT LEAST "+minSamplesGenotyped+" SAMPLES GENOTYPED");
		out.print("MAF\tTotal\t"+VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_MISSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_START_LOST+"/"+VariantFunctionalAnnotationType.ANNOTATION_NONSENSE);
		out.println("\tOther coding");
		int [] [] consolidatedDistSNVs = consolidateMAFDistributionsSNVs(annPrintOrder);
		printConsolidatedMAFDistributions(consolidatedDistSNVs, out, fmt);
		out.println();
		out.println("MAF DISTRIBUTIONS BIALLELIC INDELS AND STRs WITH AT LEAST "+minSamplesGenotyped+" SAMPLES GENOTYPED");
		out.print("MAF\tTotal indels\tTotal STRs\t"+VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL+"\t"+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT);
		out.println("\tOther coding");
		int [] [] consolidatedDistNonSNVs = consolidateMAFDistributionsNonSNVs(annPrintOrder);
		printConsolidatedMAFDistributions(consolidatedDistNonSNVs, out, fmt);
		out.println();
	}
	private int[][] consolidateMAFDistributionsSNVs(String [] printOrder) {
		double [] dist = mafDistribution[0].getDistribution();
		int [][] answer = new int [dist.length][5];
		for(int i=0;i<dist.length;i++) answer[i][0] = (int) dist[i];
		dist = mafDistAnnBiallelicSNVs.get(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][1] = (int) dist[i];
		dist = mafDistAnnBiallelicSNVs.get(VariantFunctionalAnnotationType.ANNOTATION_MISSENSE).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][2] = (int) dist[i];
		dist = mafDistAnnBiallelicSNVs.get(VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][2] += (int) dist[i];
		dist = mafDistAnnBiallelicSNVs.get(VariantFunctionalAnnotationType.ANNOTATION_NONSENSE).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][3] = (int) dist[i];
		dist = mafDistAnnBiallelicSNVs.get(VariantFunctionalAnnotationType.ANNOTATION_START_LOST).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][3] += (int) dist[i];
		//Fill others coding
		for(int i=0;i<dist.length;i++) answer[i][4] = 0;
		for(int j=8;j<=13;j++) {
			Distribution distO = mafDistAnnBiallelicSNVs.get(printOrder[j]);
			if(distO!=null) {
				dist = distO.getDistribution();
				for(int i=0;i<dist.length;i++) answer[i][4] += (int) dist[i];
			}
		}
		
		
		return answer;
	}
	private int[][] consolidateMAFDistributionsNonSNVs(String [] printOrder) {
		double [] dist = mafDistribution[1].getDistribution();
		int [][] answer = new int [dist.length][7];
		for(int i=0;i<dist.length;i++) answer[i][0] = (int) dist[i];
		dist = mafDistribution[2].getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][1] = (int) dist[i];
		dist = mafDistAnnBiallelicNonSNVs.get(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][2] = (int) dist[i];
		dist = mafDistAnnBiallelicNonSNVs.get(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][3] = (int) dist[i];
		dist = mafDistAnnBiallelicNonSNVs.get(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][4] = (int) dist[i];
		dist = mafDistAnnBiallelicNonSNVs.get(VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][5] = (int) dist[i];
		//Fill others coding
		for(int i=0;i<dist.length;i++) answer[i][6] = 0;
		for(int j=1;j<=4;j++) {
			Distribution distO = mafDistAnnBiallelicNonSNVs.get(printOrder[j]);
			if(distO!=null) {
				dist = distO.getDistribution();
				for(int i=0;i<dist.length;i++) answer[i][6] += (int) dist[i];
			}
		}
		for(int j=8;j<=13;j++) {
			Distribution distO = mafDistAnnBiallelicNonSNVs.get(printOrder[j]);
			if(distO!=null) {
				dist = distO.getDistribution();
				for(int i=0;i<dist.length;i++) answer[i][6] += (int) dist[i];
			}
		}
		
		return answer;
	}
	private void printConsolidatedMAFDistributions(int[][] consolidatedDist, PrintStream out, DecimalFormat fmt) {
		for(int i=0;i<consolidatedDist.length;i++) {
			double minMaf = 0.01*i;
			double maxMaf = 0.01*(i+1);
			out.print(fmt.format(minMaf));
			if(maxMaf<=0.5)out.print("-"+fmt.format(maxMaf));
			for(int j=0;j<consolidatedDist[i].length;j++) {
				out.print("\t"+consolidatedDist[i][j]);
			}
			out.println();
		}
	}
	public void printSamplesGenotypedDistribution(PrintStream out) {
		out.println("SAMPLES GENOTYPED DISTRIBUTIONS");
		out.print("Samples genotyped");
		for(int i=0;i<VARIANT_CATEGORIES.length;i++) out.print("\t"+VARIANT_CATEGORIES[i]);
		out.println();
		for(int i=0;i<=sampleIds.size();i++) {
			out.print(""+i);
			for(int j=0;j<VARIANT_CATEGORIES.length;j++) {
				int count = (int) Math.round(genotypedAccessionsDistribution[j].getDistribution()[i]);
				out.print("\t"+count);
			}
			out.println();
		}
		out.println();
	}
	public void printBiallelicSNVCountsPerSample(PrintStream out, DecimalFormat fmt) {
		out.println("SNP COUNTS PER SAMPLE");
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding\t"+VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_MISSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_NONSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_START_LOST);
		out.print("\tSplice donor/acceptor/region");
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_5P_UTR+"\t"+VariantFunctionalAnnotationType.ANNOTATION_3P_UTR);
		out.print("\tHeterozygous Coding\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS);
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_MISSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST);
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_NONSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_START_LOST);
		out.print("\tHeterozygous Splice donor/acceptor/region");
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_5P_UTR+"\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_3P_UTR);
		out.print("\tMinor allele known (At least "+minSamplesGenotyped+" samples genotyped)\tWith the minor allele\tWith a unique allele");
		out.print("\tTransitions\tTr/Tv\tHomozygous alternative transitions\tTr/Tv homozygous alternative\tHeterozygous transitions\tTr/Tv heterozygous");
		out.print("\tSynonymous transitions\tTr/Tv synonymous\tNon synonymous transitions\tTr/Tv non synonymous");
		out.println();
		for(int i=0;i<sampleIds.size();i++) {
			VariantsBasicCounts snpCountsSample = countsPerSample[0][i];
			out.print(sampleIds.get(i));
			out.print("\t"+snpCountsSample.getGenotyped());
			out.print("\t"+snpCountsSample.getNonReference());
			out.print("\t"+snpCountsSample.getHomozygousAlternative());
			out.print("\t"+snpCountsSample.getHeterozygous());
			out.print("\t"+snpCountsSample.getCodingTotalCount());
			out.print("\t"+snpCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS));
			out.print("\t"+(snpCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_MISSENSE)+snpCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST)));
			out.print("\t"+(snpCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_NONSENSE)+snpCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_START_LOST)));
			out.print("\t"+snpCountsSample.getTotalSpliceRegions());
			out.print("\t"+snpCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_5P_UTR));
			out.print("\t"+snpCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_3P_UTR));
			out.print("\t"+snpCountsSample.getCodingHeterozygousCount());
			out.print("\t"+snpCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS));
			out.print("\t"+(snpCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_MISSENSE)+snpCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST)));
			out.print("\t"+(snpCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_NONSENSE)+snpCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_START_LOST)));
			out.print("\t"+snpCountsSample.getHeterozygousSpliceRegions());
			out.print("\t"+snpCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_5P_UTR));
			out.print("\t"+snpCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_3P_UTR));
			out.print("\t"+snpCountsSample.getGenotypedPopCounts());
			out.print("\t"+snpCountsSample.getRareAllele());
			out.print("\t"+snpCountsSample.getUniqueAllele());
			out.print("\t"+snpCountsSample.getTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatio()));
			out.print("\t"+snpCountsSample.getHomozygousAlternativeTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioHomozygousAlternative()));
			out.print("\t"+snpCountsSample.getHeterozygousTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioHeterozygous()));
			out.print("\t"+snpCountsSample.getTransitionCount(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS));
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioSynonymous()));
			out.print("\t"+snpCountsSample.getNonSynonymousTransitionCount());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioNonSynonymous()));
			out.println();
		}
		out.println();
	}
	public void printBiallelicIndelCountsPerSample(PrintStream out, DecimalFormat fmt) {
		out.println("BIALLELIC INDEL COUNTS PER SAMPLE");
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding\t"+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL+"\t"+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT);
		out.print("\tSplice donor/acceptor/region");
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_5P_UTR+"\t"+VariantFunctionalAnnotationType.ANNOTATION_3P_UTR);
		out.print("\tHeterozygous Coding\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL+"\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS);
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT);
		out.print("\tHeterozygous Splice donor/acceptor/region");
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_5P_UTR+"\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_3P_UTR);
		out.print("\tMinor allele known (At least "+minSamplesGenotyped+" samples genotyped)\tWith the minor allele\tWith a unique allele");
		out.println();
		for(int i=0;i<sampleIds.size();i++) {
			VariantsBasicCounts indelCountsSample = countsPerSample[1][i];
			out.print(sampleIds.get(i));
			out.print("\t"+indelCountsSample.getGenotyped());
			out.print("\t"+indelCountsSample.getNonReference());
			out.print("\t"+indelCountsSample.getHomozygousAlternative());
			out.print("\t"+indelCountsSample.getHeterozygous());
			out.print("\t"+indelCountsSample.getCodingTotalCount());
			out.print("\t"+indelCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL));
			out.print("\t"+indelCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS));
			out.print("\t"+indelCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT));
			out.print("\t"+indelCountsSample.getTotalSpliceRegions());
			out.print("\t"+indelCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_5P_UTR));
			out.print("\t"+indelCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_3P_UTR));
			
			out.print("\t"+indelCountsSample.getCodingHeterozygousCount());
			out.print("\t"+indelCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL));
			out.print("\t"+indelCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS));
			out.print("\t"+indelCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT));
			out.print("\t"+indelCountsSample.getHeterozygousSpliceRegions());
			out.print("\t"+indelCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_5P_UTR));
			out.print("\t"+indelCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_3P_UTR));
			out.print("\t"+indelCountsSample.getGenotypedPopCounts());
			out.print("\t"+indelCountsSample.getRareAllele());
			out.print("\t"+indelCountsSample.getUniqueAllele());
			out.println();
			
		}
		out.println();
	}
	public void printBiallelicSTRCountsPerSample(PrintStream out, DecimalFormat fmt) {
		out.println("BIALLELIC STR COUNTS PER SAMPLE");
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding\t"+VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_MISSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_NONSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_START_LOST);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL+"\t"+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT);
		out.print("\tSplice donor/acceptor/region");
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_5P_UTR+"\t"+VariantFunctionalAnnotationType.ANNOTATION_3P_UTR);
		out.print("\tHeterozygous coding\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS);
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_MISSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST);
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_NONSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_START_LOST);
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL+"\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS);
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT);
		out.print("\tHeterozygous Splice donor/acceptor/region");
		out.print("\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_5P_UTR+"\tHeterozygous "+VariantFunctionalAnnotationType.ANNOTATION_3P_UTR);
		out.print("\tMinor allele known (At least "+minSamplesGenotyped+" samples genotyped)\tWith the minor allele\tWith a unique allele");
		out.println();
		for(int i=0;i<sampleIds.size();i++) {
			VariantsBasicCounts strCountsSample = countsPerSample[2][i];
			out.print(sampleIds.get(i));
			out.print("\t"+strCountsSample.getGenotyped());
			out.print("\t"+strCountsSample.getNonReference());
			out.print("\t"+strCountsSample.getHomozygousAlternative());
			out.print("\t"+strCountsSample.getHeterozygous());
			out.print("\t"+strCountsSample.getCodingTotalCount());
			out.print("\t"+strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS));
			out.print("\t"+(strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_MISSENSE)+strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST)));
			out.print("\t"+(strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_NONSENSE)+strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_START_LOST)));
			out.print("\t"+strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL));
			out.print("\t"+strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS));
			out.print("\t"+strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT));
			out.print("\t"+strCountsSample.getTotalSpliceRegions());
			out.print("\t"+strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_5P_UTR));
			out.print("\t"+strCountsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_3P_UTR));
			out.print("\t"+strCountsSample.getCodingHeterozygousCount());
			out.print("\t"+strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS));
			out.print("\t"+(strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_MISSENSE)+strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST)));
			out.print("\t"+(strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_NONSENSE)+strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_START_LOST)));
			out.print("\t"+strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL));
			out.print("\t"+strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS));
			out.print("\t"+strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT));
			out.print("\t"+strCountsSample.getHeterozygousSpliceRegions());
			out.print("\t"+strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_5P_UTR));
			out.print("\t"+strCountsSample.getHeterozygousCount(VariantFunctionalAnnotationType.ANNOTATION_3P_UTR));
			out.print("\t"+strCountsSample.getGenotypedPopCounts());
			out.print("\t"+strCountsSample.getRareAllele());
			out.print("\t"+strCountsSample.getUniqueAllele());
			out.println();
			
		}
		out.println();
	}
	public void genericPrintCountsPerSample(PrintStream out, String title, int index) {
		out.println(title);
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding\t"+VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_MISSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_NONSENSE+"/"+VariantFunctionalAnnotationType.ANNOTATION_START_LOST);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL+"\t"+VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS);
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT);
		out.print("\tSplice donor/acceptor/region");
		out.print("\t"+VariantFunctionalAnnotationType.ANNOTATION_5P_UTR+"\t"+VariantFunctionalAnnotationType.ANNOTATION_3P_UTR);
		out.print("\tMajor allele known (At least "+minSamplesGenotyped+" samples genotyped)\tWith a minor allele");
		out.println();
		for(int i=0;i<sampleIds.size();i++) {
			VariantsBasicCounts countsSample = countsPerSample[index][i];
			out.print(sampleIds.get(i));
			out.print("\t"+countsSample.getGenotyped());
			out.print("\t"+countsSample.getNonReference());
			out.print("\t"+countsSample.getHomozygousAlternative());
			out.print("\t"+countsSample.getHeterozygous());
			out.print("\t"+countsSample.getCodingTotalCount());
			out.print("\t"+countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS));
			out.print("\t"+(countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_MISSENSE)+countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST)));
			out.print("\t"+(countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_NONSENSE)+countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_START_LOST)));
			out.print("\t"+countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL));
			out.print("\t"+countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS));
			out.print("\t"+countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT));
			out.print("\t"+countsSample.getTotalSpliceRegions());
			out.print("\t"+countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_5P_UTR));
			out.print("\t"+countsSample.getTotalCount(VariantFunctionalAnnotationType.ANNOTATION_3P_UTR));
			out.print("\t"+countsSample.getGenotypedPopCounts());
			out.print("\t"+countsSample.getRareAllele());
			out.println();
		}
		out.println();
	}
	
	private void initPrintArraysdata(String[] annPrintOrder, boolean[] printInSNVSummary, boolean[] printInIndelSummary) {
		Arrays.fill(printInSNVSummary, true);
		Arrays.fill(printInIndelSummary, true);
		annPrintOrder[0] = VariantFunctionalAnnotationType.ANNOTATION_SYNONYMOUS;
		printInIndelSummary[0] = false;
		annPrintOrder[1] = VariantFunctionalAnnotationType.ANNOTATION_MISSENSE;
		printInIndelSummary[1] = false;
		annPrintOrder[2] = VariantFunctionalAnnotationType.ANNOTATION_STOP_LOST;
		printInIndelSummary[2] = false;
		annPrintOrder[3] = VariantFunctionalAnnotationType.ANNOTATION_NONSENSE;
		printInIndelSummary[3] = false;
		annPrintOrder[4] = VariantFunctionalAnnotationType.ANNOTATION_START_LOST;
		printInIndelSummary[4] = false;
		
		annPrintOrder[5] = VariantFunctionalAnnotationType.ANNOTATION_INFRAME_DEL;
		printInSNVSummary[5] = false;
		annPrintOrder[6] = VariantFunctionalAnnotationType.ANNOTATION_INFRAME_INS;
		printInSNVSummary[6] = false;
		annPrintOrder[7] = VariantFunctionalAnnotationType.ANNOTATION_FRAMESHIFT;
		printInSNVSummary[7] = false;
		
		annPrintOrder[8] = VariantFunctionalAnnotationType.ANNOTATION_SPLICE_DONOR;
		annPrintOrder[9] = VariantFunctionalAnnotationType.ANNOTATION_SPLICE_ACCEPTOR;
		annPrintOrder[10] = VariantFunctionalAnnotationType.ANNOTATION_EXONIC_SPLICE_REGION;
		annPrintOrder[11] = VariantFunctionalAnnotationType.ANNOTATION_SPLICE_REGION;
		
		annPrintOrder[12] = VariantFunctionalAnnotationType.ANNOTATION_5P_UTR;
		annPrintOrder[13] = VariantFunctionalAnnotationType.ANNOTATION_3P_UTR;
		annPrintOrder[14] = VariantFunctionalAnnotationType.ANNOTATION_CODING;

		annPrintOrder[15] = VariantFunctionalAnnotationType.ANNOTATION_NONCODINGRNA;
		
		annPrintOrder[16] = VariantFunctionalAnnotationType.ANNOTATION_UPSTREAM;
		annPrintOrder[17] = VariantFunctionalAnnotationType.ANNOTATION_DOWNSTREAM;
		annPrintOrder[18] = VariantFunctionalAnnotationType.ANNOTATION_INTRON;
		annPrintOrder[19] = VariantFunctionalAnnotationType.ANNOTATION_INTERGENIC;
		
	}
}
