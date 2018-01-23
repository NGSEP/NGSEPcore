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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.main.io.ParseUtils;
import ngsep.math.Distribution;
import ngsep.transcriptome.Transcriptome;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.DiversityStatistics;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.SNV;

/**
 * This class calculates summary statistics from a VCF file
 * @author Jorge Duitama
 */
public class VCFSummaryStatisticsCalculator {
	
	private Logger log = Logger.getLogger(VCFSummaryStatisticsCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;

	private static final String [] VARIANT_CATEGORIES= {"Biallelic SNVs","Biallelic Indels","Biallelic STRs","Other biallelic","Multiallelic SNVs","Multiallelic Indels","Multiallelic STRs","Other Multiallelic"};
	private int minSamplesGenotyped = 20;
	private VariantsBasicCounts [] summaryCounts = new VariantsBasicCounts[VARIANT_CATEGORIES.length];
	private List<String> sampleIds;
	private VariantsBasicCounts [][] countsPerSample;
	private int [] totalGenotypeCalls = new int [VARIANT_CATEGORIES.length];
	private Distribution [] genotypedAccessionsDistribution = new Distribution[VARIANT_CATEGORIES.length];
	private Distribution [] mafDistribution = new Distribution[VARIANT_CATEGORIES.length];
	private Map<String, Distribution> mafDistributionAnns = new HashMap<>();
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		if (args.length == 0 || args[0].equals("-h") ||args[0].equals("--help")){
			CommandsDescriptor.getInstance().printHelp(VCFSummaryStatisticsCalculator.class);
			return;
		}
		VCFSummaryStatisticsCalculator instance = new VCFSummaryStatisticsCalculator();
		boolean systemInput = false;
		int i=0;
		while(i<args.length && args[i].charAt(0)=='-') {
			if ("-m".equals(args[i])) {
				i++;
				instance.minSamplesGenotyped = Integer.parseInt(args[i]);
			} else if ("-".equals(args[i])) {
				systemInput = true;
				break;
			} else {
				System.err.println("Unrecognized option "+args[i]);
				CommandsDescriptor.getInstance().printHelp(VCFSummaryStatisticsCalculator.class);
				return;
			}
			i++;
		}
		if(systemInput) {
			instance.runStatistics(System.in, System.out);
		} else {
			String filename = args[i];
			instance.runStatistics(filename, System.out);
		}
	}
	/**
	 * Calculates summary statistics
	 * @param filename Input VCF file
	 * @param out Stream where the output VCF will be written
	 * @throws IOException If the input file can not be read
	 */
	public void runStatistics(String filename, PrintStream out) throws IOException {
		VCFFileReader in = null;
		try {
			in = new VCFFileReader(filename);
			runStatistics(in, out);
		} finally {
			if(in!=null) in.close();
		}
	}
	public void runStatistics(InputStream fis, PrintStream out) throws IOException {
		VCFFileReader in = null;
		try {
			in = new VCFFileReader(fis);
			runStatistics(in, out);
		} finally {
			if(in!=null) in.close();
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
			genotypedAccessionsDistribution[i] = new Distribution(0, sampleIds.size(), 1);
			mafDistribution[i] = new Distribution(0, 0.5, 0.01);
			
			for(int j=0;j<sampleIds.size();j++) {
				countsPerSample[i][j] = new VariantsBasicCounts();
			}
		}
		mafDistributionAnns.put(Transcriptome.ANNOTATION_SYNONYMOUS, new Distribution(0, 0.5, 0.01));
		mafDistributionAnns.put(Transcriptome.ANNOTATION_MISSENSE, new Distribution(0, 0.5, 0.01));
		mafDistributionAnns.put(Transcriptome.ANNOTATION_NONSENSE, new Distribution(0, 0.5, 0.01));
		mafDistributionAnns.put(Transcriptome.ANNOTATION_FRAMESHIFT, new Distribution(0, 0.5, 0.01));
	}
	public void processRecord(VCFRecord record) {
		GenomicVariant var = record.getVariant();
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
		int idxVarType = 3;
		if(isSNV) {
			idxVarType = 0;
		} else if(isIndel) {
			idxVarType = 1;
		} else if (isSTR) {
			idxVarType = 2;
		}
		if(!isBiallielic) idxVarType+=4;
		List<CalledGenomicVariant> varCalls = record.getCalls();
		if(varCalls.size()!=countsPerSample[0].length) {
			throw new IllegalArgumentException("Inconsistent number of calls for variant "+record.getVariant().getSequenceName()+":"+record.getVariant().getFirst());
		}
		
		boolean isTransition = isBiallelicSNV && (var instanceof SNV) && ((SNV)var).isTransition();
		String annotation = null;
		GenomicVariantAnnotation ann = record.getInfoField(GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ANNOTATION);
		if(ann!=null) annotation = ann.getValue().toString();
		int populationStatus = 0;
		if(varCalls.size()==0) {
			genotypedAccessionsDistribution[idxVarType].processDatapoint(0);
			summaryCounts[idxVarType].processGenotypeCall(VariantsBasicCounts.GENOTYPE_STATUS_HOMOALT, isTransition, annotation, populationStatus);
			return;
		}
		DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(varCalls, false);
		genotypedAccessionsDistribution[idxVarType].processDatapoint(divStats.getNumSamplesGenotyped());
		totalGenotypeCalls[idxVarType] += divStats.getNumSamplesGenotyped();
		double maf = divStats.getMaf();
		byte mafIdx = divStats.getMafIndex();
		byte wtIdx = divStats.getWtIndex();
		int [] alleleCounts = divStats.getAlleleCounts();
		if(isBiallielic && divStats.getNumSamplesGenotyped()>=minSamplesGenotyped) {
			populationStatus = VariantsBasicCounts.POPULATION_STATUS_GENPOP;
			mafDistribution[idxVarType].processDatapoint(maf);
			if(annotation!=null) {
				Distribution d = mafDistributionAnns.get(annotation);
				if(d==null) {
					d = new Distribution(0, 0.5, 0.01);
					mafDistributionAnns.put(annotation, d);
				}
				d.processDatapoint(maf);
			}
		}
				
		summaryCounts[idxVarType].processGenotypeCall(VariantsBasicCounts.GENOTYPE_STATUS_HOMOALT, isTransition, annotation, populationStatus);
		for(int i=0;i<varCalls.size();i++) {
			CalledGenomicVariant call = varCalls.get(i);
			int genotypeStatus = calculateGenotypeStatus (call);
			int popStatusSample = populationStatus;
			//if(call.getFirst()==181922)System.err.println("Call: "+i+" MAF: "+maf+". MAF idx: "+mafIdx+" wtIdx: "+wtIdx+" AC: "+alleleCounts[0]+" - "+ alleleCounts[1]+" status: "+popStatusSample);
			if(popStatusSample== VariantsBasicCounts.POPULATION_STATUS_GENPOP && wtIdx>=0) {
				byte [] calledAlleles = call.getIndexesCalledAlleles();
				byte [] allelesCN = call.getAllelesCopyNumber();
				
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

	private int calculateGenotypeStatus(CalledGenomicVariant call) {
		if(call.isUndecided()) return VariantsBasicCounts.GENOTYPE_STATUS_UNDECIDED;
		if(call.isHomozygousReference()) return VariantsBasicCounts.GENOTYPE_STATUS_HOMOREF;
		if(call.isHeterozygous()) return VariantsBasicCounts.GENOTYPE_STATUS_HETEROZYGOUS;
		return VariantsBasicCounts.GENOTYPE_STATUS_HOMOALT;
	}
	
	private void printStatistics(PrintStream out) {
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
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
		out.println("SUMMARY FUNCTIONAL BIALLELIC SNVs");
		out.println("Synonymous:\t"+summaryCounts[0].getTotalCount(Transcriptome.ANNOTATION_SYNONYMOUS));
		out.println("Start loss:\t"+summaryCounts[0].getTotalCount(Transcriptome.ANNOTATION_START_LOSS));
		out.println("Missense:\t"+summaryCounts[0].getTotalCount(Transcriptome.ANNOTATION_MISSENSE));
		out.println("Stop loss:\t"+summaryCounts[0].getTotalCount(Transcriptome.ANNOTATION_STOP_LOSS));
		out.println("Non sense:\t"+summaryCounts[0].getTotalCount(Transcriptome.ANNOTATION_NONSENSE));
		out.println("Transitions:\t"+summaryCounts[0].getTransitions()+"\tTr/Tv ratio:\t"+fmt.format(summaryCounts[0].getTrTvRatio()));
		out.println("Synonymous transitions:\t"+summaryCounts[0].getTransitionCount(Transcriptome.ANNOTATION_SYNONYMOUS)+"\tSynonymous Tr/Tv ratio:\t"+fmt.format(summaryCounts[0].getTrTvRatioSynonymous()));
		out.println("Non Synonymous transitions:\t"+summaryCounts[0].getNonSynonymousTransitionCount()+"\tNon synonymous Tr/Tv ratio:\t"+fmt.format(summaryCounts[0].getTrTvRatioNonSynonymous()));
		out.println();
		out.println("SUMMARY OTHER FUNCTIONAL BIALLELIC VARIANTS");
		out.println("Frameshift indels:\t"+summaryCounts[1].getTotalCount(Transcriptome.ANNOTATION_FRAMESHIFT)+"\tPCT:\t"+fmt.format(summaryCounts[1].getPCTFrameshift()));
		out.println("Frameshift STRs:\t"+summaryCounts[2].getTotalCount(Transcriptome.ANNOTATION_FRAMESHIFT)+"\tPCT:\t"+fmt.format(summaryCounts[2].getPCTFrameshift()));
		out.println("Frameshift other:\t"+summaryCounts[3].getTotalCount(Transcriptome.ANNOTATION_FRAMESHIFT)+"\tPCT:\t"+fmt.format(summaryCounts[3].getPCTFrameshift()));
		out.println();
		out.println("MAF DISTRIBUTIONS BIALLELIC VARIANTS WITH AT LEAST "+minSamplesGenotyped+" SAMPLES GENOTYPED");
		out.println("MAF\tBiallelic SNPs\tSynonymous\tMissense\tNonsense\tOther biallelic variants\tFrameshift indels");
		double [] [] consolidatedDist = consolidateMAFDistributions();
		for(int i=0;i<consolidatedDist.length;i++) {
			double minMaf = 0.01*i;
			double maxMaf = 0.01*(i+1);
			out.print(fmt.format(minMaf));
			if(maxMaf<=0.5)out.print("-"+fmt.format(maxMaf));
			for(int j=0;j<consolidatedDist[i].length;j++) {
				out.print("\t"+fmt.format(consolidatedDist[i][j]));
			}
			out.println();
		}
		out.println();
		out.println("SAMPLES GENOTYPED DISTRIBUTIONS");
		out.print("Samples genotyped");
		for(int i=0;i<l;i++) out.print("\t"+VARIANT_CATEGORIES[i]);
		out.println();
		for(int i=0;i<=sampleIds.size();i++) {
			out.print(""+i);
			for(int j=0;j<l;j++) {
				int count = (int) Math.round(genotypedAccessionsDistribution[j].getDistribution()[i]);
				out.print("\t"+count);
			}
			out.println();
		}
		out.println();
		out.println("SNP COUNTS PER SAMPLE");
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding\tSynonymous\tStart loss\tMissense\tStop loss\tNonsense");
		out.print("\tHeterozygous coding\tHeterozygous synonymous\tHeterozygous start loss\tHeterozygous missense\tHeterozygous stop loss\tHeterozygous nonsense");
		out.print("\tMinor allele known\tWith the minor allele\tWith a unique allele");
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
			out.print("\t"+snpCountsSample.getTotalCount(Transcriptome.ANNOTATION_SYNONYMOUS));
			out.print("\t"+snpCountsSample.getTotalCount(Transcriptome.ANNOTATION_START_LOSS));
			out.print("\t"+snpCountsSample.getTotalCount(Transcriptome.ANNOTATION_MISSENSE));
			out.print("\t"+snpCountsSample.getTotalCount(Transcriptome.ANNOTATION_STOP_LOSS));
			out.print("\t"+snpCountsSample.getTotalCount(Transcriptome.ANNOTATION_NONSENSE));
			out.print("\t"+snpCountsSample.getCodingHeterozygousCount());
			out.print("\t"+snpCountsSample.getHeterozygousCount(Transcriptome.ANNOTATION_SYNONYMOUS));
			out.print("\t"+snpCountsSample.getHeterozygousCount(Transcriptome.ANNOTATION_START_LOSS));
			out.print("\t"+snpCountsSample.getHeterozygousCount(Transcriptome.ANNOTATION_MISSENSE));
			out.print("\t"+snpCountsSample.getHeterozygousCount(Transcriptome.ANNOTATION_STOP_LOSS));
			out.print("\t"+snpCountsSample.getHeterozygousCount(Transcriptome.ANNOTATION_NONSENSE));
			out.print("\t"+snpCountsSample.getGenotypedPopCounts());
			out.print("\t"+snpCountsSample.getRareAllele());
			out.print("\t"+snpCountsSample.getUniqueAllele());
			out.print("\t"+snpCountsSample.getTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatio()));
			out.print("\t"+snpCountsSample.getHomozygousAlternativeTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioHomozygousAlternative()));
			out.print("\t"+snpCountsSample.getHeterozygousTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioHeterozygous()));
			out.print("\t"+snpCountsSample.getTransitionCount(Transcriptome.ANNOTATION_SYNONYMOUS));
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioSynonymous()));
			out.print("\t"+snpCountsSample.getNonSynonymousTransitionCount());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioNonSynonymous()));
			out.println();
		}
		out.println();
		out.println("BIALLELIC INDEL COUNTS PER SAMPLE");
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding\tFrameshift\tPercentage frameshift\tCoding heterozygous\tPercentage frameshift heterozygous");
		out.print("\tMinor allele known\tWith the minor allele\tWith a unique allele");
		out.println();
		for(int i=0;i<sampleIds.size();i++) {
			VariantsBasicCounts indelCountsSample = countsPerSample[1][i];
			out.print(sampleIds.get(i));
			out.print("\t"+indelCountsSample.getGenotyped());
			out.print("\t"+indelCountsSample.getNonReference());
			out.print("\t"+indelCountsSample.getHomozygousAlternative());
			out.print("\t"+indelCountsSample.getHeterozygous());
			out.print("\t"+indelCountsSample.getCodingTotalCount());
			out.print("\t"+indelCountsSample.getTotalCount(Transcriptome.ANNOTATION_FRAMESHIFT));
			out.print("\t"+fmt.format(indelCountsSample.getPCTFrameshift()));
			out.print("\t"+indelCountsSample.getCodingHeterozygousCount());
			out.print("\t"+indelCountsSample.getHeterozygousCount(Transcriptome.ANNOTATION_FRAMESHIFT));
			out.print("\t"+fmt.format(indelCountsSample.getPCTFrameshiftHeterozygous()));
			out.print("\t"+indelCountsSample.getGenotypedPopCounts());
			out.print("\t"+indelCountsSample.getRareAllele());
			out.print("\t"+indelCountsSample.getUniqueAllele());
			out.println();
			
		}
		out.println();
		out.println("BIALLELIC STR COUNTS PER SAMPLE");
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding\tFrameshift\tPercentage frameshift\tCoding heterozygous\tPercentage frameshift heterozygous");
		out.print("\tMinor allele known\tWith the minor allele\tWith a unique allele");
		out.println();
		for(int i=0;i<sampleIds.size();i++) {
			VariantsBasicCounts strCountsSample = countsPerSample[2][i];
			out.print(sampleIds.get(i));
			out.print("\t"+strCountsSample.getGenotyped());
			out.print("\t"+strCountsSample.getNonReference());
			out.print("\t"+strCountsSample.getHomozygousAlternative());
			out.print("\t"+strCountsSample.getHeterozygous());
			out.print("\t"+strCountsSample.getCodingTotalCount());
			out.print("\t"+strCountsSample.getTotalCount(Transcriptome.ANNOTATION_FRAMESHIFT));
			out.print("\t"+fmt.format(strCountsSample.getPCTFrameshift()));
			out.print("\t"+strCountsSample.getCodingHeterozygousCount());
			out.print("\t"+strCountsSample.getHeterozygousCount(Transcriptome.ANNOTATION_FRAMESHIFT));
			out.print("\t"+fmt.format(strCountsSample.getPCTFrameshiftHeterozygous()));
			out.print("\t"+strCountsSample.getGenotypedPopCounts());
			out.print("\t"+strCountsSample.getRareAllele());
			out.print("\t"+strCountsSample.getUniqueAllele());
			out.println();
			
		}
		out.println();
		genericPrintCountsPerSample(out,"OTHER BIALLELIC VARIANTS COUNTS PER SAMPLE",3);
		genericPrintCountsPerSample(out,"MULTIALLELIC SNP COUNTS PER SAMPLE",4);
		genericPrintCountsPerSample(out,"MULTIALLELIC INDEL COUNTS PER SAMPLE",5);
		genericPrintCountsPerSample(out,"MULTIALLELIC STR COUNTS PER SAMPLE",6);
		genericPrintCountsPerSample(out,"OTHER MULTIALLELIC VARIANTS COUNTS PER SAMPLE",7);
		
	}
	public void genericPrintCountsPerSample(PrintStream out, String title, int index) {
		out.println(title);
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding");
		out.print("\tMajor allele known\tWith a minor allele");
		out.println();
		for(int i=0;i<sampleIds.size();i++) {
			VariantsBasicCounts countsSample = countsPerSample[index][i];
			out.print(sampleIds.get(i));
			out.print("\t"+countsSample.getGenotyped());
			out.print("\t"+countsSample.getNonReference());
			out.print("\t"+countsSample.getHomozygousAlternative());
			out.print("\t"+countsSample.getHeterozygous());
			out.print("\t"+countsSample.getCodingTotalCount());
			out.print("\t"+countsSample.getGenotypedPopCounts());
			out.print("\t"+countsSample.getRareAllele());
			out.println();
		}
		out.println();
	}
	private double[][] consolidateMAFDistributions() {
		double [] dist = mafDistribution[0].getDistribution();
		double [][] answer = new double [dist.length][6];
		for(int i=0;i<dist.length;i++) answer[i][0] = dist[i];
		dist = mafDistributionAnns.get(Transcriptome.ANNOTATION_SYNONYMOUS).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][1] = dist[i];
		dist = mafDistributionAnns.get(Transcriptome.ANNOTATION_MISSENSE).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][2] = dist[i];
		dist = mafDistributionAnns.get(Transcriptome.ANNOTATION_NONSENSE).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][3] = dist[i];
		dist = mafDistribution[1].getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][4] = dist[i];
		dist = mafDistributionAnns.get(Transcriptome.ANNOTATION_FRAMESHIFT).getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][5] = dist[i];
		return answer;
	}
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	public int getMinSamplesGenotyped() {
		return minSamplesGenotyped;
	}
	public void setMinSamplesGenotyped(int minSamplesGenotyped) {
		this.minSamplesGenotyped = minSamplesGenotyped;
	}
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
}
