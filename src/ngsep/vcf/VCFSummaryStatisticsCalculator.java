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
import java.util.Iterator;
import java.util.List;
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

	public static final int VARIANT_CATEGORIES = 3;
	private int minSamplesGenotyped = 20;
	private VariantsBasicCounts [] summaryCounts = new VariantsBasicCounts[VARIANT_CATEGORIES];
	private List<String> sampleIds;
	private VariantsBasicCounts [][] countsPerSample;
	private int [] totalGenotypeCalls = new int [VARIANT_CATEGORIES];
	private Distribution [] genotypedAccessionsDistribution = new Distribution[VARIANT_CATEGORIES];
	private Distribution [] mafDistribution = new Distribution[VARIANT_CATEGORIES];
	private Distribution mafDistributionSynonymous = new Distribution(-0.05, 0.51, 0.1);
	private Distribution mafDistributionMissense = new Distribution(-0.05, 0.51, 0.1);
	private Distribution mafDistributionNonsense = new Distribution(-0.05, 0.51, 0.1);
	private Distribution  mafDistributionFrameshift = new Distribution(-0.05, 0.51, 0.1);
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
		countsPerSample = new VariantsBasicCounts[VARIANT_CATEGORIES][sampleIds.size()];
		for(int i=0;i<VARIANT_CATEGORIES;i++) {
			summaryCounts[i] = new VariantsBasicCounts();
			totalGenotypeCalls[i] = 0;
			genotypedAccessionsDistribution[i] = new Distribution(0, sampleIds.size(), 1);
			mafDistribution[i] = new Distribution(-0.05, 0.51, 0.1);
			
			for(int j=0;j<sampleIds.size();j++) {
				countsPerSample[i][j] = new VariantsBasicCounts();
			}
		}
		
	}
	public void processRecord(VCFRecord record) {
		GenomicVariant var = record.getVariant(); 
		String [] alleles = var.getAlleles();
		int idxVarType = 2;
		if(record.getVariant() instanceof SNV) {
			idxVarType = 0;
		} else if(alleles.length==2) {
			idxVarType = 1;
		}
		List<CalledGenomicVariant> varCalls = record.getCalls();
		if(varCalls.size()!=countsPerSample[0].length) {
			throw new IllegalArgumentException("Inconsistent number of calls for variant "+record.getVariant().getSequenceName()+":"+record.getVariant().getFirst());
		}
		boolean isBiallelicSNV = (idxVarType == 0);
		boolean isTransition = isBiallelicSNV && ((SNV)var).isTransition();
		int functionalStatus = calculateFunctionalStatus(record.getInfoField(GenomicVariantAnnotation.ATTRIBUTE_TRANSCRIPT_ANNOTATION));
		DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(varCalls, false);
		genotypedAccessionsDistribution[idxVarType].processDatapoint(divStats.getNumSamplesGenotyped());
		totalGenotypeCalls[idxVarType] += divStats.getNumSamplesGenotyped();
		double maf = divStats.getMaf();
		byte mafIdx = divStats.getMafIndex();
		byte wtIdx = divStats.getWtIndex();
		int [] alleleCounts = divStats.getAlleleCounts();
		int populationStatus = 0;
		if(divStats.getNumSamplesGenotyped()>=minSamplesGenotyped) {
			populationStatus = VariantsBasicCounts.POPULATION_STATUS_GENPOP;
			mafDistribution[idxVarType].processDatapoint(maf);
			if(functionalStatus==VariantsBasicCounts.FUNCTIONAL_STATUS_SYNONYMOUS) {
				mafDistributionSynonymous.processDatapoint(maf);
			} else if(functionalStatus==VariantsBasicCounts.FUNCTIONAL_STATUS_MISSENSE) {
				mafDistributionMissense.processDatapoint(maf);
			} else if(functionalStatus==VariantsBasicCounts.FUNCTIONAL_STATUS_NONSENSE) {
				mafDistributionNonsense.processDatapoint(maf);
			}  else if(functionalStatus==VariantsBasicCounts.FUNCTIONAL_STATUS_FRAMESHIFT) {
				mafDistributionFrameshift.processDatapoint(maf);
			}
		}		
		summaryCounts[idxVarType].processGenotypeCall(VariantsBasicCounts.GENOTYPE_STATUS_HOMOALT, isTransition, functionalStatus, populationStatus);
		for(int i=0;i<varCalls.size();i++) {
			CalledGenomicVariant call = varCalls.get(i);
			int genotypeStatus = calculateGenotypeStatus (call);
			int popStatusSample = populationStatus;
			if(popStatusSample== VariantsBasicCounts.POPULATION_STATUS_GENPOP && wtIdx>=0) {
				byte [] calledAlleles = call.getIndexesCalledAlleles();
				byte [] allelesCN = call.getAllelesCopyNumber();
				
				for(int j=0;j<calledAlleles.length;j++) {
					if(calledAlleles[j] != wtIdx) {
						popStatusSample = VariantsBasicCounts.POPULATION_STATUS_RARE;
						if(mafIdx >=0 && allelesCN[mafIdx] == alleleCounts[mafIdx]) {
							popStatusSample = VariantsBasicCounts.POPULATION_STATUS_UNIQUE;
						}
						break;
					}
				}
				//if(isBiallelicSNV && calledAlleles.length>0)System.out.println("MAF: "+maf+". MAF idx: "+mafIdx+" AC: "+alleleCounts[0]+" - "+ alleleCounts[1]+" called: "+calledAlleles[0]+" hetero: "+ (calledAlleles.length>1)+ " ploidy a1: "+calledAllelePloidies[0]+" status: "+populationStatus);
			}
			countsPerSample[idxVarType][i].processGenotypeCall(genotypeStatus, isTransition, functionalStatus, popStatusSample);	
		}
		
	}

	private int calculateGenotypeStatus(CalledGenomicVariant call) {
		if(call.isUndecided()) return VariantsBasicCounts.GENOTYPE_STATUS_UNDECIDED;
		if(call.isHomozygousReference()) return VariantsBasicCounts.GENOTYPE_STATUS_HOMOREF;
		if(call.isHeterozygous()) return VariantsBasicCounts.GENOTYPE_STATUS_HETEROZYGOUS;
		return VariantsBasicCounts.GENOTYPE_STATUS_HOMOALT;
	}

	private int calculateFunctionalStatus(GenomicVariantAnnotation annotation) {
		if(annotation!=null) {
			if (Transcriptome.ANNOTATION_CODING.equals(annotation.getValue())) {
				return VariantsBasicCounts.FUNCTIONAL_STATUS_GENERALCODING;
			} else if (Transcriptome.ANNOTATION_SYNONYMOUS.equals(annotation.getValue())) {
				return VariantsBasicCounts.FUNCTIONAL_STATUS_SYNONYMOUS;
			} else if (Transcriptome.ANNOTATION_MISSENSE.equals(annotation.getValue())) {
				return VariantsBasicCounts.FUNCTIONAL_STATUS_MISSENSE;
			} else if (Transcriptome.ANNOTATION_NONSENSE.equals(annotation.getValue())) {
				return VariantsBasicCounts.FUNCTIONAL_STATUS_NONSENSE;
			}  else if (Transcriptome.ANNOTATION_FRAMESHIFT.equals(annotation.getValue())) {
				return VariantsBasicCounts.FUNCTIONAL_STATUS_FRAMESHIFT;
			}
		}
		return 0;
	}
	
	private void printStatistics(PrintStream out) {
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
		out.println("SUMMARY");
		out.println("Count\tBiallelic SNPs\tBiallelic indels\tOther");
		out.println("Variants\t"+summaryCounts[0].getGenotyped()+"\t"+summaryCounts[1].getGenotyped()+"\t"+summaryCounts[2].getGenotyped());
		out.println("Genotype calls\t"+totalGenotypeCalls[0]+"\t"+totalGenotypeCalls[1]+"\t"+totalGenotypeCalls[2]);
		out.println("Coding variants\t"+summaryCounts[0].getCoding()+"\t"+summaryCounts[1].getCoding()+"\t"+summaryCounts[2].getCoding());
		out.println("Variants with known minor allele\t"+summaryCounts[0].getGenotypedPopCounts()+"\t"+summaryCounts[1].getGenotypedPopCounts()+"\t"+summaryCounts[2].getGenotypedPopCounts());
		out.println();
		out.println("Synonymous:\t"+summaryCounts[0].getSynonymous());
		out.println("Missense:\t"+summaryCounts[0].getMissense());
		out.println("Non sense:\t"+summaryCounts[0].getNonsense());
		out.println("NS/S rate:\t"+fmt.format(summaryCounts[0].getNonSynonymousToSynonymousRate()));
		out.println("Frameshift indels:\t"+summaryCounts[1].getFrameshift()+"\tPCT:\t"+fmt.format(summaryCounts[1].getPCTFrameshift()));
		out.println("Transitions:\t"+summaryCounts[0].getTransitions()+"\tTr/Tv ratio:\t"+fmt.format(summaryCounts[0].getTrTvRatio()));
		out.println("Synonymous transitions:\t"+summaryCounts[0].getSynonymousTransitions()+"\tSynonymous Tr/Tv ratio:\t"+fmt.format(summaryCounts[0].getTrTvRatioSynonymous()));
		out.println("Non Synonymous transitions:\t"+summaryCounts[0].getNonSynonymousTransitions()+"\tNon synonymous Tr/Tv ratio:\t"+fmt.format(summaryCounts[0].getTrTvRatioNonSynonymous()));
		out.println();
		out.println("MAF DISTRIBUTIONS");
		out.println("MAF\tBiallelic SNPs\tSynonymous\tMissense\tNonsense\tBiallelic indels\tFrameshift");
		double [] [] consolidatedDist = consolidateMAFDistributions();
		for(int i=0;i<consolidatedDist.length;i++) {
			double minMaf = 0;
			if(i>0) minMaf = 0.1*i-0.05;
			double maxMaf = 0.5;
			if(i<consolidatedDist.length-1) maxMaf = 0.1*i+0.05;
			out.print(fmt.format(minMaf)+"-"+fmt.format(maxMaf));
			for(int j=0;j<consolidatedDist[i].length;j++) {
				out.print("\t"+fmt.format(consolidatedDist[i][j]));
			}
			out.println();
		}
		out.println();
		out.println("SAMPLES GENOTYPED DISTRIBUTIONS");
		out.println("Samples genotyped\tBiallelic SNPs\tBiallelic indels\tOther");
		double [] covDistSNPs = genotypedAccessionsDistribution[0].getDistribution();
		double [] covDistIndels = genotypedAccessionsDistribution[1].getDistribution();
		double [] covDistOthers = genotypedAccessionsDistribution[2].getDistribution();
		for(int i=0;i<covDistSNPs.length;i++) {
			int countSNPs = (int) Math.round(covDistSNPs[i]);
			int countIndels = (int) Math.round(covDistIndels[i]);
			int countOther = (int) Math.round(covDistOthers[i]);
			out.println(""+i+"\t"+countSNPs+"\t"+countIndels+"\t"+countOther);
		}
		out.println();
		out.println("SNP COUNTS PER SAMPLE");
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding\tSynonymous\tMissense\tNonsense\tNS/S rate");
		out.print("\tHeterozygous coding\tHeterozygous synonymous\tHeterozygous missense\tHeterozygous nonsense\tNS/S rate Heterozygous");
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
			out.print("\t"+snpCountsSample.getCoding());
			out.print("\t"+snpCountsSample.getSynonymous());
			out.print("\t"+snpCountsSample.getMissense());
			out.print("\t"+snpCountsSample.getNonsense());
			out.print("\t"+fmt.format(snpCountsSample.getNonSynonymousToSynonymousRate()));
			out.print("\t"+snpCountsSample.getCodingHeterozygous());
			out.print("\t"+snpCountsSample.getSynonymousHeterozygous());
			out.print("\t"+snpCountsSample.getMissenseHeterozygous());
			out.print("\t"+snpCountsSample.getNonsenseHeterozygous());
			out.print("\t"+fmt.format(snpCountsSample.getNonSynonymousToSynonymousRateHeterozygous()));
			out.print("\t"+snpCountsSample.getGenotypedPopCounts());
			out.print("\t"+snpCountsSample.getRareAllele());
			out.print("\t"+snpCountsSample.getUniqueAllele());
			out.print("\t"+snpCountsSample.getTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatio()));
			out.print("\t"+snpCountsSample.getHomozygousAlternativeTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioHomozygousAlternative()));
			out.print("\t"+snpCountsSample.getHeterozygousTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioHeterozygous()));
			out.print("\t"+snpCountsSample.getSynonymousTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioSynonymous()));
			out.print("\t"+snpCountsSample.getNonSynonymousTransitions());
			out.print("\t"+fmt.format(snpCountsSample.getTrTvRatioNonSynonymous()));
			out.println();
		}
		out.println();
		out.println("INDEL COUNTS PER SAMPLE");
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
			out.print("\t"+indelCountsSample.getCoding());
			out.print("\t"+indelCountsSample.getFrameshift());
			out.print("\t"+fmt.format(indelCountsSample.getPCTFrameshift()));
			out.print("\t"+indelCountsSample.getCodingHeterozygous());
			out.print("\t"+indelCountsSample.getFrameshiftHeterozygous());
			out.print("\t"+fmt.format(indelCountsSample.getPCTFrameshiftHeterozygous()));
			out.print("\t"+indelCountsSample.getGenotypedPopCounts());
			out.print("\t"+indelCountsSample.getRareAllele());
			out.print("\t"+indelCountsSample.getUniqueAllele());
			out.println();
			
		}
		out.println();
		out.println("OTHER VARIANTS COUNTS PER SAMPLE");
		out.print("Sample\tGenotyped\tNon reference\tHomozygous alternative\tHeterozygous");
		out.print("\tCoding");
		out.print("\tMajor allele known\tWith a minor allele");
		out.println();
		for(int i=0;i<sampleIds.size();i++) {
			VariantsBasicCounts otherCountsSample = countsPerSample[2][i];
			out.print(sampleIds.get(i));
			out.print("\t"+otherCountsSample.getGenotyped());
			out.print("\t"+otherCountsSample.getNonReference());
			out.print("\t"+otherCountsSample.getHomozygousAlternative());
			out.print("\t"+otherCountsSample.getHeterozygous());
			out.print("\t"+otherCountsSample.getCoding());
			out.print("\t"+otherCountsSample.getGenotypedPopCounts());
			out.print("\t"+otherCountsSample.getRareAllele());
			out.println();
		}
	}
	private double[][] consolidateMAFDistributions() {
		double [] dist = mafDistribution[0].getDistribution();
		double [][] answer = new double [dist.length][6];
		for(int i=0;i<dist.length;i++) answer[i][0] = dist[i];
		dist = mafDistributionSynonymous.getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][1] = dist[i];
		dist = mafDistributionMissense.getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][2] = dist[i];
		dist = mafDistributionNonsense.getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][3] = dist[i];
		dist = mafDistribution[1].getDistribution();
		for(int i=0;i<dist.length;i++) answer[i][4] = dist[i];
		dist = mafDistributionFrameshift.getDistribution();
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
