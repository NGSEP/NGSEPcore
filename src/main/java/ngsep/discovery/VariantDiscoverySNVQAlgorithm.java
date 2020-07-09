package ngsep.discovery;

import java.util.ArrayList;
import java.util.List;

import ngsep.math.PhredScoreHelper;
import ngsep.sequences.DNASequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.variants.VariantCallReport;

public class VariantDiscoverySNVQAlgorithm {

	private static int posPrint = -1;
	
	//PRE: variant!=null
	public static CalledGenomicVariant genotypeSNV(GenomicVariant variant, CountsHelper countsHelper, double heterozygosityRate, boolean calcStrandBias) {
		CalledGenomicVariant newCall;
		if(countsHelper.getTotalCount()==0) {
			CalledGenomicVariantImpl undecidedCall = new CalledGenomicVariantImpl(variant, new byte [0]);
			undecidedCall.setAllCounts(countsHelper.getCounts());
			return undecidedCall;
		}
		int [] allCounts = countsHelper.getCounts();
		double [] [] allLogConditionals = countsHelper.getLogConditionalProbs();
		double [][] postProbs = countsHelper.getPosteriorProbabilities(heterozygosityRate);
		
		
		if(variant instanceof SNV) {
			SNV snv = (SNV)variant;
			byte indexRef = snv.getRefBaseDNAIndex();
			byte indexAlt = snv.getAltBaseDNAIndex();
			double pHomoRef = postProbs[indexRef][indexRef];
			double pMax = pHomoRef;
			byte genotype = CalledGenomicVariant.GENOTYPE_HOMOREF;
			double pHomoAlt = postProbs[indexAlt][indexAlt];
			if(pHomoAlt>pMax+0.01) {
				pMax = pHomoAlt;
				genotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
			}
			double pHetero = postProbs[indexRef][indexAlt]+postProbs[indexAlt][indexRef];
			if(pHetero>pMax+0.01) {
				pMax = pHetero;
				genotype = CalledGenomicVariant.GENOTYPE_HETERO;
			}
			short gq = PhredScoreHelper.calculatePhredScore(1-pMax);
			if(gq==0) genotype = CalledGenomicVariant.GENOTYPE_UNDECIDED;
			CalledSNV csnv = new CalledSNV((SNV) variant, genotype);
			csnv.setGenotypeQuality(gq);
			csnv.setTotalReadDepth(countsHelper.getTotalCount());
			csnv.setAllBaseCounts(allCounts);
			csnv.setAllGenotypeLogConditionals(allLogConditionals);
			if(calcStrandBias && genotype!=CalledGenomicVariant.GENOTYPE_UNDECIDED && genotype!=CalledGenomicVariant.GENOTYPE_HOMOREF) {
				csnv.setStrandBiasScore(countsHelper.getScoreStrandBiasFisher(indexRef, indexAlt));
			}
			newCall = csnv;
		} else {
			String [] alleles = variant.getAlleles();
			int [] indexes = new int [alleles.length];
			for(int i=0;i<alleles.length;i++) {
				indexes[i] = DNASequence.BASES_STRING.indexOf(alleles[i].charAt(0));
			}
			int [] reportCounts = makeReportCounts(allCounts, indexes);
			double [][] reportLogs = makeReportProbs(allLogConditionals, indexes);
			double [][] reportPosteriors = makeReportProbs(postProbs, indexes);
			int [] indexesMax = getIndexesMaxGenotype(reportPosteriors, 0);
			double maxP = reportPosteriors[indexesMax[0]][indexesMax[1]];
			byte [] genotype;
			if(indexesMax[0]!=indexesMax[1]) {
				maxP+= reportPosteriors[indexesMax[1]][indexesMax[0]];
				genotype = new byte[2];
				genotype[0] = (byte) indexesMax[0];
				genotype[1] = (byte) indexesMax[1];
			} else {
				genotype = new byte[1];
				genotype[0] = (byte) indexesMax[0];
			}
			short gq = PhredScoreHelper.calculatePhredScore(1-maxP);
			//Triallelic variant 
			CalledGenomicVariantImpl call = new CalledGenomicVariantImpl(variant, genotype);
			call.setGenotypeQuality(gq);
			call.setTotalReadDepth(countsHelper.getTotalCount());
			call.setAllCounts(allCounts);
			VariantCallReport report = new VariantCallReport(alleles, reportCounts, reportLogs);
			call.setCallReport(report);
			if(calcStrandBias) {
				if(genotype.length==1 && genotype[0]!=0) call.setStrandBiasScore(countsHelper.getScoreStrandBiasFisher(indexes[0], indexes[genotype[0]]));
				else if (genotype.length==2) call.setStrandBiasScore(countsHelper.getScoreStrandBiasFisher(indexes[genotype[0]], indexes[genotype[1]]));
			}
			newCall = call;
		}
		return newCall;
	}
	//Calls the SNVQ algorithm from the counts stored in the counts helper object
	//PRE: Reference base is uppercase
	public static CalledGenomicVariant discoverSNV(CountsHelper countsHelper, String sequenceName, int position, char refBase, double heterozygosityRate, boolean calcStrandBias) {
		if(countsHelper.getTotalCount()==0) {
			return null;
		}
		if(DNASequence.BASES_STRING.indexOf(refBase)<0) {
			//N reference can in principle be handled but it generates  many non variant sites
			return null;
		}
		String bases = DNASequence.BASES_STRING;
		byte indexRef = (byte) bases.indexOf(refBase);
		int [] counts = countsHelper.getCounts();
		double [][] postProbs = countsHelper.getPosteriorProbabilities(heterozygosityRate);
		if(position==posPrint) {
			System.out.println("Posteriors");
			countsHelper.printProbs(postProbs, false);
		}
		int [] indexesMax = getIndexesMaxGenotype(postProbs, indexRef);
		double refProb = 0;
		int indexI=indexesMax[0],indexJ=indexesMax[1];
		if(indexRef>=0) {
			refProb = postProbs[indexRef][indexRef];
		}
		double maxP = postProbs[indexI][indexJ];
		if(indexI!=indexJ) maxP += postProbs[indexJ][indexI];
		short gq = PhredScoreHelper.calculatePhredScore(1-maxP);
		if(position==posPrint) System.out.println("Max probability: "+maxP+". IndexI: "+indexI+" indexJ: "+indexJ); 
		int indexAlt=0;
		//Solve first non SNV case
		if(indexRef<0 || (indexI!=indexJ && indexI!=indexRef && indexJ!=indexRef)) {
			//Triallelic or weird reference variant
			int indexThird=-1;
			int nAlleles = 2;
			if(indexI!=indexJ) {
				if(postProbs[indexI][indexI]>postProbs[indexJ][indexJ]+0.01) {
					indexAlt = indexI;
					indexThird = indexJ;
				} else {
					indexAlt = indexJ;
					indexThird = indexI;
				}
				nAlleles = 3;
			} else {
				indexAlt = indexI;
			}
			
			List<String> alleles = new ArrayList<String>();
			byte [] idsCalledAlleles = new byte[nAlleles-1];
			idsCalledAlleles[0] = 1;
			int [] indexes = new int [nAlleles];
			indexes[0] = indexRef;
			indexes[1] = indexAlt;
			if(indexRef>=0) {
				alleles.add(DNASequence.BASES_ARRAY[indexRef]);
			} else {
				alleles.add(""+refBase);
			}
			alleles.add(DNASequence.BASES_ARRAY[indexAlt]);
			if(indexThird>=0) {
				alleles.add(DNASequence.BASES_ARRAY[indexThird]);
				idsCalledAlleles[1] = 2;
				indexes[2] = indexThird;
			}
			int[] reportCounts = makeReportCounts(counts, indexes);
			GenomicVariantImpl gv = new GenomicVariantImpl(sequenceName, position, alleles);
			gv.setType(GenomicVariant.TYPE_MULTIALLELIC_SNV);
			gv.setVariantQS(PhredScoreHelper.calculatePhredScore(refProb));
			CalledGenomicVariantImpl triallelicVar = new CalledGenomicVariantImpl(gv, idsCalledAlleles);
			triallelicVar.setGenotypeQuality(gq);
			triallelicVar.setTotalReadDepth(countsHelper.getTotalCount());
			triallelicVar.setAllCounts(counts);
			double [] [] reportLogConds = makeReportProbs(countsHelper.getLogConditionalProbs(),indexes);
			triallelicVar.setCallReport(new VariantCallReport(alleles.toArray(new String [0]), reportCounts, reportLogConds));
			if(calcStrandBias) {
				if(indexThird<0 && indexRef>=0) triallelicVar.setStrandBiasScore(countsHelper.getScoreStrandBiasFisher(indexRef, indexAlt));
				else if (indexThird>=0) triallelicVar.setStrandBiasScore(countsHelper.getScoreStrandBiasFisher(indexAlt, indexThird));
			}
			return triallelicVar;
		}
		char altBase=0;
		byte genotype=0;
		
		if (indexI!=indexJ) {
			if(indexRef!=indexI ) {
				indexAlt = indexI;
			} else {
				indexAlt = indexJ;
			}
			altBase = bases.charAt(indexAlt);
			genotype = CalledGenomicVariant.GENOTYPE_HETERO;
		} else if(indexRef!=indexI ) {
			//Homozygous non reference
			indexAlt = indexI;
			altBase = bases.charAt(indexAlt);
			genotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
		} else {
			//Homozygous reference only useful for genotypeAll mode
			List<String> alleles = new ArrayList<String>();
			alleles.add(DNASequence.BASES_ARRAY[indexRef]);
			GenomicVariantImpl gvNoVar = new GenomicVariantImpl(sequenceName, position, alleles);
			gvNoVar.setVariantQS(PhredScoreHelper.calculatePhredScore(refProb));
			byte [] idsCalledAlleles = {0};
			int [] indexes = {indexRef};
			CalledGenomicVariantImpl noVarCall = new CalledGenomicVariantImpl(gvNoVar, idsCalledAlleles);
			noVarCall.setGenotypeQuality(gq);
			noVarCall.setTotalReadDepth(countsHelper.getTotalCount());
			noVarCall.setAllCounts(counts);
			int[] reportCounts = makeReportCounts(counts, indexes);
			double [] [] reportLogConds = makeReportProbs(countsHelper.getLogConditionalProbs(),indexes);
			noVarCall.setCallReport(new VariantCallReport(alleles.toArray(new String [0]), reportCounts, reportLogConds));
			return noVarCall;
		}
		SNV snv = new SNV(sequenceName, position, refBase, altBase);
		snv.setVariantQS(PhredScoreHelper.calculatePhredScore(refProb));
		CalledSNV csnv = new CalledSNV(snv,genotype);
		csnv.setGenotypeQuality(gq);
		csnv.setTotalReadDepth(countsHelper.getTotalCount());
		csnv.setAllBaseCounts(counts);
		csnv.setAllGenotypeLogConditionals(countsHelper.getLogConditionalProbs());
		if(calcStrandBias && genotype!=CalledGenomicVariant.GENOTYPE_UNDECIDED && genotype!=CalledGenomicVariant.GENOTYPE_HOMOREF) {
			csnv.setStrandBiasScore(countsHelper.getScoreStrandBiasFisher(indexRef, indexAlt));
		}
		return csnv;
	}
	private static int [] getIndexesMaxGenotype (double [][] genotypePosteriors, int indexDefault) {
		if(indexDefault<0 || indexDefault>=genotypePosteriors.length) {
			indexDefault=0;
		}
		int [] indexes = {indexDefault,indexDefault};
		double probMax = genotypePosteriors[indexDefault][indexDefault];
		for(int i=0;i<genotypePosteriors.length;i++) {
			for(int j=i;j<genotypePosteriors[0].length;j++) {
				//For heterozygous genotypes add up the two events corresponding with the two ways to sort the alleles
				double genotypeProbability = genotypePosteriors[i][j];
				if(i!=j) genotypeProbability += genotypePosteriors[j][i];
				//Differences of less than 0.01 are considered equal
				if (genotypeProbability > probMax+0.01) {
					probMax = genotypeProbability;
					indexes[0] = i;
					indexes[1] = j;
				}
			}
		}
		return indexes;
	}
	private static int[] makeReportCounts(int[] counts, int [] indexes) {
		int [] reportCounts = new int [indexes.length];
		for(int i=0;i<indexes.length;i++) {
			if(indexes[i] <0 || indexes[i] >= counts.length) reportCounts[i] = 0;
			else reportCounts[i] = counts[indexes[i]];
		}
		return reportCounts;
	}
	private static double[][] makeReportProbs(double[][] allValues, int [] indexes) {
		double [][] answer = new double[indexes.length][indexes.length];
		for(int i=0;i<indexes.length;i++) {
			for(int j=0;j<indexes.length;j++) {
				if(indexes[i] <0 || indexes[i] >= allValues.length) answer[i][j] = 0;
				else if(indexes[j] <0 || indexes[j] >= allValues[0].length) answer[i][j] = 0;
				else answer[i][j] = allValues[indexes[i]][indexes[j]];
			}
		}
		return answer;
	}

	
	public static CalledGenomicVariant callIndel (PileupRecord pileup, CountsHelper helper, GenomicVariant variant, double heterozygosityRate, boolean calcStrandBias) {
		int [] counts = helper.getCounts();
		double [][] postProbs = helper.getPosteriorProbabilities(heterozygosityRate);
		if(pileup.getPosition()==posPrint) {
			System.out.println("Counts");
			for(int j=0;j<counts.length;j++) System.out.println("Count allele: "+helper.getAlleles()[j]+": "+counts[j]); 
			System.out.println("Probs");
			helper.printProbs(postProbs, false);
		}
		if(helper.getTotalCount()==0) {
			if(variant == null) return null;
			return new CalledGenomicVariantImpl(variant, new byte [0]);
		}
		int [] indexesMax = getIndexesMaxGenotype(postProbs, 0);
		if(pileup.getPosition()==posPrint) System.out.println("Indexes max: "+indexesMax[0]+" "+indexesMax[1]);
		byte [] calledAlleles;
		GenomicVariant gv = variant;
		int[] reportCounts = counts;
		double [][] reportLogs = helper.getLogConditionalProbs();
		if(gv == null) {
			String [] helperAlleles = helper.getAlleles();
			List<String> alleles = new ArrayList<String>();
			List<Integer> indexesList = new ArrayList<Integer>();
			//Add reference allele
			alleles.add(helperAlleles[0]);
			int referenceLength = helperAlleles[0].length();
			indexesList.add(0);
			boolean lengthChange = false;
			if(indexesMax[0] > 0 && indexesMax[0] < helperAlleles.length) {
				String allele0 = helperAlleles[indexesMax[0]]; 
				alleles.add(allele0);
				indexesList.add(indexesMax[0]);
				if(allele0.length()!=referenceLength) lengthChange = true;
			}
			if(indexesMax[1] > 0 && indexesMax[1] !=indexesMax [0] && indexesMax[1]<helperAlleles.length ) {
				String allele1 = helperAlleles[indexesMax[1]];
				alleles.add(allele1);
				indexesList.add(indexesMax[1]);
				if(allele1.length()!=referenceLength) lengthChange = true;
				if(alleles.size()==3 && allele1.length()!=alleles.get(1).length()) lengthChange = true;
			}
			//No real indel if all alleles have the same length
			if(!lengthChange && !pileup.isInputSTR()) return null;
			GenomicVariantImpl newVar = new GenomicVariantImpl(pileup.getSequenceName(), pileup.getPosition(), alleles);
			if(pileup.isSTR()) newVar.setType(GenomicVariant.TYPE_STR);
			else newVar.setType(GenomicVariant.TYPE_INDEL);
			newVar.setVariantQS(PhredScoreHelper.calculatePhredScore(postProbs[0][0]));
			gv = newVar;
			int [] indexes = new int [indexesList.size()];
			for(int i=0;i<indexes.length;i++) indexes[i] = indexesList.get(i);
			reportCounts = makeReportCounts(counts, indexes);
			reportLogs = makeReportProbs(helper.getLogConditionalProbs(), indexes);
			if(indexesMax[1] !=indexesMax [0]) {
				calledAlleles = new byte[2];
				if (alleles.size()==3) {
					calledAlleles[0] = 1;
					calledAlleles[1] = 2;
				} else  {
					calledAlleles[0] = 0;
					calledAlleles[1] = 1;
				}
			} else {
				calledAlleles = new byte[1];
				if(indexesMax[0] == 0) {
					calledAlleles[0] = 0;
				} else {
					calledAlleles[0] = 1;
				}
			}
		} else {
			if (indexesMax[0] > GenomicVariant.MAX_NUM_ALLELES || indexesMax[1]>GenomicVariant.MAX_NUM_ALLELES) {
				calledAlleles = new byte[0];
			} else if(indexesMax[1] !=indexesMax [0]) {
				calledAlleles = new byte[2];
				calledAlleles[0] = (byte) indexesMax[0];
				calledAlleles[1] = (byte) indexesMax[1];
			} else {
				calledAlleles = new byte[1];
				calledAlleles[0] = (byte) indexesMax[0];
			}
		}
		
		int totalDepth = helper.getTotalCount();
		VariantCallReport report = new VariantCallReport(gv.getAlleles(), reportCounts, reportLogs);
		CalledGenomicVariantImpl newCall = new CalledGenomicVariantImpl(gv, calledAlleles);
		double maxP = postProbs[indexesMax[0]][indexesMax[1]];
		if(indexesMax[0]!=indexesMax[1]) maxP+= postProbs[indexesMax[1]][indexesMax[0]];
		newCall.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(1-maxP));
		newCall.setTotalReadDepth(totalDepth);
		if(totalDepth>0) newCall.setCallReport(report);
		if(calcStrandBias) {
			if(calledAlleles.length==1 && calledAlleles[0]!=0) newCall.setStrandBiasScore(helper.getScoreStrandBiasFisher(0, calledAlleles[0]));
			else if (calledAlleles.length==2) newCall.setStrandBiasScore(helper.getScoreStrandBiasFisher(calledAlleles[0], calledAlleles[1]));
		}
		//if(pileup.getFirst()==82) System.out.println("Indel alleles: "+newCall.getAlleles().length+" called alleles: "+calledAlleles[0]+" "+calledAlleles[1]+" genotype prob: "+newCall.getGenotypeProbability());
		return newCall;
	}
}
