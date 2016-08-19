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
package ngsep.discovery;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.math.CountsRankHelper;
import ngsep.math.PhredScoreHelper;
import ngsep.sequences.DNASequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.variants.VariantCallReport;


public class VariantPileupListener implements PileupListener {
	
	private List<CalledGenomicVariant> calledVariants = new ArrayList<CalledGenomicVariant>();
	private GenomicRegionSortedCollection<GenomicVariant> inputVariants = new GenomicRegionSortedCollection<GenomicVariant>();
	private ReferenceGenome genome;
	private byte normalPloidy = GenomicVariant.DEFAULT_PLOIDY;
	public static final double DEF_HETEROZYGOSITY_RATE_DIPLOID = 0.001;
	public static final double DEF_HETEROZYGOSITY_RATE_HAPLOID = 0.000001;
	
	private double heterozygosityRate = DEF_HETEROZYGOSITY_RATE_DIPLOID;
	private short maxBaseQS=0; 
	private boolean ignoreLowerCaseRef = false;
	private boolean callEmbeddedSNVs = false;
	private boolean genotypeAll = false;
	
	//Filters applied on variants discovery
	public static final int DEF_MIN_ALT_COVERAGE = 0;
	public static final int DEF_MAX_ALT_COVERAGE = 0;
	public static final short DEF_MIN_QUALITY = 0;
	private int minAltCoverage=DEF_MIN_ALT_COVERAGE;
	private int maxAltCoverage=DEF_MAX_ALT_COVERAGE;
	private short minQuality = DEF_MIN_QUALITY;
	
	
	
	
	//Control attribute to avoid calling overlapping indels and to give an embedded status to SNVs within indels or STRs
	private int lastIndelEnd = 0;
	
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	
	public byte getNormalPloidy() {
		return normalPloidy;
	}
	public void setNormalPloidy(byte normalPloidy) {
		this.normalPloidy = normalPloidy;
	}
	public double getHeterozygosityRate() {
		return heterozygosityRate;
	}
	public void setHeterozygosityRate(double heterozygosityRate) {
		this.heterozygosityRate = heterozygosityRate;
	}
	public short getMaxBaseQS() {
		return maxBaseQS;
	}
	public void setMaxBaseQS(short maxBaseQS) {
		this.maxBaseQS = maxBaseQS;
	}
	public boolean isIgnoreLowerCaseRef() {
		return ignoreLowerCaseRef;
	}
	public void setIgnoreLowerCaseRef(boolean ignoreLowerCaseRef) {
		this.ignoreLowerCaseRef = ignoreLowerCaseRef;
	}
	public boolean isCallEmbeddedSNVs() {
		return callEmbeddedSNVs;
	}
	public void setCallEmbeddedSNVs(boolean callEmbeddedSNVs) {
		this.callEmbeddedSNVs = callEmbeddedSNVs;
	}
	public boolean isGenotypeAll() {
		return genotypeAll;
	}
	public void setGenotypeAll(boolean genotypeAll) {
		this.genotypeAll = genotypeAll;
	}
	public int getMinAltCoverage() {
		return minAltCoverage;
	}
	public void setMinAltCoverage(int minAltCoverage) {
		this.minAltCoverage = minAltCoverage;
	}
	public int getMaxAltCoverage() {
		return maxAltCoverage;
	}
	public void setMaxAltCoverage(int maxAltCoverage) {
		this.maxAltCoverage = maxAltCoverage;
	}
	public short getMinQuality() {
		return minQuality;
	}
	public void setMinQuality(short minQuality) {
		this.minQuality = minQuality;
	}
	public List<GenomicVariant> getInputVariants() {
		return inputVariants.asList();
	}
	public void setInputVariants(Collection<GenomicVariant> inputVariants) {
		this.inputVariants = new GenomicRegionSortedCollection<GenomicVariant>();
		if(inputVariants!=null) {
			this.inputVariants.addAll(inputVariants);
		}
	}
	public List<CalledGenomicVariant> getCalledVariants() {
		return calledVariants;
	}
	
	/**
	 * 
	 * @param pileup
	 * @param variant Variant to genotype. If null, the method tries to find a new variant given the genome and the pileup
	 * @return CalledGenomicVariant. New call of the given variant or new called variant if exists
	 */
	public CalledGenomicVariant processPileup(PileupRecord pileup, GenomicVariant variant) {
		String referenceAllele;
		if(variant!=null) {
			referenceAllele = variant.getReference();
		} else {
			if(!callEmbeddedSNVs && pileup.isEmbedded()) return null;
			int last = pileup.getPosition()+pileup.getReferenceSpan()-1;
			CharSequence seq = genome.getReference(pileup.getSequenceName(), pileup.getPosition(), last);
			if(seq == null) return null;
			referenceAllele = seq.toString();
			if(ignoreLowerCaseRef && Character.isLowerCase(referenceAllele.charAt(0))) return null;
			referenceAllele = referenceAllele.toUpperCase();
			//Avoid trying to call nested indels or SNVs within indels unless explicitly requested
			if(lastIndelEnd>=pileup.getPosition()) {
				if(!callEmbeddedSNVs) return null;
				referenceAllele = referenceAllele.substring(0,1);
				if (pileup.isSTR()) {
					pileup.setSTR(false);
				}
			}
		}
		CountsHelper helperSNV = calculateCountsSNV(pileup);
		//if(pileup.getFirst()==82) System.out.println("Pileup last: "+pileup.getLast()+" Reference allele: "+referenceAllele); 
		//if(pileup.getPosition()==9052) System.out.println("Reference allele: "+referenceAllele+". Pileup last: "+pileup.getLast());
		CalledGenomicVariant calledVar;
		if(referenceAllele.length()>1) {
			CountsHelper helperIndel = calculateCountsIndel(pileup,variant,referenceAllele); 
			calledVar = callIndel(pileup, helperIndel, variant);
			if(variant == null) {
				if(calledVar!=null && (pileup.isInputSTR() || (!calledVar.isUndecided() && !calledVar.isHomozygousReference()))) {
					//System.out.println("Called indel at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+" variant type: "+calledVar.getType());
					lastIndelEnd = calledVar.getLast();
				} else {
					if (pileup.isNewSTR()) {
						pileup.setSTR(false);
						pileup.setNewSTR(false);
					}
					//Try SNV if the indel alleles were not good to make a call
					calledVar = callSNV(pileup, helperSNV, variant, referenceAllele.charAt(0));
				}
			}
		} else {
			calledVar = callSNV(pileup, helperSNV, variant, referenceAllele.charAt(0));
		}
		if(calledVar != null && (variant!=null || genotypeAll || (!calledVar.isUndecided() && !calledVar.isHomozygousReference()) )) {
			//System.out.println("Called SNV");
			calledVar.updateAllelesCopyNumberFromCounts(normalPloidy);
			if(calledVar.isSNV() && (pileup.isEmbedded() || calledVar.getFirst()<=lastIndelEnd)) calledVar.setType(GenomicVariant.TYPE_EMBEDDED_SNV);
			//System.out.println("Called variant at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+" pileup variant type: "+pileup.getVariantType()+" variant type: "+calledVar.getType());
			//if(pileup.getPosition()==3892) System.out.println("New variant at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+". Repeat context: "+calledVar.getRepeatContext());
			return calledVar;
		}
		return null;
	}
	
	private CountsHelper calculateCountsSNV (PileupRecord pileup) {
		CountsHelper answer = new CountsHelper();
		if(maxBaseQS>0) answer.setMaxBaseQS(maxBaseQS);
		List<String> callsWithScores = pileup.getAlleleCalls(1);
		for(int i=0;i<callsWithScores.size();i+=2) {
			String call = callsWithScores.get(i);
			String scores = callsWithScores.get(i+1);
			answer.updateCounts(call.substring(0,1), scores.charAt(0), (short)255);
		}
		return answer;
	}
	//PRE: Reference allele is not null and its length is larger than 1. If variant is not null, reference allele is the reference allele of the variant
	private CountsHelper calculateCountsIndel(PileupRecord pileup, GenomicVariant variant, String referenceAllele) {
		int posPrint = -1;
		String [] indelAlleles;
		List<String> callsWithScores = pileup.getAlleleCalls(referenceAllele.length());
		if(variant!=null) {
			indelAlleles = variant.getAlleles();
		} else {
			CountsRankHelper<String> countsRank = new CountsRankHelper<String>();
			for(int i=0;i<callsWithScores.size();i+=2) {
				String allele = callsWithScores.get(i);
				countsRank.add(allele);	
			}
			if(countsRank.getNumDifferent()>100) System.err.println("WARN: Number of alleles for site at "+pileup.getSequenceName()+":"+pileup.getPosition()+" is "+countsRank.getNumDifferent()+" ref Allele: "+referenceAllele);
			Set<String> pileupAlleles = countsRank.selectBest(10);
			pileupAlleles.add(referenceAllele);
			indelAlleles = new String [pileupAlleles.size()];
			indelAlleles[0] = referenceAllele;
			if(pileup.getPosition()==posPrint) System.out.println("Reference allele for indel: "+referenceAllele);
			int i=1;
			for (String allele:pileupAlleles) {
				if(!allele.equals(referenceAllele)) {
					indelAlleles[i] = allele;
					if(pileup.getPosition()==posPrint) System.out.println("Next alternative allele for indel: "+allele);
					i++;
				}
			}
		}
		CountsHelper answer = new CountsHelper(indelAlleles);
		if(maxBaseQS>0) answer.setMaxBaseQS(maxBaseQS);
		
		for(int i=0;i<callsWithScores.size();i+=2) {
			String call = callsWithScores.get(i);
			//TODO: Make better score. By now a fixed error probability of 0.01 is used
			String scores = callsWithScores.get(i+1);
			if(pileup.getPosition()==posPrint) System.out.println("Next call for indel: "+call+" scores: "+scores);
			answer.updateCounts(call, '5', (short)255);
		}
		return answer;
	}
	//PRE: Reference base is uppercase
	public CalledGenomicVariant callSNV(PileupRecord pileup, CountsHelper countsHelper, GenomicVariant variant, char refBase) {
		CalledGenomicVariant newCall;
		int depth = countsHelper.getTotalCount();
		if(variant!=null) {
			int [] allCounts = countsHelper.getCounts();
			double [] [] allLogConditionals = countsHelper.getLogConditionalProbs();
			double [][] postProbs = countsHelper.getPosteriorProbabilities(heterozygosityRate);
			
			
			if(variant instanceof SNV) {
				byte indexRef = (byte) DNASequence.BASES_STRING.indexOf(refBase);
				int [] indexesMax = getIndexesMaxGenotype(postProbs, indexRef);
				double maxP = postProbs[indexesMax[0]][indexesMax[1]];
				if(indexesMax[0]!=indexesMax[1]) maxP+= postProbs[indexesMax[1]][indexesMax[0]];
				short gq = PhredScoreHelper.calculatePhredScore(1-maxP);
				byte genotype = 0;
				if(indexesMax[0]!=indexesMax[1]) genotype = 1;
				else if (indexesMax[0]!=indexRef) genotype = 2;
				CalledSNV csnv = new CalledSNV((SNV) variant, genotype);
				csnv.setGenotypeQuality(gq);
				csnv.setTotalReadDepth(depth);
				csnv.setAllBaseCounts(allCounts);
				csnv.setAllGenotypeLogConditionals(allLogConditionals);
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
				if(indexesMax[0]!=indexesMax[1]) maxP+= reportPosteriors[indexesMax[1]][indexesMax[0]];
				short gq = PhredScoreHelper.calculatePhredScore(1-maxP);
				//Triallelic variant
				byte [] indexesAsBytes = {(byte) indexesMax[0],(byte) indexesMax[1]}; 
				CalledGenomicVariantImpl call = new CalledGenomicVariantImpl(variant, indexesAsBytes);
				call.setGenotypeQuality(gq);
				call.setTotalReadDepth(depth);
				VariantCallReport report = new VariantCallReport(alleles, reportCounts, reportLogs);
				call.setCallReport(report);
				newCall = call;
			}
		} else {
			if(countsHelper.getTotalCount()==0) {
				return null;
			}
			if(DNASequence.BASES_STRING.indexOf(refBase)<0) {
				//N reference can in principle be handled but it generates  many non variant sites
				return null;
			}
			newCall = callSNVQ(pileup,countsHelper,refBase);
		}
		if(!passFilter(newCall)) {
			newCall.makeUndecided();
		}
		return newCall;
	}
	//Calls the SNVQ algorithm from the counts stored in the counts helper object
	private CalledGenomicVariant callSNVQ(PileupRecord pileup, CountsHelper countsHelper, char refBase) {
		int posPrint = -1;
		String bases = DNASequence.BASES_STRING;
		byte indexRef = (byte) bases.indexOf(refBase);
		int [] counts = countsHelper.getCounts();
		double [][] postProbs = countsHelper.getPosteriorProbabilities(heterozygosityRate);
		if(pileup.getPosition()==posPrint) {
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
		if(pileup.getPosition()==posPrint) System.out.println("Max probability: "+maxP+". IndexI: "+indexI+" indexJ: "+indexJ); 
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
			GenomicVariantImpl gv = new GenomicVariantImpl(pileup.getSequenceName(), pileup.getPosition(), alleles);
			gv.setType(GenomicVariant.TYPE_MULTIALLELIC_SNV);
			gv.setVariantQS(PhredScoreHelper.calculatePhredScore(refProb));
			CalledGenomicVariantImpl triallelicVar = new CalledGenomicVariantImpl(gv, idsCalledAlleles);
			triallelicVar.setGenotypeQuality(gq);
			triallelicVar.setTotalReadDepth(countsHelper.getTotalCount());
			triallelicVar.setAllCounts(counts);
			double [] [] reportLogConds = makeReportProbs(countsHelper.getLogConditionalProbs(),indexes);
			triallelicVar.setCallReport(new VariantCallReport(alleles.toArray(new String [0]), reportCounts, reportLogConds));
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
			genotype = 1;
		} else if(indexRef!=indexI ) {
			//Homozygous non reference
			indexAlt = indexI;
			altBase = bases.charAt(indexAlt);
			genotype = 2;
		} else {
			//Homozygous reference only useful for genotypeAll mode
			List<String> alleles = new ArrayList<String>();
			alleles.add(DNASequence.BASES_ARRAY[indexRef]);
			GenomicVariantImpl gvNoVar = new GenomicVariantImpl(pileup.getSequenceName(), pileup.getPosition(), alleles);
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
		SNV snv = new SNV(pileup.getSequenceName(), pileup.getPosition(), refBase, altBase);
		snv.setVariantQS(PhredScoreHelper.calculatePhredScore(refProb));
		CalledSNV csnv = new CalledSNV(snv,genotype);
		csnv.setGenotypeQuality(gq);
		csnv.setTotalReadDepth(countsHelper.getTotalCount());
		csnv.setAllBaseCounts(counts);
		csnv.setAllGenotypeLogConditionals(countsHelper.getLogConditionalProbs());
		return csnv;
	}
	private int [] getIndexesMaxGenotype (double [][] genotypePosteriors, int indexDefault) {
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
	private int[] makeReportCounts(int[] counts, int [] indexes) {
		int [] reportCounts = new int [indexes.length];
		for(int i=0;i<indexes.length;i++) {
			if(indexes[i] <0 || indexes[i] >= counts.length) reportCounts[i] = 0;
			else reportCounts[i] = counts[indexes[i]];
		}
		return reportCounts;
	}
	private double[][] makeReportProbs(double[][] allValues, int [] indexes) {
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
	
	
	private CalledGenomicVariant callIndel (PileupRecord pileup, CountsHelper helper, GenomicVariant variant) {
		int posPrint = -1;
		int [] counts = helper.getCounts();
		double [][] postProbs = helper.getPosteriorProbabilities(heterozygosityRate);
		if(pileup.getPosition()==posPrint) {
			System.out.println("Counts");
			for(int j=0;j<counts.length;j++) System.out.println("Count allele: "+helper.getAlleles()[j]+": "+counts[j]); 
			System.out.println("Probs");
			helper.printProbs(postProbs, false);
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
		//if(pileup.getFirst()==82) System.out.println("Indel alleles: "+newCall.getAlleles().length+" called alleles: "+calledAlleles[0]+" "+calledAlleles[1]+" genotype prob: "+newCall.getGenotypeProbability());
		if(!passFilter(newCall)) newCall.makeUndecided();
		return newCall;
	}
	
	/**
	 * Tells if the given variant passes the filters
	 * @param variant CalledVariant to test
	 * @return boolean true if the variant passes the filters, false otherwise
	 */
	private boolean passFilter (CalledGenomicVariant variant) {
		if(variant==null) return false;
		if(minAltCoverage != DEF_MIN_ALT_COVERAGE || maxAltCoverage != DEF_MAX_ALT_COVERAGE) {
			VariantCallReport report = variant.getCallReport();
			String [] alleles = variant.getAlleles();
			boolean passCoverage = false;
			//Starts in 1 to ignore the reference
			for(int i=1;i<alleles.length;i++) {
				int coverage = report.getCount(alleles[i]);
				boolean pass= true;
				//Min coverage filter
				if(minAltCoverage != DEF_MIN_ALT_COVERAGE && coverage<minAltCoverage) {
					pass = false;
				}
				//Max coverage filter
				if(maxAltCoverage != DEF_MAX_ALT_COVERAGE && coverage>maxAltCoverage) {
					pass = false;
				}
				if(pass) {
					passCoverage = true;
					break;
				}
			}
			if(!passCoverage) return false;
		}
		
		//Probability filter
		if(minQuality != DEF_MIN_QUALITY && variant.getGenotypeQuality()<minQuality) {
			return false;
		}
		return true;
	}
	
	private int nextSIVIndex = 0;
	private List<GenomicVariant> seqInputVariants;
	@Override
	public void onPileup(PileupRecord pileup) {
		if(inputVariants.size()==0) {
			CalledGenomicVariant calledVar = processPileup(pileup,null);
			if(calledVar!=null) {
				calledVariants.add(calledVar);
			}
		} else if(nextSIVIndex<seqInputVariants.size()) {
			GenomicVariant inputVariant = seqInputVariants.get(nextSIVIndex);
			while(inputVariant.getFirst() <= pileup.getPosition() ) {
				if(inputVariant.getFirst()==pileup.getPosition()) {
					CalledGenomicVariant calledVar = processPileup(pileup,inputVariant);
					calledVariants.add(calledVar);
				}
				nextSIVIndex++;
				if(nextSIVIndex>=seqInputVariants.size()) return;
				inputVariant = seqInputVariants.get(nextSIVIndex);
			}
		}
	}
	@Override
	public void onSequenceEnd(String sequenceName) {
		
	}
	@Override
	public void onSequenceStart(String sequenceName) {
		if(inputVariants.size()>0) seqInputVariants = inputVariants.getSequenceRegions(sequenceName).asList();
		nextSIVIndex = 0;
		lastIndelEnd = 0;
	}
	public void clear() {
		calledVariants.clear();
		
	}
}
