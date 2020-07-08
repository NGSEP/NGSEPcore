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
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.math.LogMath;
import ngsep.math.NumberArrays;
import ngsep.math.PhredScoreHelper;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.variants.VariantCallReport;


public class SingleSampleVariantPileupListener implements PileupListener {
	
	public static final double DEF_HETEROZYGOSITY_RATE_DIPLOID = CountsHelper.DEF_HETEROZYGOSITY_RATE_DIPLOID;
	public static final double DEF_HETEROZYGOSITY_RATE_HAPLOID = CountsHelper.DEF_HETEROZYGOSITY_RATE_HAPLOID;
	public static final byte DEF_MAX_BASE_QS = CountsHelper.DEF_MAX_BASE_QS;
	public static final short DEF_MIN_QUALITY = 0;
	public static final String DEF_SAMPLE_ID = "Sample";
	public static final byte DEF_MIN_PLOIDY_POOL_ALGORITHM = 20;
	
	private List<CalledGenomicVariant> calledVariants = new ArrayList<CalledGenomicVariant>();
	
	// Parameters 
	private ReferenceGenome genome;
	private double heterozygosityRate = DEF_HETEROZYGOSITY_RATE_DIPLOID;
	private boolean calcStrandBias = false;
	private byte maxBaseQS=DEF_MAX_BASE_QS; 
	private boolean ignoreLowerCaseRef = false;
	private boolean callEmbeddedSNVs = false;
	private short minQuality = DEF_MIN_QUALITY;
	
	// Parameters only for listener mode
	private GenomicRegionSortedCollection<GenomicVariant> inputVariants = new GenomicRegionSortedCollection<GenomicVariant>();
	private Sample sample = new Sample(DEF_SAMPLE_ID);
	
	private static int posPrint = -1;
	
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	
	public Sample getSample() {
		return sample;
	}
	public void setSample(Sample sample) {
		this.sample = sample;
	}
	public double getHeterozygosityRate() {
		return heterozygosityRate;
	}
	public void setHeterozygosityRate(double heterozygosityRate) {
		this.heterozygosityRate = heterozygosityRate;
	}
	/**
	 * @return the calcStrandBias
	 */
	public boolean isCalcStrandBias() {
		return calcStrandBias;
	}
	/**
	 * @param calcStrandBias the calcStrandBias to set
	 */
	public void setCalcStrandBias(boolean calcStrandBias) {
		this.calcStrandBias = calcStrandBias;
	}
	public byte getMaxBaseQS() {
		return maxBaseQS;
	}
	public void setMaxBaseQS(byte maxBaseQS) {
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
	
	/**
	 * @return the minQuality
	 */
	public short getMinQuality() {
		return minQuality;
	}
	/**
	 * @param minQuality the minQuality to set
	 */
	public void setMinQuality(short minQuality) {
		this.minQuality = minQuality;
	}
	public List<GenomicVariant> getInputVariants() {
		return inputVariants.asList();
	}
	public void setInputVariants(GenomicRegionSortedCollection<GenomicVariant> inputVariants) {
		this.inputVariants = inputVariants;
	}
	public List<CalledGenomicVariant> getCalledVariants() {
		return calledVariants;
	}
	//Control attribute to avoid calling overlapping indels and to give an embedded status to SNVs within indels or STRs
	private int lastIndelEnd = 0;
	private int nextSIVIndex = 0;
	private List<GenomicVariant> seqInputVariants;
	@Override
	public void onPileup(PileupRecord pileup) {
		if(inputVariants.size()==0) {
			if(pileup.isInputSTR()) lastIndelEnd = pileup.getPosition()+pileup.getReferenceSpan()-1;
			else if(pileup.getPosition()<=lastIndelEnd) pileup.setEmbedded(true);
			String referenceAllele = SingleSampleVariantPileupListener.calculateReferenceAlleleDiscovery(pileup,genome,callEmbeddedSNVs,ignoreLowerCaseRef);
			if(referenceAllele == null) return;
			CalledGenomicVariant calledVar = discoverVariant(pileup, referenceAllele);
			if(calledVar!=null) {
				calledVariants.add(calledVar);
				if(!calledVar.isSNV() && !calledVar.isUndecided() && !calledVar.isHomozygousReference()) {
					lastIndelEnd = calledVar.getLast();
				}
			}
		} else if(nextSIVIndex<seqInputVariants.size()) {
			GenomicVariant inputVariant = seqInputVariants.get(nextSIVIndex);
			while(inputVariant.getFirst() <= pileup.getPosition() ) {
				if(inputVariant.getFirst()==pileup.getPosition()) {
					CalledGenomicVariant calledVar = genotypeVariantSample(inputVariant, pileup, sample, heterozygosityRate);
					calledVariants.add(calledVar);
				}
				nextSIVIndex++;
				if(nextSIVIndex>=seqInputVariants.size()) return;
				inputVariant = seqInputVariants.get(nextSIVIndex);
			}
		}
	}
	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {
		
	}
	@Override
	public void onSequenceStart(QualifiedSequence sequence) {
		if(inputVariants.size()>0) seqInputVariants = inputVariants.getSequenceRegions(sequence.getName()).asList();
		nextSIVIndex = 0;
		lastIndelEnd = 0;
	}
	public void clear() {
		calledVariants.clear();
	}
	static String calculateReferenceAlleleDiscovery(PileupRecord pileup, ReferenceGenome genome, boolean callEmbeddedSNVs, boolean ignoreLowerCaseRef) {
		String referenceAllele;
		if(!callEmbeddedSNVs && pileup.isEmbedded()) return null;
		int last = pileup.getPosition()+pileup.getReferenceSpan()-1;
		CharSequence seq = genome.getReference(pileup.getSequenceName(), pileup.getPosition(), last);
		if(seq == null) return null;
		referenceAllele = seq.toString();
		if(ignoreLowerCaseRef && Character.isLowerCase(referenceAllele.charAt(0))) return null;
		referenceAllele = referenceAllele.toUpperCase();
		
		if(pileup.isEmbedded()) {
			referenceAllele = referenceAllele.substring(0,1);
			pileup.setSTR(false);
		}
		return referenceAllele;
	}
	/**
	 * 
	 * @param pileup
	 * @param variant Variant to discover and genotype. If null, the method tries to find a new variant given the genome and the pileup
	 * @return CalledGenomicVariant. New call of the given variant or new called variant if exists
	 */
	public CalledGenomicVariant discoverVariant(PileupRecord pileup, String referenceAllele) {
		
		if(pileup.getPosition()==posPrint) System.out.println("Processing pileup at "+pileup.getSequenceName()+":"+pileup.getPosition()+" span: "+pileup.getReferenceSpan()+" reference: "+referenceAllele);
		CalledGenomicVariant calledVar;
		if(referenceAllele.length()>1) {
			calledVar = discoverVariantWithSpan(pileup, referenceAllele);
		} else {
			calledVar = discoverSNV(pileup, referenceAllele.charAt(0));
		}
		//Ignore call if it is not variant with enough quality. To change if the genotype all function is brought back
		if(calledVar!=null && (calledVar.isUndecided() || calledVar.isHomozygousReference() || minQuality>calledVar.getGenotypeQuality())) calledVar = null;
		if(calledVar != null) {
			calledVar.setSampleId(sample.getId());
			calledVar.updateAllelesCopyNumberFromCounts(sample.getNormalPloidy());
			if(calledVar.isSNV() && pileup.isEmbedded()) calledVar.setType(GenomicVariant.TYPE_EMBEDDED_SNV);	
		}
		//System.out.println("Called variant at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+" pileup variant type: "+pileup.getVariantType()+" variant type: "+calledVar.getType());
		//if(pileup.getPosition()==3892) System.out.println("New variant at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+". Repeat context: "+calledVar.getRepeatContext());
		return calledVar;
	}
	
	public CalledGenomicVariant discoverSNV(PileupRecord pileup, char reference) {
		List<PileupAlleleCall> calls = pileup.getAlleleCalls(1,(String)null);
		CountsHelper helperSNV = CountsHelper.calculateCountsSNV(calls, maxBaseQS, 0.5);
		short ploidy = sample.getNormalPloidy();
		if(ploidy<DEF_MIN_PLOIDY_POOL_ALGORITHM) {
			return  VariantDiscoverySNVQAlgorithm.discoverSNV(helperSNV, pileup.getSequenceName(), pileup.getPosition(), reference, heterozygosityRate, calcStrandBias);
		} else {
			GenomicVariant variant = createSNVVariantPool(pileup, helperSNV, reference, 0.5/(double)sample.getNormalPloidy());
			if(variant == null) return null;
			CalledGenomicVariant call = genotypeVariantPool(variant, ploidy, calls, heterozygosityRate);
			if(variant.isBiallelic()) return call;
			if(call.isUndecided() || call.isHomozygousReference()) return null;
			byte [] indexes = call.getIndexesCalledAlleles();
			if(indexes.length>2) return call;
			if(indexes.length==2 && indexes[0]!=0) return call;
			Set<String> newAlleles = new LinkedHashSet<String>();
			newAlleles.add(variant.getReference());
			newAlleles.addAll(Arrays.asList(call.getCalledAlleles()));
			variant = makeNewVariant(variant, newAlleles);
			return genotypeVariantPool(variant, ploidy, calls, heterozygosityRate);
		}
	}
	
	private CalledGenomicVariant discoverVariantWithSpan(PileupRecord pileup, String referenceAllele) {
		List<PileupAlleleCall> calls = pileup.getAlleleCalls(referenceAllele.length(),(String)null);
		AlleleCallClustersBuilder acBuilder = new AlleleCallClustersBuilder(pileup.getSequenceName(),pileup.getPosition());
		String [] alleles =  acBuilder.clusterAlleleCalls(pileup, calls, referenceAllele, maxBaseQS);
		CalledGenomicVariant calledVar = discoverIndel(pileup, alleles, calls); 
		//Ignore call if it is not variant with enough quality
		if(calledVar!=null && (calledVar.isUndecided() || calledVar.isHomozygousReference() || minQuality>calledVar.getGenotypeQuality())) calledVar = null;
		if(!pileup.isInputSTR() && calledVar==null) {
			if (pileup.isNewSTR()) {
				pileup.setSTR(false);
				pileup.setNewSTR(false);
			}
			//Try SNV if the indel alleles were not good to make a call
			calledVar = discoverSNV(pileup, referenceAllele.charAt(0));
		}
		return calledVar;
	}
	
	private CalledGenomicVariant discoverIndel(PileupRecord pileup, String [] alleles, List<PileupAlleleCall> calls) {
		CountsHelper helper = CountsHelper.calculateCounts(alleles, calls, maxBaseQS, 0.5);
		short ploidy = sample.getNormalPloidy();
		if(ploidy<DEF_MIN_PLOIDY_POOL_ALGORITHM) {
			return VariantDiscoverySNVQAlgorithm.callIndel(pileup, helper, null, heterozygosityRate, calcStrandBias);
		} else {
			GenomicVariant indel = createIndelVariantPool(pileup, helper);
			if(indel == null || allelesSameLength(indel.getAlleles())) return null;
			CalledGenomicVariant call = genotypeVariantPool(indel, ploidy, calls, heterozygosityRate);
			if(indel.isBiallelic()) return call;
			if(call.isUndecided() || call.isHomozygousReference()) return null;
			byte [] indexes = call.getIndexesCalledAlleles();
			if(indexes.length>2) return call;
			if(indexes.length==2 && indexes[0]!=0) return call;
			Set<String> newAlleles = new LinkedHashSet<String>();
			newAlleles.add(indel.getReference());
			newAlleles.addAll(Arrays.asList(call.getCalledAlleles()));
			indel = makeNewVariant(indel, newAlleles);
			if(indel == null || allelesSameLength(indel.getAlleles())) return null;
			return genotypeVariantPool(indel, ploidy, calls, heterozygosityRate);
		}
	}
	public static GenomicVariant createSNVVariantPool(PileupRecord pileup, CountsHelper helper, char reference, double minAlleleDepthFrequency) {
		if(helper.getTotalCount()==0) {
			return null;
		}
		int refIdx = DNASequence.BASES_STRING.indexOf(reference);
		if(refIdx<0) {
			//N reference can in principle be handled but it generates  many non variant sites
			return null;
		}
		//Simple method based on relative counts to gather alleles.
		int [] counts = helper.getCounts();
		int sum = NumberArrays.getSum(counts); 
		double minCount = Math.max(1, minAlleleDepthFrequency*sum);
		if(pileup.getPosition()==posPrint) System.out.println("Refidx: "+refIdx+" sum: "+sum);
		boolean [] allelesSupported = new boolean [ counts.length];
		List<String> alleles = new ArrayList<>();
		alleles.add(DNASequence.BASES_ARRAY[refIdx]);
		for(int i=0;i<counts.length;i++) {
			allelesSupported[i]= counts[i] >= minCount;
			if(allelesSupported[i] && i!=refIdx) {
				alleles.add(DNASequence.BASES_ARRAY[i]);
			}
		}
		if(pileup.getPosition()==posPrint) System.out.println("Alleles: "+alleles);
		GenomicVariant variant = null;
		if(alleles.size()==2) {
			variant = new SNV(pileup.getSequenceName(), pileup.getPosition(), reference, alleles.get(1).charAt(0));
			variant.setType(GenomicVariant.TYPE_BIALLELIC_SNV);
		} else if (alleles.size()>2){
			//Perform genotyping until two alleles remain or all alleles have enough support in at least one sample
			variant = new GenomicVariantImpl(pileup.getSequenceName(), pileup.getPosition(), alleles);
			variant.setType(GenomicVariant.TYPE_MULTIALLELIC_SNV);
		}
		return variant;
	}
	public static GenomicVariant createIndelVariantPool(PileupRecord pileup, CountsHelper helperIndel) {
		List<String> alleles = helperIndel.getAllelesList();
		if(alleles.size() == 1) return null;
		if(helperIndel.getTotalCount()==0) return null;
		return new GenomicVariantImpl(pileup.getSequenceName(),pileup.getPosition(),alleles);
	}
	public static boolean allelesSameLength(String[] alleles) {
		int l = alleles[0].length();
		for(String allele:alleles) {
			if(allele.length()!=l) return false;
		}
		return true;
	}
	public static GenomicVariant makeNewVariant(GenomicVariant variant, Set<String> newAlleles) {
		if(variant.getFirst()==posPrint) System.out.println("Recoding alleles for "+variant.getFirst()+" alleles: "+Arrays.asList(variant.getAlleles()));
		List<String> alleles = new ArrayList<>(newAlleles.size());
		String reference = variant.getReference(); 
		alleles.add(reference);
		for(String allele:newAlleles) {
			if(!allele.equals(reference)) alleles.add(allele);
		}
		if(variant.getFirst()==posPrint) System.out.println("New alleles: "+alleles);
		if(variant.isSNV() && alleles.size()==2) {
			return new SNV(variant.getSequenceName(), variant.getFirst(), reference.charAt(0), alleles.get(1).charAt(0));
		}
		return new GenomicVariantImpl(variant.getSequenceName(), variant.getFirst(), alleles);
	}
	
	public CalledGenomicVariant genotypeVariantSample(GenomicVariant variant, PileupRecord pileup,  Sample sample, double h) {
		String referenceAllele = variant.getReference();
		short ploidy = sample.getNormalPloidy();
		CalledGenomicVariant calledVar = null;
		List<PileupAlleleCall> calls = pileup.getAlleleCalls(referenceAllele.length(),sample.getReadGroups());
		if(variant.isSNV()) {
			if(ploidy>=DEF_MIN_PLOIDY_POOL_ALGORITHM) {
				calledVar = genotypeVariantPool(variant, ploidy, calls, h);
			} else {
				CountsHelper helperSNV = CountsHelper.calculateCountsSNV(calls, maxBaseQS, 0.5);
				calledVar = VariantDiscoverySNVQAlgorithm.genotypeSNV(variant, helperSNV, h, false);
				calledVar.updateAllelesCopyNumberFromCounts(ploidy);
			}
		} else {
			if(ploidy>=DEF_MIN_PLOIDY_POOL_ALGORITHM) {
				calledVar = genotypeVariantPool(variant, ploidy, calls, h);
			} else {
				CountsHelper helperIndel = CountsHelper.calculateCounts(variant.getAlleles(), calls, maxBaseQS, 0.5);
				calledVar = VariantDiscoverySNVQAlgorithm.callIndel(pileup, helperIndel, variant, h, false);
				calledVar.updateAllelesCopyNumberFromCounts(ploidy);
			}
		}
		if(calledVar==null) calledVar = new CalledGenomicVariantImpl(variant, new byte[0]);
		else if(minQuality>calledVar.getGenotypeQuality()) calledVar.makeUndecided();
		calledVar.setSampleId(sample.getId());		
		return calledVar;
	}
	
	/**
	 * Simple algorithm based on pool thresholds to genotype a variant in a pool
	 * @param variant to genotype
	 * @param sample to genotype as a pool
	 * @param ploidy
	 * @param calls
	 * @param h prior heterozygosity rate
	 * @return CalledGenomicVariant Genotype call for the given pool at the given variant
	 */
	private CalledGenomicVariant genotypeVariantPool(GenomicVariant variant, short haplotypes, List<PileupAlleleCall> calls, double h) {
		double step = 1.0/(double)haplotypes;
		List<Byte> selectedAlleles = new ArrayList<Byte>();
		String [] alleles = variant.getAlleles();
		//Calculate counts for different hypotheses of heterozygosity
		//TODO: Make it more smart based on total read depth
		List<Double> freqs = new ArrayList<Double>();
		List<CountsHelper> helpers = new ArrayList<CountsHelper>();
		for(double freq = step;freq<0.51;freq+=step) {
			freqs.add(freq);
			helpers.add(CountsHelper.calculateCounts(alleles, calls, maxBaseQS, freq));
		}
		//Select the first to obtain counts and most frequent allele
		CountsHelper helper = helpers.get(0);
		int [] counts = helper.getCounts();
		int idxMaxFreq = NumberArrays.getIndexMaximum(counts);
		if(counts[idxMaxFreq]<haplotypes) return new CalledGenomicVariantImpl(variant, new byte[0]);
		//Save most frequent allele
		selectedAlleles.add((byte)idxMaxFreq);
		//Calculate priors
		int heteroGenotypes = helpers.size();
		double logPriorHetero = Math.log10(h/heteroGenotypes);
		double logPriorHomo = Math.log10((1-h));
		int numHypotheses = helpers.size()+1;
		double [] terms = new double[numHypotheses];
		double termHomozygous = LogMath.logProduct(helper.getLogConditionalProbs()[idxMaxFreq][idxMaxFreq],logPriorHomo);
		double maxHetPosterior = 0;
		double minHomoPosterior = 1;
		int maxFreqIdx = 0;
		double maxFreq = 0;
		byte maxAltAllele = -1;
		for(byte i=0;i<alleles.length;i++) {
			if (i==idxMaxFreq) continue;
			//Recover conditional different hypothesis
			terms[0] = termHomozygous;
			for(int j=0;j<freqs.size();j++) {
				CountsHelper helperF = helpers.get(j);
				terms[j+1] = LogMath.logProduct(helperF.getLogConditionalProbs()[i][idxMaxFreq],logPriorHetero);
			}
			CountsHelper.calculatePosteriorProbabilities(terms);
			int idxMax = NumberArrays.getIndexMaximum(terms);
			if(idxMax==0) {
				minHomoPosterior = Math.min(minHomoPosterior, terms[0]);
				continue;
			}
			if(maxAltAllele==-1 || maxHetPosterior<terms[idxMax]) {
				maxHetPosterior = terms[idxMax];
				maxFreqIdx = idxMax-1;
				maxFreq = freqs.get(idxMax-1);
				maxAltAllele = i;
			}
		}
		if(maxAltAllele!=-1) selectedAlleles.add((byte)maxAltAllele);
		Collections.sort(selectedAlleles);
		byte [] calledAlleles = new byte [selectedAlleles.size()];
		
		for(int i =0;i<selectedAlleles.size();i++) {
			calledAlleles[i] = selectedAlleles.get(i);
		}
		CalledGenomicVariantImpl call = new CalledGenomicVariantImpl(variant, calledAlleles);
		call.setTotalReadDepth(helper.getTotalCount());
		short [] acn = new short [alleles.length];
		Arrays.fill(acn, (short)0);
		if(maxAltAllele==-1) {
			//Homozygous genotype call
			call.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(1-minHomoPosterior));
			acn[idxMaxFreq] = haplotypes;
		} else {
			//Heterozygous genotype call
			call.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(1-maxHetPosterior));
			helper = helpers.get(maxFreqIdx);
			
			short altCN = (short) Math.round(maxFreq*haplotypes);
			if(altCN==0) altCN++;
			else if (altCN == haplotypes) altCN--;
			acn[maxAltAllele] = altCN;
			acn[idxMaxFreq] = (short) (haplotypes-altCN);
		}
		call.setAllelesCopyNumber(acn);
		VariantCallReport report = new VariantCallReport(alleles, counts, helper.getLogConditionalProbs());
		call.setCallReport(report);
		
		return call;
	}

}
