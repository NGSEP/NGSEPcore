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
import java.util.List;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;


public class SingleSampleVariantPileupListener implements PileupListener {
	
	public static final double DEF_HETEROZYGOSITY_RATE_DIPLOID = CountsHelper.DEF_HETEROZYGOSITY_RATE_DIPLOID;
	public static final double DEF_HETEROZYGOSITY_RATE_HAPLOID = CountsHelper.DEF_HETEROZYGOSITY_RATE_HAPLOID;
	public static final byte DEF_MAX_BASE_QS = CountsHelper.DEF_MAX_BASE_QS;
	public static final short DEF_MIN_QUALITY = 0;
	public static final String DEF_SAMPLE_ID = "Sample";
	
	private List<CalledGenomicVariant> calledVariants = new ArrayList<CalledGenomicVariant>();
	
	// Parameters 
	private ReferenceGenome genome;
	private double heterozygosityRate = DEF_HETEROZYGOSITY_RATE_DIPLOID;
	private boolean calcStrandBias = false;
	private byte maxBaseQS=DEF_MAX_BASE_QS; 
	private boolean ignoreLowerCaseRef = false;
	private boolean callEmbeddedSNVs = false;
	private short minQuality = DEF_MIN_QUALITY;
	private boolean pool = false;
	
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
	public boolean isPool() {
		return pool;
	}
	public void setPool(boolean pool) {
		this.pool = pool;
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
	 * @param variant Variant to genotype. If null, the method tries to find a new variant given the genome and the pileup
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
		CountsHelper helperSNV = CountsHelper.calculateCountsSNV(calls, maxBaseQS);
		return  VariantDiscoverySNVQAlgorithm.discoverSNV(helperSNV, pileup.getSequenceName(), pileup.getPosition(), reference, heterozygosityRate, calcStrandBias);
	}
	
	private CalledGenomicVariant discoverVariantWithSpan(PileupRecord pileup, String referenceAllele) {
		List<PileupAlleleCall> calls = pileup.getAlleleCalls(referenceAllele.length(),(String)null);
		AlleleCallClustersBuilder acBuilder = new AlleleCallClustersBuilder(pileup.getSequenceName(),pileup.getPosition());
		String [] alleles =  acBuilder.clusterAlleleCalls(pileup, calls, referenceAllele, maxBaseQS);	
		CountsHelper helper = CountsHelper.calculateCounts(alleles,calls, maxBaseQS); 
		CalledGenomicVariant calledVar = VariantDiscoverySNVQAlgorithm.callIndel(pileup, helper, null, heterozygosityRate, calcStrandBias);
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
	
	public CalledGenomicVariant genotypeVariantSample(GenomicVariant variant, PileupRecord pileup,  Sample sample, double h) {
		String referenceAllele = variant.getReference();
		
		CalledGenomicVariant calledVar = null;
		List<PileupAlleleCall> calls = pileup.getAlleleCalls(referenceAllele.length(),sample.getReadGroups());
		if(variant.isSNV()) {
			CountsHelper helperSNV = CountsHelper.calculateCountsSNV(calls, maxBaseQS);
			if(pool) {
				calledVar = genotypeVariantPool(variant, helperSNV);
			} else {
				calledVar = VariantDiscoverySNVQAlgorithm.genotypeSNV(variant, helperSNV, h, false);
			}
		} else {
			CountsHelper helperIndel = CountsHelper.calculateCounts(variant.getAlleles(),calls,maxBaseQS);
			if(pool) {
				calledVar = genotypeVariantPool(variant, helperIndel);
			} else {
				calledVar = VariantDiscoverySNVQAlgorithm.callIndel(pileup, helperIndel, variant, h, false);
			}
			
		}
		if(calledVar==null) calledVar = new CalledGenomicVariantImpl(variant, new byte[0]);
		else if(minQuality>calledVar.getGenotypeQuality()) calledVar.makeUndecided();
		calledVar.setSampleId(sample.getId());
		calledVar.updateAllelesCopyNumberFromCounts(sample.getNormalPloidy());
		return calledVar;
	}
	
	

	private CalledGenomicVariant genotypeVariantPool(GenomicVariant variant, CountsHelper helper) {
		List<Byte> selectedAlleles = new ArrayList<Byte>();
		String [] alleles = variant.getAlleles();
		int total = helper.getTotalCount();
		double threshold = 0.005*total;
		for(int i=0;i<alleles.length;i++) {
			String allele = alleles[i];
			if(helper.getCount(allele)>=threshold) {
				selectedAlleles.add((byte) i);
			}
		}
		byte [] calledAlleles = new byte [selectedAlleles.size()];
		for(int i =0;i<selectedAlleles.size();i++) {
			calledAlleles[i] = selectedAlleles.get(i);
		}
		CalledGenomicVariantImpl call = new CalledGenomicVariantImpl(variant, calledAlleles);
		if(variant.isSNV()) call.setAllCounts(helper.getCounts());
		return call;
	}

}
