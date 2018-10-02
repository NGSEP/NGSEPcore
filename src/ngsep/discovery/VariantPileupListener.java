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
import java.util.Set;

import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;


public class VariantPileupListener implements PileupListener {
	
	public static final double DEF_HETEROZYGOSITY_RATE_DIPLOID = 0.001;
	public static final double DEF_HETEROZYGOSITY_RATE_HAPLOID = 0.000001;
	public static final byte DEF_MAX_BASE_QS = 100;
	public static final short DEF_MIN_QUALITY = 0;
	
	private List<CalledGenomicVariant> calledVariants = new ArrayList<CalledGenomicVariant>();
	private GenomicRegionSortedCollection<GenomicVariant> inputVariants = new GenomicRegionSortedCollection<GenomicVariant>();
	private ReferenceGenome genome;
	private byte normalPloidy = GenomicVariant.DEFAULT_PLOIDY;
	
	
	private double heterozygosityRate = DEF_HETEROZYGOSITY_RATE_DIPLOID;
	private boolean calcStrandBias = false;
	private byte maxBaseQS=DEF_MAX_BASE_QS; 
	private boolean ignoreLowerCaseRef = false;
	private boolean callEmbeddedSNVs = false;
	private Set<String> readGroups;
	
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
	
	/**
	 * @return the readGroups
	 */
	public Set<String> getReadGroups() {
		return readGroups;
	}
	/**
	 * @param readGroups the readGroups to set
	 */
	public void setReadGroups(Set<String> readGroups) {
		this.readGroups = readGroups;
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
		CountsHelper helperSNV = VariantDiscoverySNVQAlgorithm.calculateCountsSNV(pileup, maxBaseQS, readGroups);
		//if(pileup.getFirst()==82) System.out.println("Pileup last: "+pileup.getLast()+" Reference allele: "+referenceAllele); 
		//if(pileup.getPosition()==9052) System.out.println("Reference allele: "+referenceAllele+". Pileup last: "+pileup.getLast());
		CalledGenomicVariant calledVar;
		if(referenceAllele.length()>1) {
			CountsHelper helperIndel = VariantDiscoverySNVQAlgorithm.calculateCountsIndel(pileup,variant,referenceAllele, readGroups); 
			calledVar = VariantDiscoverySNVQAlgorithm.callIndel(pileup, helperIndel, variant, heterozygosityRate, calcStrandBias);
			if(variant == null) {
				if(calledVar!=null && minQuality!=DEF_MIN_QUALITY && minQuality>calledVar.getGenotypeQuality()) calledVar.makeUndecided();
				if(calledVar!=null && (pileup.isInputSTR() || (!calledVar.isUndecided() && !calledVar.isHomozygousReference()))) {
					//System.out.println("Called indel at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+" variant type: "+calledVar.getType());
					lastIndelEnd = calledVar.getLast();
				} else {
					if (pileup.isNewSTR()) {
						pileup.setSTR(false);
						pileup.setNewSTR(false);
					}
					//Try SNV if the indel alleles were not good to make a call
					calledVar = VariantDiscoverySNVQAlgorithm.callSNV(pileup, helperSNV, variant, referenceAllele.charAt(0), heterozygosityRate, calcStrandBias);
				}
			}
		} else {
			calledVar = VariantDiscoverySNVQAlgorithm.callSNV(pileup, helperSNV, variant, referenceAllele.charAt(0), heterozygosityRate, calcStrandBias);
		}
		if(calledVar!=null && minQuality!=DEF_MIN_QUALITY && minQuality>calledVar.getGenotypeQuality()) calledVar.makeUndecided();
		if(calledVar != null && (variant!=null || (!calledVar.isUndecided() && !calledVar.isHomozygousReference()) )) {
			//System.out.println("Called SNV");
			calledVar.updateAllelesCopyNumberFromCounts(normalPloidy);
			if(calledVar.isSNV() && (pileup.isEmbedded() || calledVar.getFirst()<=lastIndelEnd)) calledVar.setType(GenomicVariant.TYPE_EMBEDDED_SNV);
			//System.out.println("Called variant at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+" pileup variant type: "+pileup.getVariantType()+" variant type: "+calledVar.getType());
			//if(pileup.getPosition()==3892) System.out.println("New variant at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+". Repeat context: "+calledVar.getRepeatContext());
			return calledVar;
		}
		return null;
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
}
