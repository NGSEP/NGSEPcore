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
package ngsep.variants;

import java.util.Arrays;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegionSpanComparator;

public class CalledCNV implements CalledGenomicVariant {
	
	private static final byte GENOTYPE_UNDEFINED = -1;
	public static final String TEXT_GEN_DEL = "DEL";
	public static final String TEXT_GEN_TANDEMDUP = "TANDEMDUP";
	public static final String TEXT_GEN_TRANSDUP = "TRANSDUP";
	
	private GenomicVariant var;
	private String sampleId;
	private byte genotype = GENOTYPE_UNDEFINED;
	private short genotypeQuality=0;
	private int uniqueAlns=0;
	private int nonUniqueAlns=0;
	private int depth=0;
	private float numCopies=0;
	
	private int heterozygousVariants=0;
	private int tandemFragments = 0;
	private int transDupFragments = 0;
	private String textGenotype = null;
	private String source = null;
	
	/**
	 * Creates a called CNV with an undefined genotype
	 * @param var CNV called by this object
	 */
	public CalledCNV (GenomicVariant var) {
		this.var = var;
		this.genotype = GENOTYPE_UNDEFINED;
	}
	public CalledCNV (GenomicVariant var, byte genotype) {
		this.var = var;
		if(genotype!=GENOTYPE_UNDEFINED) {
			this.setNumCopies(genotype,false);
			//Changed if the genotype is larger than the maximum number of copies
			genotype = (byte)this.getNumCopies();
		}
		this.genotype = genotype;
	}
	
	public CalledCNV (GenomicVariant var, float numCopies) {
		this.var = var;
		this.setNumCopies(numCopies,true);
	}
	

	public String getId() {
		return var.getId();
	}

	@Override
	public String getSequenceName() {
		return var.getSequenceName();
	}

	@Override
	public int getFirst() {
		return var.getFirst();
	}

	@Override
	public int getLast() {
		return var.getLast();
	}

	@Override
	public int length() {
		return var.length();
	}

	@Override
	public boolean isPositiveStrand() {
		return var.isPositiveStrand();
	}
	@Override
	public String[] getAlleles() {
		return var.getAlleles();
	}

	@Override
	public String getReference() {
		return var.getReference();
	}

	public void setId(String id) {
		var.setId(id);
	}

	public short getVariantQS() {
		return var.getVariantQS();
	}

	public void setVariantQS(short qualityScore) {
		var.setVariantQS(qualityScore);
	}

	public boolean isCompatible(GenomicVariant variant) {
		return var.isCompatible(variant);
	}

	public boolean isBiallelic() {
		return var.isBiallelic();
	}

	public boolean isSNV() {
		return var.isSNV();
	}

	@Override
	public boolean isNegativeStrand() {
		return var.isNegativeStrand();
	}

	public String getSampleId() {
		return sampleId;
	}

	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}
	

	@Override
	public short getGenotypeQuality() {
		return genotypeQuality;
	}
	@Override
	public void setGenotypeQuality(short genotypeQuality) {
		this.genotypeQuality = genotypeQuality;
	}

	public int getUniqueAlns() {
		return uniqueAlns;
	}

	public void setUniqueAlns(int uniqueAlns) {
		this.uniqueAlns = uniqueAlns;
	}

	public int getNonUniqueAlns() {
		return nonUniqueAlns;
	}

	public void setNonUniqueAlns(int nonUniqueAlns) {
		this.nonUniqueAlns = nonUniqueAlns;
	}
	
	public float getNumCopies() {
		return numCopies;
	}

	public void setNumCopies(float numCopies, boolean updateGenotype) {
		if(numCopies<0) this.numCopies = 0;
		if(numCopies>CalledGenomicVariant.MAX_NUM_COPIES) this.numCopies = CalledGenomicVariant.MAX_NUM_COPIES; 
		else this.numCopies = numCopies;
		if(updateGenotype) this.genotype = (byte)Math.round(this.getNumCopies());
	}

	public int getHeterozygousVariants() {
		return heterozygousVariants;
	}

	public void setHeterozygousVariants(int heterozygousVariants) {
		this.heterozygousVariants = heterozygousVariants;
	}

	public void addHeterozygousVariant() {
		heterozygousVariants++;
	}

	

	public byte getType() {
		return var.getType();
	}
	
	public void setType(byte type) {
		var.setType(type);
	}

	public String getSource() {
		return source;
	}

	public void setSource(String source) {
		this.source = source;
	}

	public int getTandemFragments() {
		return tandemFragments;
	}

	public void setTandemFragments(int tandemFragments) {
		this.tandemFragments = tandemFragments;
	}

	public int getTransDupFragments() {
		return transDupFragments;
	}

	public void setTransDupFragments(int transDupFragments) {
		this.transDupFragments = transDupFragments;
	}
	
	public String getTextGenotype() {
		return textGenotype;
	}

	public void setTextGenotype(String textGenotype) {
		this.textGenotype = textGenotype;
	}

	/**
	 * Updates counts of fragments within this CNV given the alignment
	 * @param aln Alignment to process. PRE: The alignment must be located in the same sequence as this CNV
	 * @param avgInsertLength Average insert length of the library corresponding to the input alignment
	 * @return boolean true if the alignment was accepted
	 */
	public boolean addPairedFragment (ReadAlignment aln, int avgInsertLength) {
		if(aln.isReadUnmapped()) return false;
		if(!GenomicRegionSpanComparator.getInstance().span(this, aln)) return false;
		if(aln.isMateUnmapped()) return false;
		int maxInsLenTandem = 3*Math.max(this.length(), avgInsertLength);
		int diffStartAlns = Math.abs(aln.getFirst()-aln.getMateFirst());
		if(aln.isMateSameSequence() && GenomicRegionSpanComparator.getInstance().span(this, aln.getMateFirst(), aln.getMateFirst()+aln.getReadLength())) {
			if(aln.isPositiveStrand() && aln.getFirst()>aln.getMateFirst()) {
				tandemFragments++;
				return true;
			}
			if(aln.isNegativeStrand() && aln.getFirst()<aln.getMateFirst()) {
				tandemFragments++;
				return true;
			}
		} else if(aln.isMateDifferentSequence() || diffStartAlns>maxInsLenTandem) {
			transDupFragments++;
			return true;
		}
		return false;
	}

	@Override
	public byte[] getIndexesCalledAlleles() {
		if(genotype == GENOTYPE_UNDEFINED) return new byte[0];
		byte [] answer = new byte[1];
		answer[0] = genotype;
		return answer;
	}
	@Override
	public String[] getCalledAlleles() {
		if(genotype == GENOTYPE_UNDEFINED) return new String[0];
		String [] alleles = new String[1];
		alleles[0] = ""+genotype;
		return alleles;
	}
	@Override
	public short[] getAllelesCopyNumber() {
		short [] answer = new short[getAlleles().length];
		Arrays.fill(answer, (byte)0);
		if(genotype !=GENOTYPE_UNDEFINED) answer[genotype] = 1;
		return answer;
	}

	@Override
	public short getCopyNumber() {
		return (short) Math.round(numCopies);
	}

	@Override
	public int getTotalReadDepth() {
		return depth;
	}

	@Override
	public void setTotalReadDepth(int depth) {
		this.depth = depth;
	}

	@Override
	public void makeUndecided() {
		genotype = GENOTYPE_UNDEFINED;
		genotypeQuality = 0;
	}

	@Override
	public boolean isUndecided() {
		return genotype == GENOTYPE_UNDEFINED;
	}

	@Override
	public boolean isHeterozygous() {
		return false;
	}

	@Override
	public boolean isHomozygous() {
		return !isUndecided();
	}

	@Override
	public boolean isHomozygousReference() {
		return genotype==DEFAULT_PLOIDY;
	}

	@Override
	public void updateAllelesCopyNumberFromCounts(short totalCopyNumber) {
		// TODO Auto-generated method stub
	}

	@Override
	public void setAllelesCopyNumber(short[] allelesCN) {
		// TODO Auto-generated method stub
	}

	@Override
	public VariantCallReport getCallReport() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int[] getAllCounts() {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public boolean isPhased() {
		// TODO Auto-generated method stub
		return false;
	}
	@Override
	public String[] getPhasedAlleles() {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public byte[] getIndexesPhasedAlleles() {
		// TODO Auto-generated method stub
		return null;
	}
	@Override
	public byte getStrandBiasScore() {
		return 0;
	}
}
