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

public class ReadPairCalledGenomicVariant implements CalledGenomicVariant {

	public static final String SOURCE_READPAIR = "Readpairs";
	//Used implementation to allow changing the first and last coordinate
	private GenomicVariantImpl var;
	private String sampleId;
	private byte genotype = CalledSNV.GENOTYPE_UNDECIDED;
	private int predictedLength = 0;
	private int readDepth = 0;
	private int supportingFragments = 0;
	private short genotypeQuality = 0;
	private int numSplitReads = 0;
	
	
	
	
	public ReadPairCalledGenomicVariant(GenomicVariantImpl var, byte genotype, int predictedLength) {
		super();
		this.var = var;
		this.genotype = genotype;
		this.predictedLength = predictedLength;
	}

	
	
	public String getSampleId() {
		return sampleId;
	}



	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}



	public byte getGenotype() {
		return genotype;
	}



	public void setGenotype(byte genotype) {
		this.genotype = genotype;
	}



	public int getPredictedLength() {
		return predictedLength;
	}



	public void setPredictedLength(int predictedLength) {
		this.predictedLength = predictedLength;
	}



	public int getReadDepth() {
		return readDepth;
	}



	public void setReadDepth(int readDepth) {
		this.readDepth = readDepth;
	}



	public int getSupportingFragments() {
		return supportingFragments;
	}



	public void setSupportingFragments(int supportingFragments) {
		this.supportingFragments = supportingFragments;
	}



	public int getNumSplitReads() {
		return numSplitReads;
	}



	public void setNumSplitReads(int numSplitReads) {
		this.numSplitReads = numSplitReads;
	}



	public void setFirst(int first) {
		var.setFirst(first);
	}



	public void setLast(int last) {
		var.setLast(last);
	}
	
	public String getSource() {
		return SOURCE_READPAIR;
	}



	@Override
	public String[] getAlleles() {
		return var.getAlleles();
	}

	@Override
	public String getReference() {
		return var.getReference();
	}

	@Override
	public String getId() {
		return var.getId();
	}

	@Override
	public void setId(String id) {
		var.setId(id);
	}

	@Override
	public short getVariantQS() {
		return var.getVariantQS();
	}

	@Override
	public void setVariantQS(short qualityScore) {
		var.setVariantQS(qualityScore);
	}

	@Override
	public boolean isCompatible(GenomicVariant variant) {
		//TODO: Implement better
		return GenomicVariantImpl.testCompatibility(this, variant, true, false);
	}

	@Override
	public boolean isBiallelic() {
		return true;
	}

	@Override
	public boolean isSNV() {
		return false;
	}

	@Override
	public byte getType() {
		return var.getType();
	}

	@Override
	public void setType(byte type) {
		var.setType(type);
		
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
	public boolean isNegativeStrand() {
		return var.isNegativeStrand();
	}

	@Override
	public String[] getCalledAlleles() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public byte[] getIndexesCalledAlleles() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public short[] getAllelesCopyNumber() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public short getCopyNumber() {
		// TODO Auto-generated method stub
		return 2;
	}

	@Override
	public int getTotalReadDepth() {
		return readDepth;
	}

	@Override
	public void setTotalReadDepth(int depth) {
		readDepth = depth;
	}

	@Override
	public short getGenotypeQuality() {
		return genotypeQuality;
	}

	@Override
	public void setGenotypeQuality(short genotypeQuality) {
		this.genotypeQuality = genotypeQuality;
	}

	@Override
	public void makeUndecided() {
		this.genotype = GENOTYPE_UNDECIDED;
		
	}

	@Override
	public boolean isUndecided() {
		return genotype == GENOTYPE_UNDECIDED;
	}

	@Override
	public boolean isHeterozygous() {
		return genotype == GENOTYPE_HETERO;
	}

	@Override
	public boolean isHomozygous() {
		return genotype == GENOTYPE_HOMOREF || genotype == GENOTYPE_HOMOALT;
	}

	@Override
	public boolean isHomozygousReference() {
		return genotype == GENOTYPE_HOMOREF;
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
		// TODO Auto-generated method stub
		return INVALID_STRAND_BIAS_SCORE;
	}


}
