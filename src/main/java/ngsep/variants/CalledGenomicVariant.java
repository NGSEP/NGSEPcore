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

public interface CalledGenomicVariant extends GenomicVariant {

	public static final int MAX_NUM_COPIES = 100;
	public static final int MAX_PLOIDY_SAMPLE = Short.MAX_VALUE;
	public static final byte GENOTYPE_UNDECIDED=-1;
	public static final byte GENOTYPE_HOMOREF=0;
	public static final byte GENOTYPE_HETERO=1;
	public static final byte GENOTYPE_HOMOALT=2;
	
	public static final byte ALLELE_UNDECIDED = -1;
	public static final byte ALLELE_REFERENCE = 0;
	public static final byte ALLELE_ALTERNATIVE = 1;
	
	public static final byte MAX_STRAND_BIAS_SCORE = 100;
	public static final byte INVALID_STRAND_BIAS_SCORE = -1;

	/**
	 * @return String sampleId for this call
	 */
	public String getSampleId();

	/**
	 * Changes the sample id
	 * @param sampleId New sample id
	 */
	public void setSampleId(String sampleId);
	
	/**
	 * @return String [] the alleles called for this sample. If the call is undecided, the array is empty
	 */
	public String [] getCalledAlleles();
	/**
	 * @return byte [] the indexes in the array returned by getAlleles() of the alleles called for this sample.
	 * If the call is undecided, the array is empty
	 */
	public byte [] getIndexesCalledAlleles();
	
	/**
	 * @return short [] copy number of each of the alleles in the variant call. 
	 * If the call is undecided, returns an array of zeroes
	 */
	public short[] getAllelesCopyNumber();
	
	/**
	 * Returns the predicted local copy number of the region surrounding the variant.
	 * Unless the call is undecided, it should be the sum of the elements of the array returned by
	 * the method getAllelesCopyNumber
	 * @return short Copy number of the region surrounding this variant 
	 */
	public short getCopyNumber();
	
	/**
	 * @return int total read depth of the variant in this sample
	 */
	public int getTotalReadDepth ();
	
	/**
	 * Changes the total read depth
	 * @param depth New read depth
	 */
	public void setTotalReadDepth (int depth);
	
	/**
	 * @return short quality score of this genotype call
	 */
	public short getGenotypeQuality();
	
	/**
	 * Changes the quality score of this genotype call
	 * @param qualityScore New quality score
	 */
	public void setGenotypeQuality(short genotypeQuality);
	
	/**
	 * Switches the status of the genotype call to an undecided state keeping information from counts and probabilities
	 */
	public void makeUndecided();
	
	/**
	 * @return true if the call is undecided, false otherwise
	 */
	public boolean isUndecided();
	
	/**
	 * @return true if the call is heterozygous, false otherwise
	 */
	public boolean isHeterozygous();
	
	/**
	 * @return true if the call is homozygous, false otherwise
	 */
	public boolean isHomozygous();
	
	/**
	 * @return true if the call is homozygous reference, false otherwise
	 */
	public boolean isHomozygousReference();
	/**
	 * Estimates the copy numbers of the called alleles based on the allele read depth and the given totalCopyNumber
	 * POST: The copy number of each allele is updated with the estimated values. 
	 * The phasing information should be unset because it can become inconsistent with the new total copy number
	 * @param totalCopyNumber predicted copy number of the region surrounding the variant 
	 */
	public void updateAllelesCopyNumberFromCounts(short totalCopyNumber);
	/**
	 * Changes the alleles copy number
	 * PRE: The copy number per allele should be consistent with this genotype
	 * Only indexes of called alleles should contain values different than zero
	 * POST: The copy number of each allele is updated. The phasing information should be
	 * unset because it can become inconsistent with the new total copy number
	 * @param allelesCN Array with the copy number of each allele 
	 */
	public void setAllelesCopyNumber(short [] allelesCN);
	
	/**
	 * @return VariantCallReport with the information gathered while calling the variant in this sample 
	 */
	public VariantCallReport getCallReport();
	
	/**
	 * @return int [] Array with counts which do not need to match the number of alleles of the variant
	 */
	public int[] getAllCounts();
	/**
	 * @return true if the variant is phased
	 */
	public boolean isPhased();
	
	/**
	 * Returns the phasing of the alleles called for this sample
	 * @return String [] Alleles observed in the sample, sorted according with the phasing
	 */
	public String [] getPhasedAlleles();
	/**
	 * Returns the phased alleles as a byte array
	 * @return byte [] indexes in the alleles array of the called alleles sorted according to their phase
	 */
	public byte [] getIndexesPhasedAlleles ();
	/**
	 * Returns the strand bias score
	 * @return byte Phred-scaled probability of strand bias according to the fisher exact test
	 */
	public byte getStrandBiasScore();
}
