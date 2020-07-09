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

import ngsep.genome.GenomicRegion;

/**
 * Interface for methods that a variant against a reference genome should implement
 * @author Jorge Duitama
 *
 */
public interface GenomicVariant extends GenomicRegion {
	public static final short DEFAULT_PLOIDY = 2;
	
	public static final byte TYPE_UNDETERMINED = 0;
	public static final byte TYPE_BIALLELIC_SNV = 1;
	public static final byte TYPE_MULTIALLELIC_SNV = 2;
	public static final byte TYPE_EMBEDDED_SNV = 3;
	public static final byte TYPE_INDEL = 4;
	public static final byte TYPE_STR = 5;
	public static final byte TYPE_CNV = 10;
	public static final byte TYPE_REPEAT = 11;
	public static final byte TYPE_LARGEDEL = 12;
	public static final byte TYPE_LARGEINS = 13;
	public static final byte TYPE_INVERSION = 14;
	
	public static final String TYPENAME_BIALLELIC_SNV = "SNV";
	public static final String TYPENAME_MULTIALLELIC_SNV = "MULTISNV";
	public static final String TYPENAME_EMBEDDED_SNV = "EMBEDDED";
	public static final String TYPENAME_INDEL = "INDEL";
	public static final String TYPENAME_STR = "STR";
	public static final String TYPENAME_CNV = "CNV";
	public static final String TYPENAME_REPEAT = "REPEAT";
	public static final String TYPENAME_LARGEDEL = "DEL";
	public static final String TYPENAME_LARGEINS = "INS";
	public static final String TYPENAME_INVERSION = "INV";
	
	
	public static final int MAX_NUM_ALLELES = 100;
	
	/**
	 * Returns the alleles observed for this variant. The reference allele must be the first and all alleles must be
	 * uppercase
	 * @return String [] Alleles observed at this locus. 
	 */
	public String [] getAlleles();
	/**
	 * @return String Uppercase Reference allele for this variant. 
	 */
	public String getReference();
	/**
	 * @return String Variant id
	 */
	public String getId ();
	/**
	 * Changes the id of the variant
	 * @param id New id
	 */
	public void setId (String id);
	/**
	 * @return short Probability of existence of this variant in Phred format
	 */
	public short getVariantQS();
	/**
	 * Changes the value of the quality score
	 * @param qualityScore new quality score
	 */
	public void setVariantQS(short qualityScore);
	/**
	 * Compares this variant with the given for compatibility of position and alleles
	 * @param variant Variant to compare with
	 * @return true if this called variant is located in the same position as the given one and have the same alleles 
	 */
	public boolean isCompatible(GenomicVariant variant);
	/**
	 * @return true if the variant has only two alleles. False otherwise
	 */
	public boolean isBiallelic ();
	/**
	 * @return true if the variant is a SNV, false otherwise
	 */
	public boolean isSNV ();
	/**
	 * @return byte Type of variant
	 */
	public byte getType();
	/**
	 * Changes the variant type
	 * @param type New type. See constants in GenomicVariant
	 */
	public void setType(byte type);
	
}
