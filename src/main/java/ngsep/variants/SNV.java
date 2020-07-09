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

import ngsep.sequences.DNASequence;

/**
 * Information for a single nucleotide variant (SNV).
 * Implements the GenomicRegion interface by assigning position to both the
 * first and the last coordiantes
 * This implementation is meant to be space efficient at the expense of usefulness. 
 * The main limitation is that it can only handle one reference and one alternative allele
 * @author Jorge Duitama
 */
public class SNV implements GenomicVariant { 
	private String id;
	private String sequenceName;
	private int position;
	private byte refBaseIdx; //0 for A, 1 for C, 2 for G and 3 for T 
	private byte altBaseIdx; //0 for A, 1 for C, 2 for G and 3 for T
	private boolean negativeStrand = false;
	private short variantQS=0;
	private byte type = GenomicVariant.TYPE_BIALLELIC_SNV;
	
	/**
	 * Creates a new SNV with the given information
	 * @param sequenceName Name of the sequence where the SNV is located
	 * @param position Position in the sequence where the SNV is located
	 * @param refBase Base present in the reference sequence in the given position
	 * @param altBase Alternative base for the given position 
	 */
	public SNV (String sequenceName, int position, char refBase, char altBase) {
		this.setSequenceName(sequenceName);
		this.setPosition(position);
		this.setReferenceBase(refBase);
		this.setAlternativeBase(altBase);
	}
	/**
	 * @return String Name of the sequence where the SNV is located
	 */
	public String getSequenceName() {
		return sequenceName;
	}
	/**
	 * Changes the name of the sequence where the SNV is located
	 * @param sequenceName New name of the sequence
	 */
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
	/**
	 * @return int Position in the sequence where the SNV is located
	 */
	public int getPosition() {
		return position;
	}
	/**
	 * Changes the position in the sequence where the SNV is located
	 * @param position New coordinate
	 */
	public void setPosition(int position) {
		this.position = position;
	}
	/**
	 * @return char Reference base in upper case letters
	 */
	public char getReferenceBase() {
		return DNASequence.BASES_STRING.charAt(refBaseIdx);
	}
	/**
	 * @return int 0 if the reference base is A, 1 if it is C, 2 if it is G, and 3 if it is T
	 */
	public byte getRefBaseDNAIndex() {
		return refBaseIdx;
	}
	/**
	 * Changes the reference base
	 * @param refBase New reference base
	 */
	public void setReferenceBase(char refBase) {
		char upRef = Character.toUpperCase(refBase);
		byte idx = (byte)DNASequence.BASES_STRING.indexOf(upRef);
		if(idx<0) throw new IllegalArgumentException("Allele "+refBase+ " is not supported as reference allele in a SNV");
		this.refBaseIdx = idx;
	}
	
	/**
	 * @return char Alternative base in uppercase.
	 */
	public char getAlternativeBase() {
		return DNASequence.BASES_STRING.charAt(altBaseIdx);
	}
	/**
	 * @return byte 0 if the alternative base is A, 1 if it is C, 2 if it is G, and 3 if it is T
	 */
	public byte getAltBaseDNAIndex() {
		return altBaseIdx;
	}
	/**
	 * Changes the alternative base
	 * @param altBase New alternative base
	 */
	public void setAlternativeBase(char altBase) {
		char upAlt = Character.toUpperCase(altBase);
		byte idx = (byte)DNASequence.BASES_STRING.indexOf(upAlt);
		if(idx<0) throw new IllegalArgumentException("Allele "+altBase+ " is not supported as alternative allele in a SNV");
		this.altBaseIdx = idx;
	}
	
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public short getVariantQS() {
		return variantQS;
	}
	/**
	 * Changes the quality score of the variant
	 * @param variantQS New variant quality score
	 */
	public void setVariantQS(short variantQS) {
		this.variantQS = variantQS;
	}
	public boolean isTransition() {
		return Math.abs(getRefBaseDNAIndex()-getAltBaseDNAIndex())==2;
	}
	public boolean isTransversion() {
		return !isTransition();
	}
	
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
	@Override
	public int getFirst() {
		return position;
	}
	@Override
	public int getLast() {
		return position;
	}
	@Override
	public int length() {
		return 1;
	}
	@Override
	public boolean isPositiveStrand() {
		return !negativeStrand;
	}
	@Override
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	@Override
	public String [] getAlleles() {
		String [] answer = new String [2];
		answer[0] = DNASequence.BASES_ARRAY[refBaseIdx];
		answer[1] = DNASequence.BASES_ARRAY[altBaseIdx];
		return answer;
	}
	@Override
	public String getReference() {
		return DNASequence.BASES_ARRAY[refBaseIdx];
	}
	@Override
	public boolean isCompatible(GenomicVariant variant) {
		return GenomicVariantImpl.testCompatibility(this, variant,true, true);
	}
	@Override
	public boolean isBiallelic() {
		return true;
	}
	@Override
	public boolean isSNV() {
		return true;
	}
	@Override
	public byte getType() {
		return type;
	}
	@Override
	public void setType(byte type) {
		this.type = type;
	}
}
