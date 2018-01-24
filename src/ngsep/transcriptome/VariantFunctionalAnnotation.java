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
package ngsep.transcriptome;

import ngsep.variants.GenomicVariant;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class VariantFunctionalAnnotation {
	
	
	
	
	
	private GenomicVariant variant;
	private VariantFunctionalAnnotationType type;
	private Transcript transcript;
	private byte altAlleleIdx=-1;
	private int codonNumber=0;
	private byte codonPosition=0;
	private String aminoacidChange;
	
	//PRE: typeName is a type supported by the class VariantFunctionalAnnotationType
	public VariantFunctionalAnnotation(GenomicVariant variant, String typeName) {
		this(variant,VariantFunctionalAnnotationType.getTypeByName(typeName));
	}
	//PRE: type!=null
	public VariantFunctionalAnnotation(GenomicVariant variant, VariantFunctionalAnnotationType type) {
		super();
		if(type==null) throw new NullPointerException("The type of an annotation can not be null");
		this.variant = variant;
		this.type = type;
	}
	/**
	 * @return the variant
	 */
	public GenomicVariant getVariant() {
		return variant;
	}
	
	/**
	 * @return the type
	 */
	public VariantFunctionalAnnotationType getType() {
		return type;
	}
	
	public String getTypeName() {
		return type.getName();
	}
	/**
	 * @return the transcript
	 */
	public Transcript getTranscript() {
		return transcript;
	}
	/**
	 * @param transcript the transcript to set
	 */
	public void setTranscript(Transcript transcript) {
		this.transcript = transcript;
	}

	/**
	 * @param altAlleleIdx the altAlleleIdx to set
	 */
	public void setAltAlleleIdx(byte altAlleleIdx) {
		this.altAlleleIdx = altAlleleIdx;
	}
	/**
	 * @return the altAlleleIdx
	 */
	public int getAltAlleleIdx() {
		return altAlleleIdx;
	}
	/**
	 * @return the codonNumber
	 */
	public int getCodonNumber() {
		return codonNumber;
	}
	/**
	 * @param codonNumber the codonNumber to set
	 */
	public void setCodonNumber(int codonNumber) {
		this.codonNumber = codonNumber;
	}
	/**
	 * @return the codonPosition
	 */
	public byte getCodonPosition() {
		return codonPosition;
	}
	/**
	 * @param codonPosition the codonPosition to set
	 */
	public void setCodonPosition(byte codonPosition) {
		this.codonPosition = codonPosition;
	}
	/**
	 * @return the aminoacidChange
	 */
	public String getAminoacidChange() {
		return aminoacidChange;
	}
	/**
	 * @param aminoacidChange the aminoacidChange to set
	 */
	public void setAminoacidChange(String aminoacidChange) {
		this.aminoacidChange = aminoacidChange;
	}
	/**
	 * @return
	 * @see ngsep.transcriptome.VariantFunctionalAnnotationType#getSoAccession()
	 */
	public String getSoAccession() {
		return type.getSoAccession();
	}
	/**
	 * @return
	 * @see ngsep.transcriptome.VariantFunctionalAnnotationType#isCoding()
	 */
	public boolean isCoding() {
		return type.isCoding();
	}
	
}
