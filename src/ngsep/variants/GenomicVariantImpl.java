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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.sequences.DNASequence;

public class GenomicVariantImpl implements GenomicVariant {
	private String id;
	private String sequenceName;
	private int first;
	private int last;
	private int length;
	private List<String> alleles=new ArrayList<String>();
	private boolean negativeStrand = false;
	private short variantQS=0;
	private byte type = GenomicVariant.TYPE_UNDETERMINED;
	
	private static Map<String, Byte> variantTypesByName = null;
	private static Map<Byte, String> variantTypesById = null;
	private static void buildVariantTypeMaps () {
		variantTypesByName = new HashMap<String, Byte>(10);	
		variantTypesByName.put(TYPENAME_BIALLELIC_SNV, TYPE_BIALLELIC_SNV);
		variantTypesByName.put(TYPENAME_MULTIALLELIC_SNV, TYPE_MULTIALLELIC_SNV);
		variantTypesByName.put(TYPENAME_EMBEDDED_SNV, TYPE_EMBEDDED_SNV);
		variantTypesByName.put(TYPENAME_INDEL, TYPE_INDEL);
		variantTypesByName.put(TYPENAME_STR, TYPE_STR);
		variantTypesByName.put(TYPENAME_CNV, TYPE_CNV);
		variantTypesByName.put(TYPENAME_REPEAT, TYPE_REPEAT);
		variantTypesByName.put(TYPENAME_LARGEDEL, TYPE_LARGEDEL);
		variantTypesByName.put(TYPENAME_LARGEINS, TYPE_LARGEINS);
		variantTypesByName.put(TYPENAME_INVERSION, TYPE_INVERSION);
		//Backwards compatibility
		variantTypesByName.put("Deletion", TYPE_LARGEDEL);
		variantTypesByName.put("Insertion", TYPE_LARGEINS);
		
		variantTypesById = new HashMap<Byte, String>(10);
		variantTypesById.put(TYPE_BIALLELIC_SNV, TYPENAME_BIALLELIC_SNV);
		variantTypesById.put(TYPE_MULTIALLELIC_SNV, TYPENAME_MULTIALLELIC_SNV);
		variantTypesById.put(TYPE_EMBEDDED_SNV, TYPENAME_EMBEDDED_SNV);
		variantTypesById.put(TYPE_INDEL, TYPENAME_INDEL);
		variantTypesById.put(TYPE_STR, TYPENAME_STR);
		variantTypesById.put(TYPE_CNV, TYPENAME_CNV);
		variantTypesById.put(TYPE_REPEAT, TYPENAME_REPEAT);
		variantTypesById.put(TYPE_LARGEDEL, TYPENAME_LARGEDEL);
		variantTypesById.put(TYPE_LARGEINS, TYPENAME_LARGEINS);
		variantTypesById.put(TYPE_INVERSION, TYPENAME_INVERSION);
	}
	
	
	
	public static String getVariantTypeName (byte variantTypeId) {
		if(variantTypesById==null) buildVariantTypeMaps();
		return variantTypesById.get(variantTypeId);
	}
	
	private static final List<String> DEFAULT_ALLELES_BIALLELIC = buildDefaultAllelesBiallelic();
	private static List<String> buildDefaultAllelesBiallelic() {
		List<String> alleles = new ArrayList<String>();
		alleles.add("REF");
		alleles.add("ALT");
		return Collections.unmodifiableList(alleles);
	}

	public static byte getVariantTypeId (String variantTypeName) {
		if(variantTypesByName==null) buildVariantTypeMaps();
		if(variantTypeName==null) return TYPE_UNDETERMINED;
		Byte id = variantTypesByName.get(variantTypeName);
		if(id == null) return TYPE_UNDETERMINED;
		return id;
	}
	
	
	public GenomicVariantImpl(String sequenceName, int first, List<String> alleles) {
		this.setSequenceName(sequenceName);
		this.setFirst(first);
		this.length = alleles.get(0).length();
		int last = first+length-1;
		this.setLast(last);
		this.setAlleles(alleles);
	}
	public GenomicVariantImpl(String sequenceName, int first, int last, List<String> alleles) {
		this.setSequenceName(sequenceName);
		this.setFirst(first);
		this.length = alleles.get(0).length();
		this.setLast(last);
		this.setAlleles(alleles);
	}
	/**
	 * Creates a biallelic genomic variant with default allele names
	 * @param sequenceName
	 * @param first
	 * @param last
	 * @param type Type of the variant
	 */
	public GenomicVariantImpl(String sequenceName, int first, int last, byte type) {
		this.setSequenceName(sequenceName);
		this.setFirst(first);
		this.setLast(last);
		this.setType(type);
		this.setLength(last-first+1);
		alleles = DEFAULT_ALLELES_BIALLELIC;
	}
	public String getSequenceName() {
		return sequenceName;
	}
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}

	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}

	public void setAlleles(List<String> alleles) {
		if(alleles.size()==0) throw new IllegalArgumentException("Alelle list can not be empty");
		this.alleles.clear();
		for(String allele:alleles) {
			addAllele(allele);
		}
	}
	public void addAllele(String allele) {
		alleles.add(allele.toUpperCase());
	}
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
	@Override
	public int length() {
		return length;
	}
	@Override
	public int getFirst() {
		return first;
	}
	@Override
	public int getLast() {
		return last;
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
	public String[] getAlleles() {
		return alleles.toArray(new String[0]);
	}
	@Override
	public String getReference() {
		return alleles.get(0);
	}
	@Override
	public short getVariantQS() {
		return variantQS;
	}
	@Override
	public byte getType() {
		return type;
	}
	@Override
	public void setType(byte type) {
		this.type = type;
	}
	/**
	 * Changes the quality score of the variant
	 * @param variantQS New variant quality score
	 */
	public void setVariantQS(short variantQS) {
		this.variantQS = variantQS;
	}
	
	
	
	public void setFirst(int first) {
		this.first = first;
	}

	public void setLast(int last) {
		this.last = last;
	}
	
	/**
	 * Changes the length of the genomic variant
	 * @param length New length
	 */
	public void setLength(int length) {
		this.length = length;
	}



	@Override
	public boolean isCompatible(GenomicVariant variant) {
		return testCompatibility(this, variant,true,true);
	}
	public static boolean testCompatibility(GenomicVariant v1, GenomicVariant v2, boolean testLocation, boolean testAlternative) {
		if (v1 == v2) return true;
		if(testLocation) {
			if(v1.getFirst() != v2.getFirst()) return false;
			if(v1.getLast() != v2.getLast()) return false;
			if(!v1.getSequenceName().equals(v2.getSequenceName())) return false;
		}
		if(testAlternative) {
			
			
			List<String> v1Alleles = Arrays.asList(v1.getAlleles());
			Collections.sort(v1Alleles);
			List<String> v2Alleles = Arrays.asList(v2.getAlleles());
			Collections.sort(v2Alleles);
			if(v1Alleles.size()!=v2Alleles.size()) return false;
			for(int i=0;i<v1Alleles.size();i++) {
				if(!v1Alleles.get(i).equals(v2Alleles.get(i))) return false;
			}
		} else if(!v1.getReference().equals(v2.getReference())) return false;
		return true;
	}
	@Override
	public boolean isBiallelic() {
		return alleles.size()==2;
	}
	@Override
	public boolean isSNV() {
		if(alleles.size()<2) return false;
		for(String allele:alleles) {
			if(allele.length()!=1) return false;
			if(!DNASequence.isInAlphabeth(allele.charAt(0))) return false;
		}
		return true;
	}

	/**
	 * Merges two informations about type of variant by priority, having STRs larger priority than indels and SNVs
	 * @param c1 First type
	 * @param c2 Second type
	 * @return byte Merged type
	 */
	public static byte mergeType(byte c1, byte c2) {
		return (byte)Math.max(c1, c2);
	}
	
}
