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
package ngsep.sequences;


/**
 * Implementation of AbstractSequence for DNA Sequences masked with lower case characters to 
 * flag repeats
 * @author Jorge Duitama
 */
public class DNAMaskedSequence extends AbstractLimitedSequence{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static final String BASES = "AaCcNngGtT";
	public static final String [] BASES_ARRAY = {"A","a","C","c","N","n","g","G","t","T"};
	/**
	 * Creates a new DNA masked sequence with the given bases
	 * @param sequence String with the bases for the new sequence
	 */
	public DNAMaskedSequence (CharSequence sequence) {
		this.setSequence(sequence);
	}
	/**
	 * Creates a new empty DNA Masked Sequence
	 */
	public DNAMaskedSequence () {
		
	}
	@Override
	public String getAlphabet() {
		return BASES;
	}
	@Override
	protected int getBitsPerCharacter() {
		return 4;
	}
	@Override
	protected int getDefaultIndex(){
		return 4;
	}
	/**
	 * @return DNAMaskedSequence Reverse complement of this sequence
	 */
	public DNAMaskedSequence getReverseComplement() {
		StringBuffer sb = new StringBuffer();
		DNAMaskedSequence answer = new DNAMaskedSequence();
		int l = length();
		for(int i=l-1;i>=0;i--) {
			sb.append(getComplement(charAt(i)));
			if(sb.length()%100000==0) {
				answer.append(sb);
				sb =new StringBuffer();
			}
		}
		answer.append(sb);
		return answer;
	}
	/**
	 * Gets the complement of the given base
	 * @param base to complement
	 * @return char Complement of the given base
	 */
	public static char getComplement(char base) {
		int index = BASES.indexOf(base);
		if(index>=0 && base != 'N' && base !='n') {
			return BASES.charAt(BASES.length()-index-1);
		}
		return base;
	}
	/**
	 * Returns the reverse complement of a DNA sequence encoded in the given String
	 * @param sequence Sequence to be translated
	 * @return CharSequence reverse complement of the given sequence
	 */
	public static CharSequence getReverseComplement(CharSequence sequence) {
		if(sequence instanceof DNASequence) {
			return ((DNASequence)sequence).getReverseComplement();
		}
		if(sequence instanceof DNAMaskedSequence) {
			return ((DNAMaskedSequence)sequence).getReverseComplement();
		}
		return new DNAMaskedSequence(sequence).getReverseComplement();
	}
}
