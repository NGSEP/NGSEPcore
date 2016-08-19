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
 * Implementation of AbstractSequence for DNA Sequences without degenerate bases
 * @author Jorge Duitama
 */
public class DNASequence extends AbstractLimitedSequence{
	public static final String [] BASES_ARRAY = {"A","C","G","T"};
	public static final String BASES_STRING = "ACGT";
	/**
	 * Creates a new DNA sequence with the given bases
	 * @param sequence String with the bases for the new sequence
	 */
	public DNASequence(String sequence) {
		this.setSequence(sequence);
	}
	/**
	 * Creates a new empty DNA Sequence
	 */
	public DNASequence() {
	}
	@Override
	public String getAlphabet() {
		return BASES_STRING;
	}
	@Override
	protected int getBitsPerCharacter() {
		return 2;
	}
	public static boolean isInAlphabeth(char base) {
		return BASES_STRING.indexOf(base)>=0;
	}
	
	
	/**
	 * @return DNASequence Reverse complement of this sequence
	 */
	public DNASequence getReverseComplement() {
		StringBuffer sb = new StringBuffer();
		DNASequence answer = new DNASequence();
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
		int index = BASES_STRING.indexOf(base);
		if(index>=0) {
			return BASES_STRING.charAt(BASES_STRING.length()-index-1);
		}
		return base;
	}
	/**
	 * Returns the reverse complement of a DNA sequence encoded in the given String
	 * @param sequence Sequence to be translated
	 * @return String reverse complement of the given sequence
	 */
	public static String getReverseComplement(String sequence) {
		return new DNASequence(sequence).getReverseComplement().toString();
	}
}
