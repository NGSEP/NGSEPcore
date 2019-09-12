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
 * Implementation of AbstractSequence for DNA Sequences with degenerate bases
 * @author Jorge Duitama
 *
 */
public class DegenerateSequence extends AbstractLimitedSequence {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static final String bases = "ACGT";
	public static final String [] iubCodes2 = {"AMRW","MCSY","RSGK","WYKT"};
	public static final String iubCodes3 = "BDHV";
	public static final String alphabet = "ACMRBDWNSHVYKGT";
	//W and S appear two times because they are self-complementary
	private static final String alphaForComplement = "ACMRBDWSNSWHVYKGT";
	
	/**
	 * Creates a new degenerate sequence with the given bases
	 * @param sequence String with the bases for the new sequence
	 */
	public DegenerateSequence (String sequence) {
		this.setSequence(sequence);
	}
	/**
	 * Creates a new empty degenerate sequence
	 */
	public DegenerateSequence() {
		
	}
	
	@Override
	public String getAlphabet() {
		return alphabet;
	}
	@Override
	protected int getBitsPerCharacter() {
		return 4;
	}
	@Override
	protected int getDefaultIndex(){
		return 7;
	}
	/**
	 * Given a degenerate base calculate the bases that it represents
	 * @param base Degenerate base
	 * @return String The same base if it is not degenerate. The two or three bases concatenated
	 * if it is a degenerate base or null if it is not found
	 */
	public static String getExtendedBases(char base) {
		if(base == 'N') return bases;
		int index = bases.indexOf(base); 
		if(index>=0) {
			return ""+base;
		}
		index = iubCodes3.indexOf(base);
		if(index>=0) {
			StringBuffer answer = new StringBuffer();
			for(int i=0;i<bases.length();i++) {
				if(i!=index) {
					answer.append(bases.charAt(i));
				}
			}
			return answer.toString();
		}
		for(int i=0;i<iubCodes2.length;i++) {
			index = iubCodes2[i].indexOf(base);
			if(index>=0) {
				return""+bases.charAt(i)+""+bases.charAt(index);
			}
		}
		return null;
	}
	
	/**
	 * @return String Reverse complement of this sequence as a String object
	 */
	public DegenerateSequence getReverseComplement() {
		StringBuffer sb = new StringBuffer();
		DegenerateSequence answer = new DegenerateSequence();
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
	 * Gets the complement of the given base, which can be degenerate
	 * @param base to complement
	 * @return char Complement of the given base
	 */
	public static char getComplement(char base) {
		int index = alphaForComplement.indexOf(base);
		if(index>=0) {
			return alphaForComplement.charAt(alphaForComplement.length()-index-1);
		}
		return base;
	}
	
	/**
	 * Returns the reverse complement of a Degenerate sequence encoded in the given String
	 * @param sequence Sequence to be translated
	 * @return String reverse complement of the given sequence
	 */
	public static String getReverseComplement(String sequence) {
		return new DegenerateSequence(sequence).getReverseComplement().toString();
	}
	
	/**
	 * Returns the degenerate base corresponding with the two non degenerate given bases
	 * @param base1 First non degenerate base
	 * @param base2 Second non degenerate base
	 * @return char Degenerate base corresponding with the given two bases
	 */
	public static char getDegenerateBase(char base1, char base2) {
		int i = bases.indexOf(base1);
		int j = bases.indexOf(base2);
		return iubCodes2[i].charAt(j);
	}
	/**
	 * Builds a java regular expression for the given degenerate sequence 
	 * @param degenerateSequence Sequence to build a regular expression
	 * @return String regular expression representing the sequence
	 */
	public static String makeRegularExpression(String degenerateSequence) {
		StringBuilder regexp = new StringBuilder();
		for(int i=0;i<degenerateSequence.length();i++) {
			char degBase = degenerateSequence.charAt(i);
			String extBases = DegenerateSequence.getExtendedBases(degBase);
			if(extBases.length()==0) throw new IllegalArgumentException("Unrecognized degenerate base "+degBase+" in sequence "+degenerateSequence);
			if(extBases.length()==1) regexp.append(extBases);
			else {
				regexp.append("(");
				for(int j=0;j<extBases.length();j++) {
					if(j>0) regexp.append("|");
					regexp.append(extBases.charAt(j));
				}
				regexp.append(")");
			}
		}
		return regexp.toString();
	}
}
