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

import java.util.HashMap;
import java.util.ResourceBundle;
/**
 * Class that implements translation from CDNA sequences to proteins
 * @author Jorge Duitama
 *
 */
public class ProteinTranslator {
	public static final String DEFAULT_BUNDLE="ngsep.transcriptome.ProteinTranslatorDefaultBundle";
	public HashMap<String, Codon> translateTable = new HashMap<String, Codon>();
	private static ProteinTranslator instance = new ProteinTranslator();
	/**
	 * Creates a new Protein translator with default translation
	 */
	public ProteinTranslator() {
		loadBundle(DEFAULT_BUNDLE);
	}
	/**
	 * Creates a new Protein translator translating according with the given resource information
	 * @param resourceName Name of the resource describing the rules for translation
	 */
	public ProteinTranslator(String resourceName) {
		loadBundle(resourceName);
	}
	/**
	 * @return ProteinTranslator Instance with default translation
	 */
	public static ProteinTranslator getInstance() {
		return instance;
	}
	/**
	 * Changes translation rules according with the given bundle
	 * @param resourceName Name of the resource describing the rules for translation
	 */
	public void loadBundle(String resourceName) {
		ResourceBundle bundle = ResourceBundle.getBundle(resourceName);
		String bases = "ACGU";
		for(int i=0;i<4;i++) {
			for(int j=0;j<4;j++) {
				for(int k=0;k<4;k++) {
					String codonStr =getCodonString(bases.charAt(i), bases.charAt(j), bases.charAt(k)); 
					String key="codon."+codonStr;
					String aminoacid=bundle.getString(key);
					boolean stop = "STOP".equalsIgnoreCase(aminoacid);
					boolean start = aminoacid.contains("START");
					Codon c;
					if(!stop) {
						c = new Codon(codonStr,aminoacid.charAt(0),start,stop);
					} else {
						c = new Codon(codonStr,'\0',start,stop);
					}
					translateTable.put(codonStr, c);
				}
			}
		}
	}
	/**
	 * Translates the given sequence
	 * @param sequence to translate
	 * @return String Corresponding aminoacid chain
	 */
	public String getProteinSequence(CharSequence sequence) {
		return new String(getProteinSequence(sequence.toString().toCharArray(), 0, false));
	}
	
	/**
	 * Translates the given sequence from the given start
	 * @param sequence to translate
	 * @param beginIndex First index to start translating
	 * @return String Corresponding aminoacid chain
	 */
	public String getProteinSequence(CharSequence sequence, int beginIndex) {
		return new String(getProteinSequence(sequence.toString().toCharArray(), beginIndex, false));
	}

	/**
	 * Translates the subsequence of the given sequence according with the given settings
	 * @param sequence Sequence to translate
	 * @param beginIndex First index to start translating
	 * @param lookForStartCodon If true, translation will start just when a start codon is observed
	 * @return char [] Aminoacid chain as a char array
	 */
	public char [] getProteinSequence(char [] sequence, int beginIndex, boolean lookForStartCodon) {
		return getProteinSequence(sequence, beginIndex, sequence.length, lookForStartCodon);
	}
	/**
	 * Translates the subsequence of the given sequence according with the given settings
	 * @param sequence Sequence to translate
	 * @param beginIndex First index to start translating
	 * @param endIndex Last index to finish translating
	 * @param lookForStartCodon If true, translation will start just when a start codon is observed
	 * @return char [] Aminoacid chain as a char array
	 */
	public char [] getProteinSequence(char [] sequence, int beginIndex, int endIndex, boolean lookForStartCodon) {
		int length = sequence.length;
		StringBuffer sb = new StringBuffer();
		boolean translate = !lookForStartCodon;
		for(int i=beginIndex;i<length-2 && i<endIndex-2;i+=3) {
			Codon c = getCodon(sequence[i],sequence[i+1],sequence[i+2]);
			if (c == null) {
				//System.err.println("Invalid codon: "+sequence[i]+sequence[i+1]+sequence[i+2]+" in sequence: "+String.valueOf(sequence));
				break;
			}
			if(c.isStop()) {
				break;
			}
			if(c.isStart()) {
				translate = true;
			}
			if(translate) {
				sb.append(c.getAminoacid());
			}
		}
		char [] answer = new char [sb.length()];
		sb.getChars(0, sb.length(), answer, 0);
		return answer;
	}
	/**
	 * Builds a codon for the given bases
	 * @param base1 First codon base in CDNA alphabet
	 * @param base2 Second codon base in CDNA alphabet
	 * @param base3 Third codon base in CDNA alphabet
	 * @return Codon Object with the codon information
	 */
	public Codon getCodon(char base1, char base2, char base3) {
		return translateTable.get(getCodonString(base1, base2, base3));
	}
	/**
	 * Calculates the aminoacid translation for the Codon made with the given bases
	 * @param base1 First codon base in CDNA alphabet
	 * @param base2 Second codon base in CDNA alphabet
	 * @param base3 Third codon base in CDNA alphabet
	 * @return char Aminoacid after translation
	 */
	public char getAminoacid(char base1, char base2, char base3) {
		
		Codon c = translateTable.get(getCodonString(base1, base2, base3));
		return c.getAminoacid();
	}
	/**
	 * Creates an codon RNA sequence from the given bases
	 * @param base1 First codon base in CDNA alphabet
	 * @param base2 Second codon base in CDNA alphabet
	 * @param base3 Third codon base in CDNA alphabet
	 * @return String codon sequence in RNA
	 */
	public String getCodonString(char base1, char base2, char base3) {
		base1=Character.toUpperCase(base1);
		base2=Character.toUpperCase(base2);
		base3=Character.toUpperCase(base3);
		if (base1 == 'T') base1='U';
		if (base2 == 'T') base2='U';
		if (base3 == 'T') base3='U';
		StringBuffer sb = new StringBuffer();
		sb.append(base1);
		sb.append(base2);
		sb.append(base3);
		return sb.toString();
	}
	

}
