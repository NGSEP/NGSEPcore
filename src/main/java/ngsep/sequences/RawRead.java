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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Represents a raw sequencing read including an id, the sequence and quality scores
 * @author Jorge Duitama
 */
public class RawRead extends QualifiedSequence {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public RawRead(String id, CharSequence sequence, String qualityScores) {
		super(id,sequence,qualityScores);
	}
	/**
	 * Loads a raw read from the given BufferedReader
	 * @param in buffer to read
	 * @return RawRead
	 * @throws IOException if the buffer can not be read
	 * @deprecated See ngsep.sequences.FastqFileReader
	 */
	public static RawRead load(BufferedReader in) throws IOException {
		String id = in.readLine();
		if(id==null) return null;
		String seq = in.readLine();
		if(seq==null) return null;
		String plus = in.readLine();
		if(plus==null) return null;
		String qs = in.readLine();
		if(qs==null) return null;
		return new RawRead(id, seq, qs);
	}
	public String getSequenceString() {
		return getCharacters().toString();
	}
	public void trimFirstNucleotides (int n) {
		setCharacters(getSequenceString().substring(n));
		setQualityScores(getQualityScores().substring(n));
	}
	/**
	 * Trims the end of the read starting from the first appearance of the given pattern 
	 * @param regexp Sequence to trim reads
	 * @return boolean true if the read was trimmed
	 */
	public boolean trimFromSequence(Pattern regexp) {
		String sequence = getSequenceString();
		if(regexp==null) return false;
		Matcher m = regexp.matcher(sequence);
		if(m.find()) {
			int start = m.start();
			trimToLength(start);
			return true;
		}
		return false;
	}
	public void trimToLength(int newLength) {
		String sequence = getSequenceString();
		if(newLength<sequence.length()) {
			setCharacters(sequence.substring(0,newLength));
			setQualityScores(getQualityScores().substring(0,newLength));
		}
		
	}
	public void save (PrintStream out) {
		out.println("@"+this.getName());
		out.println(this.getCharacters());
		out.println("+");
		out.println(this.getQualityScores());
	}
	
	private static Map<Character,String> largeFixedQSStrings = new HashMap<Character, String>();
	
	public static final String generateFixedQSString(char qs, int length) {
		String charsStr = largeFixedQSStrings.get(qs);
		if(charsStr==null || charsStr.length()<length) {
			char [] qualities = new char[Math.max(length,500000)];
			Arrays.fill(qualities, qs);
			charsStr = new String(qualities);
			largeFixedQSStrings.put(qs, charsStr);
		}
		return charsStr.substring(0,length);
	}
}
