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
package ngsep.main.io;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class ParseUtils {
	public static final DecimalFormat ENGLISHFMT = new DecimalFormat("##0.0#",DecimalFormatSymbols.getInstance(Locale.ENGLISH));
	public static final DecimalFormat ENGLISHFMT_PROBABILITIES = new DecimalFormat("0.0###",DecimalFormatSymbols.getInstance(Locale.ENGLISH));
	
	/**
	 * Parse the given string sequence looking for a single delimiter
	 * @param s The string to parse
	 * @param delim Parsing delimiter
	 * @return String [] Array of strings with the identified tokens
	 */
	public static String[] parseString(String s,char delim) {
		//return s.split(""+delim);
		List<String> answer = new ArrayList<String>();
		char [] a = s.toCharArray();
		StringBuilder current = new StringBuilder(a.length);
		for (int i=0;i<a.length;i++) {
			if(a[i] == delim) {
				answer.add(current.toString());
				current = new StringBuilder(a.length);
				continue;
			}
			current.append(a[i]);
		}
		answer.add(current.toString());
		return answer.toArray(new String[0]);
	}
	/**
	 * Parse the given string with possible text that should not be parsed, 
	 * looking for a single delimiter 
	 * @param s The string to parse
	 * @param delim Parsing delimiter
	 * @param textDelimiter Character used to identify the beginning and end of texts that should not be parsed
	 * @return String [] Array of strings with the identified tokens
	 */
	public static String[] parseStringWithText(String s, char delim, char textDelimiter) {
		//return s.split(""+delim);
		List<String> answer = new ArrayList<String>();
		char [] a = s.toCharArray();
		StringBuilder current = new StringBuilder(a.length);
		boolean inText = false;
		for (int i=0;i<a.length;i++) {
			if(a[i] == textDelimiter) {
				inText = !inText;
			} else if( !inText && a[i] == delim) {
				answer.add(current.toString());
				current = new StringBuilder(a.length);
				continue;
			}
			current.append(a[i]);
		}
		answer.add(current.toString());
		return answer.toArray(new String[0]);
	}
	
	/**
	 * Parse the given string sequence looking for two delimiters
	 * @param s The string to parse
	 * @param delim1 First delimiter
	 * @param delim2 Second delimiter
	 * @return String [] Array of strings with the identified tokens 
	 */
	public static String[] parseString(String s,char delim1, char delim2) {
		//return s.split(""+delim);
		char [] a = s.toCharArray();
		List<String> answer = new ArrayList<String>();
		StringBuilder current = new StringBuilder(a.length);
		for (int i=0;i<a.length;i++) {
			if(a[i] == delim1 || a[i] == delim2) {
				answer.add(current.toString());
				current = new StringBuilder(a.length);
				continue;
			}
			current.append(a[i]);
		}
		answer.add(current.toString());
		return answer.toArray(new String[0]);
		
	}
}
