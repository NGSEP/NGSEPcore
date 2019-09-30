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
 * Interface for sequences with a small defined character set
 * @author Jorge Duitama
 */
public interface LimitedSequence extends CharSequence,Comparable<LimitedSequence> {
	
	public static final char GAP_CHARACTER = '-';
	/**
	 * Makes this sequence equals to the given sequence
	 * @param sequence Sequence with data to copy
	 */
	public void setSequence(CharSequence sequence);
	/**
	 * Append the given CharSequence at the end of this
	 * @param sequence Sequence to append
	 */
	public void append(CharSequence sequence);
	/**
	 * Sets the given character at the given index
	 * @param index Zero based position to modify
	 * @param ch Character to set at index
	 */
	public void setCharAt(int index, char ch);
	/**
	 * Returns a String with the characters allowed in this sequence
	 * @return String
	 */
	public String getAlphabet();
	/**
	 * Returns the character of the alphabet with the given index in the alphabet String.
	 * @param index Zero based index to look for in the alphabet
	 * @return char Character at the given index
	 */
	public char getAlphabetCharacter(int index);
	/**
	 * Gets the index for the given character in the alphabet.
	 * @param ch character to look for
	 * @return in Zero based index or -1 if the character is not in the alphabet
	 */
	public int getAlphabetIndex(char ch);
	/**
	 * Tells if the character is in the alphabet
	 * @param ch character to look for
	 * @return boolean true if ch is in the alphabet, false otherwise
	 */
	public boolean isInAlphabet(char ch);
	/**
	 * Returns the size of the alphabet
	 * @return int Size of the alphabet
	 */
	public int getAlphabetSize();
	
}
