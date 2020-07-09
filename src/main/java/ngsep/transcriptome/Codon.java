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
/**
 * Information for a codon that is translated into an aminoacid
 * @author Jorge Duitama
 *
 */
public class Codon {
	private String rnaSequence;
	private char aminoacid;
	private boolean start;
	private boolean stop;
	/**
	 * Creates a codon with the given information
	 * @param rnaSequence Sequence of RNA that becomes the given aminoacid
	 * @param aminoacid Aminoacid after the codon translation
	 */
	public Codon (String rnaSequence, char aminoacid) {
		this(rnaSequence,aminoacid,false,false);
	}
	/**
	 * Creates a codon with the given information
	 * @param rnaSequence Sequence of RNA that becomes the given aminoacid
	 * @param aminoacid Aminoacid after the codon translation
	 * @param start true if the codon is start for translation
	 * @param stop true if the codon is stop for translation
	 */
	public Codon (String rnaSequence, char aminoacid, boolean start, boolean stop) {
		this.setRnaSequence(rnaSequence);
		this.setAminoacid(aminoacid);
		this.setStart(start);
		this.setStop(stop);
	}
	/**
	 * @return String Sequence of RNA that becomes the given aminoacid
	 */
	public String getRnaSequence() {
		return rnaSequence;
	}
	/**
	 * Changes the RNA sequence
	 * @param rnaSequence New RNA sequence
	 */
	public void setRnaSequence(String rnaSequence) {
		this.rnaSequence = rnaSequence;
	}
	/**
	 * @return char Aminoacid after this codon translation
	 */
	public char getAminoacid() {
		return aminoacid;
	}
	/**
	 * Changes the aminoacid
	 * @param aminoacid new aminoacid
	 */
	public void setAminoacid(char aminoacid) {
		this.aminoacid = aminoacid;
	}
	/**
	 * @return boolean true if this codon is start for translation
	 */
	public boolean isStart() {
		return start;
	}
	/**
	 * Sets the start status
	 * @param start New status
	 */
	public void setStart(boolean start) {
		this.start = start;
	}
	/**
	 * @return boolean true if this codon is stop for translation
	 */
	public boolean isStop() {
		return stop;
	}
	/**
	 * Sets the stop status
	 * @param stop New status
	 */
	public void setStop(boolean stop) {
		this.stop = stop;
	}
}
