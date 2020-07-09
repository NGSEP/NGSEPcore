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
package ngsep.discovery;

public class PileupAlleleCall {
	private CharSequence sequence;
	private String qualityScores;
	private String readGroup;
	private boolean negativeStrand = false;
	
	public PileupAlleleCall(CharSequence sequence, String qualityScores) {
		super();
		this.sequence = sequence;
		this.qualityScores = qualityScores;
	}
	
	/**
	 * @return the alleleCall
	 */
	public CharSequence getSequence() {
		return sequence;
	}
	
	public String getAlleleString () {
		return sequence.toString().toUpperCase();
	}

	/**
	 * @return the qualityScores
	 */
	public String getQualityScores() {
		return qualityScores;
	}
	
	/**
	 * @return the readGroup
	 */
	public String getReadGroup() {
		return readGroup;
	}
	/**
	 * @param readGroup the readGroup to set
	 */
	public void setReadGroup(String readGroup) {
		this.readGroup = readGroup;
	}

	/**
	 * @return the negativeStrand
	 */
	public boolean isNegativeStrand() {
		return negativeStrand;
	}

	/**
	 * @param negativeStrand the negativeStrand to set
	 */
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
}
