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
package ngsep.alignments;

/**
 * Represents a paired read alignment
 * @author German Andrade
 * @author Jorge Duitama
 *
 */
public class ReadAlignmentPair{
	private ReadAlignment aln1;
	public ReadAlignment getAln1() {
		return aln1;
	}
	public ReadAlignment getAln2() {
		return aln2;
	}
	private ReadAlignment aln2;
	public ReadAlignmentPair(ReadAlignment pAln1, ReadAlignment pAln2) {
		aln1=pAln1;
		aln2=pAln2;
	}

}