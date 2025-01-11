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
package ngsep.transposons;

import ngsep.hmm.ProfileAlignmentDomain;
import ngsep.sequences.QualifiedSequence;

/**
 * @author Leonidas Villamil
 * @author Jorge Duitama
 */
public class TransposonDomainAlignment {
	private QualifiedSequence qseq;
	private int start;
	private int end;
	private boolean reverse = false;
	private ProfileAlignmentDomain alnDomain;
	public TransposonDomainAlignment(QualifiedSequence qseq, int start, int end, ProfileAlignmentDomain alnDomain) {
		this.qseq = qseq;
		this.start = start;
		this.end = end;
		this.alnDomain = alnDomain;
	}
	public boolean isReverse() {
		return reverse;
	}
	public void setReverse(boolean reverse) {
		this.reverse = reverse;
	}
	public QualifiedSequence getQseq() {
		return qseq;
	}
	public int getStart() {
		return start;
	}
	public int getEnd() {
		return end;
	}
	public int getLength() {
		return end-start;
	}
	public double getEvalue() {
		return alnDomain.getEvalue();
	}
	public String getHmmID() {
		return alnDomain.getHmmID();
	}
	public String getDomainCode() {
		return alnDomain.getDomainCode();
	}
}
