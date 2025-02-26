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
 * 
 * @author Jorge Duitama
 *
 */
public class UngappedSearchHit {
	
	private int queryStart=-1;
	private int subjectIdx=-1;
	//Zero based start of this hit within the sequence
	private int subjectStart;
	private short hitLength;
	//Weight of this hit encoded as a byte
	private byte weight = 100;
	
	/**
	 * Creates a new ungapped hit to a subject sequence
	 * @param sequenceIdx id of the subject sequence
	 * @param subjectStart zero based first position of the hit within the subject
	 */
	public UngappedSearchHit(int subjectIdx, int subjectStart) {
		this.subjectIdx = subjectIdx;
		this.subjectStart = subjectStart;
	}
	
	public int getSubjectIdx() {
		return subjectIdx;
	}
	
	public int getSubjectStart() {
		return subjectStart;
	}
	
	public int getQueryStart() {
		return queryStart;
	}

	public void setQueryStart(int queryStart) {
		this.queryStart = queryStart;
	}
	
	public short getHitLength() {
		return hitLength;
	}

	public void setHitLength(short hitLength) {
		this.hitLength = hitLength;
	}

	public void setWeight (double weight) {
		if(weight>1) weight=1;
		if(weight<0) weight = 0;
		this.weight = (byte)Math.round(100.0*weight);
	}
	
	
	public double getWeight () {
		return 0.01*weight;
	}
	public int estimateSubjectStart() {
		return getSubjectStart() - getQueryStart();
	}
	public int estimateQueryStart() {
		return getQueryStart() - getSubjectStart();
	}
}
