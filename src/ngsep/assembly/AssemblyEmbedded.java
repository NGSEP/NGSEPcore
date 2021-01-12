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
package ngsep.assembly;

import ngsep.sequences.QualifiedSequence;

/**
 * @author Jorge Duitama
 */
public class AssemblyEmbedded {
	
	private int sequenceId;
	private QualifiedSequence read;
	private boolean isReverse;
	private int hostId;
	private int hostStart;
	private int hostEnd;
	private int hostStartStandardDeviation;
	private int rawKmerHits = 0;
	private int rawKmerHitsSubjectStartSD = 0;
	private int numSharedKmers;
	private int coverageSharedKmers;
	private int weightedCoverageSharedKmers;
	private int numIndels = 0;
	private int hostEvidenceStart;
	private int hostEvidenceEnd;
	private int sequenceEvidenceStart;
	private int sequenceEvidenceEnd;
	private int numMismatches = -1;
	

	public AssemblyEmbedded(int sequenceId, QualifiedSequence read, boolean isReverse, int hostId, int hostStart, int hostEnd) {
		this.sequenceId = sequenceId;
		this.read = read;
		this.isReverse = isReverse;
		this.hostId = hostId;
		this.hostStart = hostStart;
		this.hostEnd = hostEnd;
	}

	
	/**
	 * @return the sequenceId
	 */
	public int getSequenceId() {
		return sequenceId;
	}

	public QualifiedSequence getRead() {
		return read;
	}

	public boolean isReverse() {
		return isReverse;
	}
	public void setReverse(boolean isReverse) {
		this.isReverse = isReverse;
	}
	
	public int getHostId() {
		return hostId;
	}

	public void setHostId(int hostId) {
		this.hostId = hostId;
	}

	public int getHostStart() {
		return hostStart;
	}
	public void setHostStart(int hostStart) {
		this.hostStart = hostStart;
	}

	public int getHostEnd() {
		return hostEnd;
	}
	public void setHostEnd(int hostEnd) {
		this.hostEnd = hostEnd;
	}
	
	public int getHostStartStandardDeviation() {
		return hostStartStandardDeviation;
	}
	public void setHostStartStandardDeviation(int hostStartStandardDeviation) {
		this.hostStartStandardDeviation = hostStartStandardDeviation;
	}
	public int getRawKmerHits() {
		return rawKmerHits;
	}

	public void setRawKmerHits(int rawKmerHits) {
		this.rawKmerHits = rawKmerHits;
	}

	public int getRawKmerHitsSubjectStartSD() {
		return rawKmerHitsSubjectStartSD;
	}
	public void setRawKmerHitsSubjectStartSD(int rawKmerHitsSubjectStartSD) {
		this.rawKmerHitsSubjectStartSD = rawKmerHitsSubjectStartSD;
	}
	
	public int getNumSharedKmers() {
		return numSharedKmers;
	}
	public void setNumSharedKmers(int numSharedKmers) {
		this.numSharedKmers = numSharedKmers;
	}

	public int getCoverageSharedKmers() {
		return coverageSharedKmers;
	}
	public void setCoverageSharedKmers(int coverageSharedKmers) {
		this.coverageSharedKmers = coverageSharedKmers;
	}
	
	public int getWeightedCoverageSharedKmers() {
		return weightedCoverageSharedKmers;
	}
	public void setWeightedCoverageSharedKmers(int weightedCoverageSharedKmers) {
		this.weightedCoverageSharedKmers = weightedCoverageSharedKmers;
	}
	
	public int getNumIndels() {
		return numIndels;
	}
	public void setNumIndels(int numIndels) {
		this.numIndels = numIndels;
	}

	public int getHostEvidenceStart() {
		return hostEvidenceStart;
	}
	public void setHostEvidenceStart(int hostEvidenceStart) {
		this.hostEvidenceStart = hostEvidenceStart;
	}

	public int getHostEvidenceEnd() {
		return hostEvidenceEnd;
	}
	public void setHostEvidenceEnd(int hostEvidenceEnd) {
		this.hostEvidenceEnd = hostEvidenceEnd;
	}
	
	public int getSequenceEvidenceStart() {
		return sequenceEvidenceStart;
	}
	public void setSequenceEvidenceStart(int sequenceEvidenceStart) {
		this.sequenceEvidenceStart = sequenceEvidenceStart;
	}

	public int getSequenceEvidenceEnd() {
		return sequenceEvidenceEnd;
	}
	public void setSequenceEvidenceEnd(int sequenceEvidenceEnd) {
		this.sequenceEvidenceEnd = sequenceEvidenceEnd;
	}
	
	public int getNumMismatches() {
		return numMismatches;
	}

	public void setNumMismatches(int numMismatches) {
		this.numMismatches = numMismatches;
	}


	public double calculateEvidenceProportion () {
		double evidenceProp = hostEvidenceEnd-hostEvidenceStart;
		evidenceProp += (sequenceEvidenceEnd-sequenceEvidenceStart);
		evidenceProp/=(2*read.getLength());
		if(evidenceProp>1) evidenceProp = 2 - evidenceProp;
		return evidenceProp;
	}

	public String toString () {
		return ""+sequenceId+"_"+read.getName()+"_"+read.getLength()+" "+isReverse+"_"+hostId+"_"+hostStart+"_"+hostEnd;
	}
	
}
