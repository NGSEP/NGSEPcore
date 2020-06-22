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

import java.io.Serializable;

import ngsep.sequences.KmerHitsCluster;

/**
 * @author Jorge Duitama
 */
public class AssemblyEmbedded implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1120528617112863153L;
	
	private int sequenceId;
	private CharSequence read;
	private boolean isReverse;
	private int hostId;
	private int startPosition;
	private KmerHitsCluster evidence;
	private int coverageSharedKmers;
	private int mismatches;
	

	public AssemblyEmbedded(int sequenceId, CharSequence read, boolean isReverse, int hostId, int startPosition) {
		this.sequenceId = sequenceId;
		this.read = read;
		this.isReverse = isReverse;
		this.hostId = hostId;
		this.startPosition = startPosition;
	}

	
	/**
	 * @return the sequenceId
	 */
	public int getSequenceId() {
		return sequenceId;
	}


	public CharSequence getRead() {
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


	public int getStartPosition() {
		return startPosition;
	}

	public void setStartPosition(int startPosition) {
		this.startPosition = startPosition;
	}


	public KmerHitsCluster getEvidence() {
		return evidence;
	}

	public void setEvidence(KmerHitsCluster evidence) {
		this.evidence = evidence;
	}
	
	
	public int getCoverageSharedKmers() {
		return coverageSharedKmers;
	}


	public void setCoverageSharedKmers(int coverageSharedKmers) {
		this.coverageSharedKmers = coverageSharedKmers;
	}


	public int getMismatches() {
		return mismatches;
	}


	public void setMismatches(int mismatches) {
		this.mismatches = mismatches;
	}


	public String toString () {
		return ""+sequenceId+"_"+isReverse+"_"+hostId+"_"+startPosition;
	}
	
}
