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

/**
 * @author Jorge Duitama
 */
public class AssemblyEdge implements AssemblySequencesRelationship {
	
	private AssemblyVertex vertex1;
	private AssemblyVertex vertex2;
	private int overlap;
	private int averageOverlap;
	private int medianOverlap;
	private int fromLimitsOverlap;
	private int overlapStandardDeviation;
	private int rawKmerHits = 0;
	private int rawKmerHitsSubjectStartSD = 0;
	private int numSharedKmers;
	private int coverageSharedKmers;
	private int weightedCoverageSharedKmers;
	private int numIndels = 0;
	private int vertex1EvidenceStart;
	private int vertex1EvidenceEnd;
	private int vertex2EvidenceStart;
	private int vertex2EvidenceEnd;
	private int numMismatches = -1;
	private boolean layoutEdge = false;
	

	public AssemblyEdge(AssemblyVertex vertex1, AssemblyVertex vertex2, int overlap) {
		this.vertex1 = vertex1;
		this.vertex2 = vertex2;
		this.overlap = overlap;
	}

	/**
	 * @return the vertex1
	 */
	public AssemblyVertex getVertex1() {
		return vertex1;
	}

	/**
	 * @return the vertex2
	 */
	public AssemblyVertex getVertex2() {
		return vertex2;
	}
	
	/**
	 * @return the overlap
	 */
	public int getOverlap() {
		return overlap;
		//return fromLimitsOverlap;
		//return averageOverlap;
	}

	/**
	 * @param overlap the overlap to set
	 */
	public void setOverlap(int overlap) {
		this.overlap = overlap;
	}
	
	

	public int getAverageOverlap() {
		return averageOverlap;
	}

	public void setAverageOverlap(int averageOverlap) {
		this.averageOverlap = averageOverlap;
	}

	public int getMedianOverlap() {
		return medianOverlap;
	}

	public void setMedianOverlap(int medianOverlap) {
		this.medianOverlap = medianOverlap;
	}

	public int getFromLimitsOverlap() {
		return fromLimitsOverlap;
	}

	public void setFromLimitsOverlap(int fromLimitsOverlap) {
		this.fromLimitsOverlap = fromLimitsOverlap;
	}

	/**
	 * @return the cost
	 */
	public int getCost() {
		int l1 = vertex1.getRead().getLength();
		int l2 = vertex2.getRead().getLength();
		if(isSameSequenceEdge()) return l1;
		int cost = l1 + l2;
		int toSubstract = Math.min(l1, l2)-1;
		toSubstract = Math.min(toSubstract, overlap);
		cost-= toSubstract;
		return cost;
	}

	public int getOverlapStandardDeviation() {
		return overlapStandardDeviation;
	}

	public void setOverlapStandardDeviation(int overlapStandardDeviation) {
		this.overlapStandardDeviation = Math.round(overlapStandardDeviation);
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
	
	public double getIndelsPerKbp () {
		return 1000.0*numIndels / (double)overlap;
	}

	public int getVertex1EvidenceStart() {
		return vertex1EvidenceStart;
	}

	public void setVertex1EvidenceStart(int vertex1EvidenceStart) {
		this.vertex1EvidenceStart = vertex1EvidenceStart;
	}

	public int getVertex1EvidenceEnd() {
		return vertex1EvidenceEnd;
	}

	public void setVertex1EvidenceEnd(int vertex1EvidenceEnd) {
		this.vertex1EvidenceEnd = vertex1EvidenceEnd;
	}

	public int getVertex2EvidenceStart() {
		return vertex2EvidenceStart;
	}

	public void setVertex2EvidenceStart(int vertex2EvidenceStart) {
		this.vertex2EvidenceStart = vertex2EvidenceStart;
	}

	public int getVertex2EvidenceEnd() {
		return vertex2EvidenceEnd;
	}

	public void setVertex2EvidenceEnd(int vertex2EvidenceEnd) {
		this.vertex2EvidenceEnd = vertex2EvidenceEnd;
	}
	
	public int getNumMismatches() {
		return numMismatches;
	}

	public void setNumMismatches(int numMismatches) {
		this.numMismatches = numMismatches;
	}

	public AssemblyVertex getConnectingVertex(AssemblyVertex vertex) {
		if(vertex1==vertex) return vertex2;
		if(vertex2==vertex) return vertex1;
		return null;
	}
	
	public AssemblyVertex getSharedVertex (AssemblyEdge edge2) {
		if(vertex1==edge2.vertex1 || vertex1 == edge2.vertex2) return vertex1;
		if(vertex2==edge2.vertex1 || vertex2 == edge2.vertex2) return vertex2;
		return null;
	}
	
	public boolean isSameSequenceEdge() {
		return vertex1.getSequenceIndex() == vertex2.getSequenceIndex();
	}
	
	public double calculateEvidenceProportion() {
		double evidenceProp = vertex1EvidenceEnd-vertex1EvidenceStart;
		evidenceProp += vertex2EvidenceEnd-vertex2EvidenceStart;
		evidenceProp/=(2*overlap);
		if(evidenceProp>1) evidenceProp = 2 - evidenceProp;
		return evidenceProp;
	}

	public boolean isLayoutEdge() {
		return layoutEdge;
	}

	public void setLayoutEdge(boolean layoutEdge) {
		this.layoutEdge = layoutEdge;
	}
	
	public String toString() {
		return System.lineSeparator()+"v1 "+getVertex1()+" v2: "+getVertex2()+" overlap: "+getOverlap()+" CSK: "+getCoverageSharedKmers()+" WCSK: "+getWeightedCoverageSharedKmers()+" Ev1: "+vertex1EvidenceStart+" "+vertex1EvidenceEnd+" "+((double)(vertex1EvidenceEnd-vertex1EvidenceStart)/(overlap+1))+" Ev2: "+vertex2EvidenceStart+" "+vertex2EvidenceEnd+" "+((double)(vertex2EvidenceEnd-vertex2EvidenceStart)/(overlap+1))+" Indels: "+numIndels+" Mismatches: "+numMismatches+" layout: "+layoutEdge;
	}
	

}
