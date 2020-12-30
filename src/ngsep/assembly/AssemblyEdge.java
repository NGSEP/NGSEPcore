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
public class AssemblyEdge {
	
	private AssemblyVertex vertex1;
	private AssemblyVertex vertex2;
	private int overlap;
	private int averageOverlap;
	private int medianOverlap;
	private int fromLimitsOverlap;
	private int overlapStandardDeviation;
	private int numSharedKmers;
	private int coverageSharedKmers;
	private int weightedCoverageSharedKmers;
	private int rawKmerHits = 0;
	private int rawKmerHitsSubjectStartSD = 0;
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
	
	public int getNumSharedKmers() {
		return numSharedKmers;
	}

	public void setNumSharedKmers(int numSharedKmers) {
		this.numSharedKmers = numSharedKmers;
	}

	public int getOverlapStandardDeviation() {
		return overlapStandardDeviation;
	}

	public void setOverlapStandardDeviation(int overlapStandardDeviation) {
		this.overlapStandardDeviation = Math.round(overlapStandardDeviation);
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

	public boolean isLayoutEdge() {
		return layoutEdge;
	}

	public void setLayoutEdge(boolean layoutEdge) {
		this.layoutEdge = layoutEdge;
	}
	
	public String toString() {
		return System.lineSeparator()+"v1 "+getVertex1()+" v2: "+getVertex2()+" overlap: "+getOverlap()+" coverage shared kmers: "+getCoverageSharedKmers()+" weighted cov: "+getWeightedCoverageSharedKmers()+" layout: "+layoutEdge;
	}
	

}
