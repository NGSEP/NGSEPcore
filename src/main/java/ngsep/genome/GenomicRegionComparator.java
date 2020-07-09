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
package ngsep.genome;

import java.util.Comparator;

import ngsep.sequences.QualifiedSequenceList;

/**
 * Comparator for two genomic regions. Compares sequence names according with a given sorted list
 * of possible sequence names.
 * @author Jorge Duitama
 */
public class GenomicRegionComparator implements Comparator<GenomicRegion> {
	private QualifiedSequenceList sequences;
	
	public GenomicRegionComparator(QualifiedSequenceList sequences) {
		this.sequences = sequences; 
	}
	public QualifiedSequenceList getSequences() {
		return sequences;
	}
	/**
	 * Compare the two regions. In general it follows the contract of the compare method in Comparator, which means that it 
	 * returns a negative number if r1 goes before r2 and a positive number if r1 goes after r2. However, this comparator
	 * provides additional information with the following codes (in absolute value):
	 * 3: Regions are located in different sequences
	 * 2: Regions are located in the same sequence but they do not overlap
	 * 1: Regions are located in the same sequence and they overlap
	 * To improve efficiency avoiding String comparisons, if sequence indexes are set for the genomic
	 * regions, it assumes that these indexes are consistent with the indexes in the sequences list 
	 */
	@Override
	public int compare(GenomicRegion r1, GenomicRegion r2) {
		int p1 = sequences.indexOf(r1.getSequenceName());
		int p2 = sequences.indexOf(r2.getSequenceName());
		if(p1 < p2) return -3;
		if (p1>p2) return 3;
		boolean overlap = GenomicRegionSpanComparator.getInstance().span(r1, r2);
		if(r1.getFirst()<r2.getFirst()) {
			return overlap?-1:-2;
		}
		if(r2.getFirst()<r1.getFirst()) {
			return overlap?1:2;
		}
		//In these two cases the regions for sure overlap
		if(r1.getLast()<r2.getLast()) return -1;
		if(r2.getLast()<r1.getLast()) return 1;
		return 0;
	}
}
