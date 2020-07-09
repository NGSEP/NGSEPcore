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

/**
 * Comparator for genomic regions which makes two regions equal if they overlap.
 * Assumes that genomic regions to compare are located in the same sequence
 * @author Jorge Duitama
 *
 */
public class GenomicRegionSpanComparator implements Comparator<GenomicRegion> {
	private static GenomicRegionSpanComparator instance = new GenomicRegionSpanComparator();
	private GenomicRegionSpanComparator() {
	}
	
	@Override
	public int compare(GenomicRegion r1, GenomicRegion r2) {
		return compare(r1,r2.getFirst(),r2.getLast());
	}
	public int compare(GenomicRegion r1, int first, int last) {
		if (span(r1,first,last)) {
			return 0;
		}
		return r1.getFirst()-first;
	}
	public boolean span(GenomicRegion r1, GenomicRegion r2) {
		return span(r1.getFirst(),r1.getLast(),r2.getFirst(),r2.getLast());
	}
	public boolean span(GenomicRegion r, int first, int last) {
		return span(r.getFirst(),r.getLast(),first,last);
	}
	public boolean span(int first1, int last1, int first2, int last2) {
		return first1 <= last2 && first2 <= last1;
	}
	public int getSpanLength(int first1, int last1, int first2, int last2) {
		if(!span(first1,last1,first2,last2)) return 0;
		int span = Math.min(last1-first2+1, last2-first1+1);
		span = Math.min(span, last1-first1+1);
		span = Math.min(span, last2-first2+1);
		return span;
	}
	public static GenomicRegionSpanComparator getInstance() {
		return instance;
	}
}
