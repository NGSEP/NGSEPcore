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
 * Comparator that assumes that both regions belong to the same sequence name. Compares the two
 * starts first and for equal starts compares the two ends. Since this comparator does not need
 * additional information, a static instance is kept following the Singleton pattern
 * @author Jorge Duitama
 */
public class GenomicRegionPositionComparator implements Comparator<GenomicRegion> {
	private static GenomicRegionPositionComparator instance = new GenomicRegionPositionComparator();
	private GenomicRegionPositionComparator() {
	}
	@Override
	public int compare(GenomicRegion r1, GenomicRegion r2) {
		if(r1.getFirst()!=r2.getFirst()) return r1.getFirst() - r2.getFirst();
		return r1.getLast() - r2.getLast();
	}
	public static GenomicRegionPositionComparator getInstance() {
		return instance;
	}
}
