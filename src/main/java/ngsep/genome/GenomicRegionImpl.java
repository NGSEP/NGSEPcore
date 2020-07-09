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

/**
 * Basic implementation of genomic region useful for searching
 * @author Jorge Duitama
 *
 */
public class GenomicRegionImpl implements GenomicRegion {
	private String sequenceName;
	private int first;
	private int last;
	private boolean negativeStrand = false;
	
	public GenomicRegionImpl(String sequenceName, int first, int last) {
		super();
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
	}

	@Override
	public String getSequenceName() {
		return sequenceName;
	}

	@Override
	public int getFirst() {
		return first;
	}
	@Override
	public int getLast() {
		return last;
	}

	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}
	
	public void setFirst(int first) {
		this.first = first;
	}

	public void setLast(int last) {
		this.last = last;
	}

	@Override
	public int length() {
		return last - first +1;
	}

	@Override
	public boolean isPositiveStrand() {
		return !negativeStrand;
	}

	@Override
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	
	public void setNegativeStrand (boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
	
	

}
