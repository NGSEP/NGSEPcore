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

import ngsep.genome.GenomicRegion;

/**
 * @author Jorge Duitama
 */
public class TransposableElementAnnotation implements GenomicRegion{
	private String sequenceName;
	private int first;
	private int last;
	private TransposableElement source;
	private TransposableElementFamily inferredFamily;
	private int count = 0;
	private boolean negativeStrand=false;
	private boolean bordersFixed = false;
	
	/**
	 * Build New Transposable Element 
	 * @param sequenceName
	 * @param first
	 * @param last
	 */
	public TransposableElementAnnotation(String sequenceName, int first, int last) {
		super();
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
	}

	public String getSequenceName() {
		return sequenceName;
	}
	public int getFirst() {
		return first;
	}
	
	public void setFirst(int first) {
		this.first = first;
	}

	public int getLast() {
		return last;
	}
	public void setLast(int last) {
		this.last = last;
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return !negativeStrand;
	}
	@Override
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
	
	public TransposableElement getSource() {
		return source;
	}

	public void setSource(TransposableElement source) {
		this.source = source;
	}
	
	public TransposableElementFamily getInferredFamily() {
		return inferredFamily;
	}

	public void setInferredFamily(TransposableElementFamily inferredFamily) {
		this.inferredFamily = inferredFamily;
	}

	public String getQueryName() {
		return source.getId();
	}

	public TransposableElementFamily getSourceFamily() {
		return source.getFamily();
	}

	public String getTaxonomy() {
		return source.getTaxonomy();
	}
	public int getCount() {
		return count;
	}
	public void setCount(int count) {
		this.count = count;
	}

	public boolean isBordersFixed() {
		return bordersFixed;
	}

	public void setBordersFixed(boolean bordersFixed) {
		this.bordersFixed = bordersFixed;
	}

	public boolean isValidated() {
		return isBordersFixed() && getSourceFamily()!=null;
	}
	
}
