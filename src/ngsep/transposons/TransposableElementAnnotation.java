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
	private String queryName;
	private TransposableElementFamily family;
	private String taxonomy;
	private double count = 0;
	private boolean negativeStrand=false;
	
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
	
	public String getQueryName() {
		return queryName;
	}

	public TransposableElementFamily getFamily() {
		return family;
	}

	public String getTaxonomy() {
		return taxonomy;
	}
	
	public void setQueryName(String queryName) {
		this.queryName = queryName;
	}

	public void setTaxonomy(String taxonomy) {
		this.taxonomy = taxonomy;
	}

	public void setSourceInfo(String sourceInfo) {
		int i = sourceInfo.indexOf('#');
		if(i<0) {
			queryName = sourceInfo;
			return;
		}
		queryName = sourceInfo.substring(0,i);
		taxonomy = sourceInfo.substring(i+1);
		i=taxonomy.indexOf('/');
		if(i<0) i=taxonomy.length();
		String orderStr = taxonomy.substring(0,i);
		if (i<taxonomy.length()) {
			String familyInfo2 = taxonomy.substring(i+1);
			i=familyInfo2.indexOf('/');
			if(i<0) i=familyInfo2.length();
			String familyStr = familyInfo2.substring(0,i);
			family = TransposableElementFamily.findFamily(orderStr, familyStr);
		} else {
			family = TransposableElementFamily.findUnknown(orderStr);
		}
	}
	public String getSequenceName() {
		return sequenceName;
	}
	public int getFirst() {
		return first;
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
	public double getCount() {
		return count;
	}
	public void setCount(double count) {
		this.count = count;
	}
}
