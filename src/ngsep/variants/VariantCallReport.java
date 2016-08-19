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
package ngsep.variants;

import java.util.Map;
import java.util.TreeMap;

public class VariantCallReport {
	private Map<String, Integer> allelesMap;
	private int [] counts;
	private double [][] logConditionals;
	public VariantCallReport(String[] alleles, int[] counts, double[][] logCond) {
		super();
		if(alleles == null) throw new IllegalArgumentException("Alleles in a variant report can not be null");
		allelesMap = new TreeMap<String, Integer>();
		for(int i=0;i<alleles.length;i++) {
			allelesMap.put(alleles[i], i);
		}
		this.counts = counts;
		logConditionals = logCond;
	}
	public int getCount(String allele) {
		if(!countsPresent()) return 0;
		Integer i = allelesMap.get(allele);
		if(i==null) return 0;
		return counts[i];
	}
	public double getLogConditionalProbability (String allele1, String allele2) {
		if(!logConditionalsPresent()) return 0;
		Integer i = allelesMap.get(allele1);
		Integer j = allelesMap.get(allele2);
		if(i==null || j==null) return 0;
		return logConditionals[i][j];
	}
	public boolean countsPresent() {
		return counts!=null;
	}
	public boolean logConditionalsPresent() {
		return logConditionals!=null;
	}
}
