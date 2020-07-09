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


import ngsep.genome.GenomicRegionImpl;


public class CalledLargeIndel extends GenomicRegionImpl {
	public static final String TYPE_INSERTION ="Insertion";
	public static final String TYPE_DELETION ="Deletion";
	private String id;
	private String sampleId;
	private int length = 0;
	private int numAlns = 0;
	private boolean deletion = true;
	private double genotypeProbability = 0;
	private int numSplitReads = 0;
	
	public CalledLargeIndel(String seqName, int first, int last, boolean deletion, int length, int numAlns) {
		super(seqName,first,last);
		this.length = length;
		this.numAlns = numAlns;
		this.deletion = deletion;
	}
	
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	
	public boolean isDeletion() {
		return deletion;
	}
	public int getNumSplitReads() {
		return numSplitReads;
	}
	public void setNumSplitReads(int numSplitReads) {
		this.numSplitReads = numSplitReads;
	}
	
	public String getSampleId() {
		return sampleId;
	}

	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}

	public double getGenotypeProbability() {
		return genotypeProbability;
	}
	public void setGenotypeProbability(double genotypeProbability) {
		this.genotypeProbability = genotypeProbability;
	}
	
	public int getSupportingFragments() {
		return numAlns;
	}
	
	public byte getType() {
		if(deletion) return GenomicVariant.TYPE_LARGEDEL;
		return GenomicVariant.TYPE_LARGEINS;
	}
	@Override
	public int length() {
		return length;
	}
	/*public String getSource() {
		return ImpreciseCalledGenomicVariant.SOURCE_READPAIR;
	}*/
}
