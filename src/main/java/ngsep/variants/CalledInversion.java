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

public class CalledInversion extends GenomicRegionImpl  {
	public static final String TYPE_INVERSION ="Inversion";
	private String id;
	private String sampleId;
	private int numAlns = 0;
	private double genotypeProbability = 0;
	
	public CalledInversion (String seqName, int first, int last, int numAlns) {
		super(seqName,first,last);
		this.numAlns = numAlns;
		
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
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
		return GenomicVariant.TYPE_INVERSION;
	}
	
	/*public String getSource() {
		return ImpreciseCalledGenomicVariant.SOURCE_READPAIR;
	}*/
}