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
package ngsep.discovery.rd;

import ngsep.genome.GenomicRegionImpl;

public class ReadDepthBin extends GenomicRegionImpl {
	
	private double gcContent;
	private double rawReadDepth=0;
	private double correctedReadDepth=0;
	private double readDepthLevel = 0;
	private boolean inRepetitiveRegion = false;
	public ReadDepthBin(String sequenceName, int first, int last,double gcContent) {
		super(sequenceName, first, last);
		this.gcContent = gcContent;
	}
	public double getGcContent() {
		return gcContent;
	}
	public void setGcContent(double gcContent) {
		this.gcContent = gcContent;
	}
	public double getRawReadDepth() {
		return rawReadDepth;
	}
	public void setRawReadDepth(double rawReadDepth) {
		this.rawReadDepth = rawReadDepth;
	}
	public double getCorrectedReadDepth() {
		return correctedReadDepth;
	}
	public void setCorrectedReadDepth(double correctedReadDepth) {
		this.correctedReadDepth = correctedReadDepth;
	}
	public void addRead() {
		rawReadDepth++;
	}
	public double getReadDepthLevel() {
		return readDepthLevel;
	}
	public void setReadDepthLevel(double level) {
		this.readDepthLevel = level;
	}
	public boolean isInRepetitiveRegion() {
		return inRepetitiveRegion;
	}
	public void setInRepetitiveRegion(boolean inRepetitiveRegion) {
		this.inRepetitiveRegion = inRepetitiveRegion;
	}
	public boolean isGoodForAverage () {
		return !inRepetitiveRegion && gcContent>=0;
	}
}
