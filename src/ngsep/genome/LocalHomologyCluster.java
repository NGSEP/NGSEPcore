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

import java.util.ArrayList;
import java.util.List;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class LocalHomologyCluster implements GenomicRegion {
	private List<HomologyUnit> units = new ArrayList<HomologyUnit>(); 
	private int genomeId;
	private HomologyCluster parent;
	private String seqName;
	private int first;
	private int last;
	public LocalHomologyCluster (HomologyCluster parent, HomologyUnit unit) {
		this.parent = parent;
		genomeId = unit.getGenomeId();
		units.add(unit);
		seqName = unit.getSequenceName();
		first = unit.getFirst();
		last = unit.getLast();
	}
	public void addUnit(HomologyUnit unit) {
		if(genomeId!=unit.getGenomeId()) throw new RuntimeException("Error merging local homology unit. Current genome id: "+genomeId+" new genome id: "+unit.getGenomeId());
		if(seqName!=unit.getSequenceName()) throw new RuntimeException("Error merging local homology unit. Current sequence name: "+seqName+" new sequence name: "+unit.getSequenceName());
		units.add(unit);
		first = Math.min(first, unit.getFirst());
		last = Math.max(last, unit.getLast());
	}
	
	
	public List<HomologyUnit> getUnits() {
		return units;
	}
	public int getGenomeId() {
		return genomeId;
	}
	public List<HomologyUnit> getHomologyUnitsCluster() {
		return units;
	}
	
	public HomologyCluster getParent() {
		return parent;
	}
	public void setParent(HomologyCluster parent) {
		this.parent = parent;
	}
	@Override
	public String getSequenceName() {
		return seqName;
	}
	@Override
	public int getFirst() {
		return first;
	}
	@Override
	public int getLast() {
		return last;
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return true;
	}
	@Override
	public boolean isNegativeStrand() {
		return false;
	}
}
