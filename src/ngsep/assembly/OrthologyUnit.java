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
package ngsep.assembly;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.genome.GenomicRegion;

/**
 * @author Daniel Tello
 * @author Jorge Duitama
 */
class OrthologyUnit implements GenomicRegion {
	private String id;
	private int genomeId;
	private String sequenceName;
	private int first;
	private int last;
	private boolean negativeStrand = false;
	private String proteinSequence;
	private List<OrthologyUnit> paralogs = new ArrayList<>();
	//Orthologs of other genomes
	private Map<Integer, List<OrthologyUnit>> orthologsMap = new HashMap<>();
	
	//Genomes for which this orthology unit is in LCS
	private Set<Integer> genomesInLCS = new HashSet<>();
	
	public OrthologyUnit(int genomeId, String id, String sequenceName, int first, int last) {
		super();
		this.genomeId = genomeId;
		this.id = id;
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
	}
	

	public String getId() {
		return id;
	}

	/**
	 * @return the genomeId
	 */
	public int getGenomeId() {
		return genomeId;
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
	/**
	 * @param negativeStrand the negativeStrand to set
	 */
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
	/**
	 * @return the proteinSequence
	 */
	public String getProteinSequence() {
		return proteinSequence;
	}
	/**
	 * @param proteinSequence the proteinSequence to set
	 */
	public void setProteinSequence(String proteinSequence) {
		this.proteinSequence = proteinSequence;
	}


	/**
	 * @return the paralogs
	 */
	public List<OrthologyUnit> getParalogs() {
		return paralogs;
	}
	
	/**
	 * Adds a new paralog to this unit
	 * @param unit to associate as paralog
	 */
	public void addParalog (OrthologyUnit unit) {
		if(unit.getGenomeId()!=genomeId) throw new RuntimeException("Can not add a paralog from a different genome. This genome id: "+genomeId+" unit id: "+unit.getGenomeId());
		paralogs.add(unit);
	}


	/**
	 * Searches orthologs from the given genome id
	 * @param id of the genome to query
	 * @return Orthologs from the given id
	 */
	public List<OrthologyUnit> getOrthologs(int genomeId) {
		List<OrthologyUnit> orthologs = orthologsMap.get(genomeId);
		if(orthologs==null) orthologs = new ArrayList<>();
		return orthologs;
	}
	
	/**
	 * Returns the ortholog of this unit with the given genome id only if it is unique
	 * @param genomeId Id of the genome to query
	 * @return OrthologyUnit ortholog of this unit or null if there are not orthologs or if the ortholog is not unique
	 */
	public OrthologyUnit getUniqueOrtholog(int genomeId) {
		List<OrthologyUnit> orthologs = getOrthologs(genomeId);
		if(orthologs.size()!=1) return null;
		return orthologs.get(0);
		
	}
	
	/**
	 * Adds a new ortholog to this unit
	 * @param unit to associate as ortholog
	 */
	public void addOrtholog (OrthologyUnit unit) {
		int genomeId = unit.getGenomeId();
		List<OrthologyUnit> unitsGenome = orthologsMap.get(genomeId);
		if(unitsGenome == null) {
			unitsGenome = new ArrayList<>();
			orthologsMap.put(genomeId, unitsGenome);
		}
		unitsGenome.add(unit);
	}
	
	public boolean isUnique() {
		return paralogs.size()==0;
	}
	
	public boolean isInLCS (int genomeId) {
		return genomesInLCS.contains(genomeId);
	}
	
	public void setInLCS (int genomeId) {
		genomesInLCS.add(genomeId);
	}
}