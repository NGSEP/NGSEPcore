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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Daniel Tello
 * @author Jorge Duitama
 */
public class OrthologyUnit implements GenomicRegion {
	private String id;
	private int genomeId;
	private String sequenceName;
	private int first;
	private int last;
	private boolean negativeStrand = false;
	private String unitSequence;
	//Orthologs of other genomes
	private Map<Integer, Map<String,OrthologyUnit>> orthologsMap = new HashMap<>();
	private int totalOrthologs = 0;
	
	//LCS mates of this orthology unit
	private Map<Integer, OrthologyUnit> matesInLCS = new HashMap<>();
	
	public OrthologyUnit(int genomeId, String id, String sequenceName, int first, int last) {
		super();
		this.genomeId = genomeId;
		this.id = id;
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		orthologsMap.put(genomeId, new HashMap<>());
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
	 * @return the unitSequence
	 */
	public String getUnitSequence() {
		return unitSequence;
	}


	/**
	 * @param unitSequence the unitSequence to set
	 */
	public void setUnitSequence(String unitSequence) {
		this.unitSequence = unitSequence;
	}


	/**
	 * @return the paralogs
	 */
	public Collection<OrthologyUnit> getParalogs() {
		return Collections.unmodifiableCollection(orthologsMap.get(genomeId).values());
	}
	
	/**
	 * Adds a new paralog to this unit
	 * @param unit to associate as paralog
	 */
	public void addParalog (OrthologyUnit unit) {
		if(unit.getGenomeId()!=genomeId) throw new RuntimeException("Can not add a paralog from a different genome. This genome id: "+genomeId+" unit id: "+unit.getGenomeId());
		addOrtholog(unit);
	}


	/**
	 * Searches orthologs from the given genome id
	 * @param id of the genome to query
	 * @return Orthologs from the given id
	 */
	public Collection<OrthologyUnit> getOrthologs(int genomeId) {
		Map<String,OrthologyUnit> orthologs = orthologsMap.get(genomeId);
		if(orthologs==null) return new ArrayList<>();
		return Collections.unmodifiableCollection(orthologs.values());
	}
	
	/**
	 * Get all orthologs from genomes different than that of this orthology unit
	 * @return List<OrthologyUnit> 
	 */
	public List<OrthologyUnit> getOrthologsOtherGenomes() {
		List<OrthologyUnit> answer = new ArrayList<>();
		for(int id:orthologsMap.keySet()) {
			if(id!=genomeId) answer.addAll(orthologsMap.get(id).values());
		}
		return answer;
	}
	/**
	 * Get all orthologs from all genomes including the genome of this unit
	 * @return List<OrthologyUnit> 
	 */
	public List<OrthologyUnit> getOrthologsAllGenomes() {
		List<OrthologyUnit> answer = new ArrayList<>();
		for(int id:orthologsMap.keySet()) {
			answer.addAll(orthologsMap.get(id).values());
		}
		return answer;
	}
	
	/**
	 * Returns the ortholog of this unit with the given genome id only if it is unique
	 * @param genomeId Id of the genome to query
	 * @return OrthologyUnit ortholog of this unit or null if there are not orthologs or if the ortholog is not unique
	 */
	public OrthologyUnit getUniqueOrtholog(int genomeId) {
		Collection<OrthologyUnit> orthologs = getOrthologs(genomeId);
		if(orthologs.size()!=1) return null;
		return orthologs.iterator().next();
		
	}
	
	/**
	 * Returns the ortholog of this unit with the given genome id and the given id
	 * @param genomeId Id of the genome to query
	 * @param id of the unit to look for
	 * @return OrthologyUnit ortholog of this unit or null if the ortholog was not found
	 */
	public OrthologyUnit getOrtholog(int genomeId, String id) {
		Map<String,OrthologyUnit> orthologs = orthologsMap.get(genomeId);
		if(orthologs==null) return null;
		return orthologs.get(id);
		
	}
	/**
	 * Tells if the given unit is an ortholog of this unit
	 * @param unit to search
	 * @return boolean true if the given ortholog was found in the list of orthologs of this unit 
	 */
	public boolean isOrtholog (OrthologyUnit unit) {
		return getOrtholog(unit.getGenomeId(), unit.getId())!=null;
	}

	/**
	 * @return the total number of orthologs
	 */
	public int getTotalOrthologs() {
		return totalOrthologs;
	}


	/**
	 * Adds a new ortholog to this unit
	 * @param unit to associate as ortholog
	 */
	public void addOrtholog (OrthologyUnit unit) {
		int genomeId = unit.getGenomeId();
		Map<String,OrthologyUnit> unitsGenome = orthologsMap.get(genomeId);
		if(unitsGenome == null) {
			unitsGenome = new HashMap<>();
			orthologsMap.put(genomeId, unitsGenome);
		}
		if(!unitsGenome.containsKey(unit.getId())) totalOrthologs++;
		unitsGenome.put(unit.getId(),unit);
	}
	
	public boolean isUnique() {
		return orthologsMap.get(genomeId).size()==0;
	}
	
	public OrthologyUnit getLCSMate (int genomeId) {
		return matesInLCS.get(genomeId);
	}
	
	public void setMateInLCS (OrthologyUnit mate) {
		matesInLCS.put(mate.genomeId,mate);
	}
	public String getUniqueKey() {
		return ""+genomeId+"\t"+id;
	}
}