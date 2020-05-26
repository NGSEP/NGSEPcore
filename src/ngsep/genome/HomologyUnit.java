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
public class HomologyUnit implements GenomicRegion {
	private String id;
	//TODO: Change to entityID, more generic? GenomeID when built for a genome, OrganismID when built for organisms
	private int genomeId;
	private String sequenceName;
	private int first;
	private int last;
	private boolean negativeStrand = false;
	private String unitSequence;
	
	//Homolog relationships
	private Map<Integer, Map<String,HomologyEdge>> homologsMap = new HashMap<>();
	private int totalHomologs = 0;
	
	//LCS mates of this orthology unit
	private Map<Integer, HomologyUnit> matesInLCS = new HashMap<>();
	
	public HomologyUnit(int genomeId, String id, String sequenceName, int first, int last) {
		super();
		this.genomeId = genomeId;
		this.id = id;
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		homologsMap.put(genomeId, new HashMap<>());
	}
	public HomologyUnit(int genomeId, String id, String unitSequence) {
		this.genomeId = genomeId;
		this.id = id;
		this.unitSequence = unitSequence;
		homologsMap.put(genomeId, new HashMap<>());
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
	 * Adds a new ortholog to this unit
	 * @param unit to associate as ortholog
	 */
	public void addHomologRelationship (HomologyEdge edge) {
		HomologyUnit subjectUnit = edge.getSubjectUnit();
		int genomeId = subjectUnit.getGenomeId();
		String unitId = subjectUnit.getId();
		Map<String,HomologyEdge> unitsGenome = homologsMap.computeIfAbsent(genomeId, V -> new HashMap<String, HomologyEdge>());
		if(!unitsGenome.containsKey(unitId)) totalHomologs++;
		unitsGenome.put(unitId, edge);
	}
	
	
	public Collection<HomologyEdge> getAllHomologyRelationships() {
		ArrayList<HomologyEdge> relationships = new ArrayList<>();
		for(Integer id : homologsMap.keySet()) {
			Map<String, HomologyEdge> genomeEdges = homologsMap.get(id);
			for(String geneId : genomeEdges.keySet()) {
				relationships.add(genomeEdges.get(geneId));
			}
		}
		return relationships;
	}
	
	public void removeAllHomologyRelationships() {
		homologsMap = new HashMap<Integer, Map<String,HomologyEdge>>();
		homologsMap.put(genomeId, new HashMap<String, HomologyEdge>());
	}

	/**
	 * @return the paralogs
	 */
	public Collection<HomologyEdge> getParalogRelationships() {
		Map<String,HomologyEdge> map = homologsMap.get(genomeId);
		Collection<HomologyEdge> collection = map.values();
		return Collections.unmodifiableCollection(collection);
	}

	/**
	 * Searches ortholog relationships to the given genome id
	 * @param id of the genome to query
	 * @return Homolog relationships to units in the given genome id
	 */
	public Collection<HomologyEdge> getOrthologRelationships(int genomeId) {
		Map<String,HomologyEdge> orthologs = homologsMap.get(genomeId);
		if(orthologs==null) return new ArrayList<>();
		return Collections.unmodifiableCollection(orthologs.values());
	}
	
	public Collection<HomologyUnit> getOrthologUnits(int genomeId) {
		List<HomologyUnit> answer = new ArrayList<>();
		Map<String,HomologyEdge> orthologRelationships = homologsMap.get(genomeId);
		if(orthologRelationships==null) return answer;
		for(HomologyEdge edge:orthologRelationships.values()) {
			answer.add(edge.getSubjectUnit());
		}
		return answer;
	}
	
	/**
	 * Returns the ortholog of this unit with the given genome id only if it is unique
	 * @param genomeId Id of the genome to query
	 * @return OrthologyUnit ortholog of this unit or null if there are not orthologs or if the ortholog is not unique
	 */
	public HomologyUnit getUniqueOrtholog(int genomeId) {
		Collection<HomologyEdge> orthologRelationships = getOrthologRelationships(genomeId);
		if(orthologRelationships.size()!=1) return null;
		return orthologRelationships.iterator().next().getSubjectUnit();
		
	}
	
	/**
	 * Returns the ortholog of this unit with the given genome id and the given id
	 * @param genomeId Id of the genome to query
	 * @param id of the unit to look for
	 * @return OrthologyUnit ortholog of this unit or null if the ortholog was not found
	 */
	public HomologyUnit getOrtholog(int genomeId, String id) {
		Map<String,HomologyEdge> orthologRelationships = homologsMap.get(genomeId);
		if(orthologRelationships==null) return null;
		HomologyEdge edge = orthologRelationships.get(id);
		if(edge == null) return null;
		return edge.getSubjectUnit();
		
	}
	/**
	 * Tells if the given unit is an ortholog of this unit
	 * @param unit to search
	 * @return boolean true if the given ortholog was found in the list of orthologs of this unit 
	 */
	public boolean isOrtholog (HomologyUnit unit) {
		return getOrtholog(unit.getGenomeId(), unit.getId())!=null;
	}

	/**
	 * @return the total number of orthologs
	 */
	public int getTotalHomologs() {
		return totalHomologs;
	}
	
	public boolean isUnique() {
		return homologsMap.get(genomeId).size()==0;
	}
	
	public HomologyUnit getLCSMate (int genomeId) {
		return matesInLCS.get(genomeId);
	}
	
	public void setMateInLCS (HomologyUnit mate) {
		matesInLCS.put(mate.genomeId,mate);
	}
	public String getUniqueKey() {
		return ""+genomeId+"\t"+id;
	}
}