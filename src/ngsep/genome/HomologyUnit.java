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

import ngsep.math.ShannonEntropyCalculator;
import ngsep.sequences.AbstractLimitedSequence;
import ngsep.sequences.AminoacidSequence;

/**
 * @author Daniel Tello
 * @author Jorge Duitama
 */
public class HomologyUnit implements GenomicRegion {
	private String id;
	//TODO: Change to entityID, more generic? GenomeID when built for a genome, OrganismID when built for organisms
	private int genomeId;
	private String uniqueKey;
	private String sequenceName;
	private int first;
	private int last;
	private boolean negativeStrand = false;
	private CharSequence unitSequence;
	private CharSequence cdsSequence;
	private ShannonEntropyCalculator entropyCalculator = new ShannonEntropyCalculator(0);
	
	//Homolog relationships
	private Map<Integer, Map<String,HomologyEdge>> homologsMap = new HashMap<>();
	private int totalHomologs = 0;
	
	public HomologyUnit(int genomeId, String id, String sequenceName, int first, int last) {
		super();
		this.genomeId = genomeId;
		this.id = id;
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		homologsMap.put(genomeId, new HashMap<>());
		uniqueKey = ""+genomeId+"\t"+id;
	}
	public HomologyUnit(int genomeId, String id, CharSequence unitSequence) {
		this.genomeId = genomeId;
		this.id = id;
		this.unitSequence = unitSequence;
		homologsMap.put(genomeId, new HashMap<>());
		uniqueKey = ""+genomeId+"\t"+id;
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
	public CharSequence getUnitSequence() {
		return unitSequence;
	}


	/**
	 * @param unitSequence the unitSequence to set
	 */
	public void setUnitSequence(CharSequence unitSequence) {
		this.unitSequence = unitSequence;
	}
	
	
	
	public CharSequence getCdsSequence() {
		return cdsSequence;
	}
	public void setCdsSequence(CharSequence cdsSequence) {
		this.cdsSequence = cdsSequence;
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
		totalHomologs = 0;
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
	 * Searches homolog relationships to the given genome id
	 * @param id of the genome to query
	 * @return Homolog relationships to units in the given genome id
	 */
	public Collection<HomologyEdge> getHomologRelationships(int genomeId) {
		Map<String,HomologyEdge> homologs = homologsMap.get(genomeId);
		if(homologs==null) return new ArrayList<>();
		return Collections.unmodifiableCollection(homologs.values());
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
	 * Returns the homolog of this unit with the given genome id only if it is unique
	 * @param genomeId Id of the genome to query
	 * @return OrthologyUnit homolog of this unit or null if there are not homologs or if the homolog is not unique
	 */
	public HomologyUnit getUniqueHomolog(int genomeId) {
		Collection<HomologyEdge> homologRelationships = getHomologRelationships(genomeId);
		if(homologRelationships.size()!=1) return null;
		return homologRelationships.iterator().next().getSubjectUnit();
		
	}
	
	/**
	 * Returns the ortholog of this unit with the given genome id and the given id
	 * @param genomeId Id of the genome to query
	 * @param id of the unit to look for
	 * @return HomologyUnit ortholog of this unit or null if the ortholog was not found
	 */
	public HomologyUnit getOrtholog(int genomeId, String id) {
		Map<String,HomologyEdge> orthologRelationships = homologsMap.get(genomeId);
		if(orthologRelationships==null) return null;
		HomologyEdge edge = orthologRelationships.get(id);
		if(edge == null) return null;
		return edge.getSubjectUnit();
		
	}
	
	/**
	 * Returns the edge connecting this unit with the given unit
	 * @param unit Homology unit related to this unit
	 * @return HomologyEdge Relationship between units
	 */
	public HomologyEdge getHomologyEdge(HomologyUnit unit) {
		Map<String,HomologyEdge> orthologRelationships = homologsMap.get(unit.getGenomeId());
		if(orthologRelationships==null) return null;
		HomologyEdge edge = orthologRelationships.get(unit.getId());
		return edge;
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
	
	public String getUniqueKey() {
		return uniqueKey;
	}
	public Map<Long,Double> getKmerCodesWithEntropies(int kmerLength,int kmerOffset) {
		Map<Long,Double> answer = new HashMap<>();
		for(int i=0; i<unitSequence.length()-kmerLength+1; i+=kmerOffset) { 
			long kmerCode = AbstractLimitedSequence.getHash(unitSequence, i, i+kmerLength, (AbstractLimitedSequence)unitSequence);
			String kmer = new String(AbstractLimitedSequence.getSequence(kmerCode, kmerLength, AminoacidSequence.EMPTY_AA_SEQUENCE));
			double entropy = entropyCalculator.calculateEntropy(kmer);
			//System.out.println("code: "+kmerCode+" kmer: "+kmer+" weight: "+weight );
			answer.put(kmerCode,entropy);
		}
		return answer;
	}
}