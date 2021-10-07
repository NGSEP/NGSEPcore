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
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * @author Laura Gonzalez
 */
public class HomologyCluster 
{

	// Constants for default values
	public final static String SCORE = "soft_core";
	public final static String SACCESORY = "soft_accesory";
	public final static String ECORE = "exact_core";
	public final static String EACCESORY = "exact_accesory";
	public final static String PERSISTENT = "persistent";
	public final static String SHELL = "shell";
	public final static String CLOUD = "cloud";
	
	// Attributes
	private List<HomologyUnit> homologyUnitsCluster = new ArrayList<>();
	private Map<Integer,GenomicRegionSortedCollection<LocalHomologyCluster>> homologyUnitsByRegion = new HashMap<>();
	private int clusterId;
	private double frequency;
	private String exactCat;
	private String softCat;
	
	public HomologyCluster(int clusterId, List<HomologyUnit> homologyUnitsCluster) {
		this.clusterId = clusterId;
		this.homologyUnitsCluster = homologyUnitsCluster;
		Map<Integer,GenomicRegionSortedCollection<HomologyUnit>> rawUnitsByRegion = new HashMap<>();
		for(HomologyUnit unit:homologyUnitsCluster) {
			GenomicRegionSortedCollection<HomologyUnit> genomeUnits = rawUnitsByRegion.computeIfAbsent(unit.getGenomeId(), v->new GenomicRegionSortedCollection<>());
			genomeUnits.add(unit);
		}
		for(Map.Entry<Integer,GenomicRegionSortedCollection<HomologyUnit>> entry:rawUnitsByRegion.entrySet()) {
			int genomeId = entry.getKey();
			List<HomologyUnit> current = entry.getValue().asList(); 
			GenomicRegionSortedCollection<LocalHomologyCluster> merged = new GenomicRegionSortedCollection<>();
			LocalHomologyCluster nextCluster = null;
			for(HomologyUnit unit:current) {
				if(nextCluster==null || nextCluster.getSequenceName()!=unit.getSequenceName() || nextCluster.getLast()<unit.getFirst()+50000) {
					nextCluster = new LocalHomologyCluster(this,unit);
					merged.add(nextCluster);
				} else {
					nextCluster.addUnit(unit);
				}
			}
			homologyUnitsByRegion.put(genomeId, merged);
		}
	}
	
	public List<HomologyUnit> getHomologyUnitsCluster() {
		return this.homologyUnitsCluster;
	}
	
	public int getClusterId() {
		return this.clusterId;
	}
	
	public double getFrequency() {
		return this.frequency;
	}
	
	public void setFrequency(double frequency) {
		this.frequency=frequency;
	}
	
	public String getExactCategory() {
		return this.exactCat;
	}
	
	public void setExactCategory(String exactCat) {
		this.exactCat = exactCat;
	}
	
	public String getSoftCategory() {
		return this.softCat;
	}
	
	public void setSoftCategory(String softCat) {
		this.softCat = softCat;
	}

	public HomologyUnit findUniqueHomologyUnit(int genomeId) {
		GenomicRegionSortedCollection<LocalHomologyCluster> subclustersGenome = homologyUnitsByRegion.get(genomeId);
		if(subclustersGenome==null) return null;
		//TODO: define better
		return subclustersGenome.asList().get(0).getHomologyUnitsCluster().get(0);
	}
	public List<LocalHomologyCluster> getLocalClusters(int genomeId, String sequenceName) {
		List<LocalHomologyCluster> answer = new ArrayList<LocalHomologyCluster>();
		GenomicRegionSortedCollection<LocalHomologyCluster> subclustersGenome = homologyUnitsByRegion.get(genomeId);
		if(subclustersGenome==null) return answer;
		return subclustersGenome.getSequenceRegions(sequenceName).asList();
	}
	public LocalHomologyCluster getLocalCluster(HomologyUnit unit) {
		GenomicRegionSortedCollection<LocalHomologyCluster> subclustersGenome = homologyUnitsByRegion.get(unit.getGenomeId());
		if(subclustersGenome==null) return null;
		GenomicRegionSortedCollection<LocalHomologyCluster> intersecting = subclustersGenome.findSpanningRegions(unit);
		if(intersecting.size()==0) return null;
		return intersecting.asList().get(0);
	}
}
