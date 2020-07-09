package ngsep.genome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class HomologyCatalog {
	private Map<String, HomologyUnit> homologyUnitsMap= new HashMap<String, HomologyUnit>();
	private FMIndex indexHomologyUnits=null;
	
	public HomologyCatalog (List<HomologyUnit> units) {
		for(HomologyUnit unit: units) {
			homologyUnitsMap.put(unit.getId(), unit);
		}
		buildFMIndex();
	}
	private void buildFMIndex() {
		indexHomologyUnits = new FMIndex();
		QualifiedSequenceList unitSequences = new QualifiedSequenceList();
		for (HomologyUnit ql:homologyUnitsMap.values()) {
			String unitSequence = ql.getUnitSequence();
			String unitId = ql.getId();
			QualifiedSequence qualifiedSequence = new QualifiedSequence(unitId, unitSequence);
			unitSequences.add(qualifiedSequence);
		}
		indexHomologyUnits.loadQualifiedSequences(unitSequences);
	}
	
	/**
	 * Returns an FM-index of the homology units
	 * @return FMIndex to search for homologs
	 */
	public FMIndex getIndexHomologyUnits() {
		return indexHomologyUnits;
	}
	
	public List<HomologyUnit> getHomologyUnits() {
		return new ArrayList<>(homologyUnitsMap.values());
	}
	
	public HomologyUnit getHomologyUnit(String unitId) {
		return homologyUnitsMap.get(unitId);
	}
}
