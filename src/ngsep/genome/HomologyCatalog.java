package ngsep.genome;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.ProteinTranslator;

public class HomologyCatalog {
	public static final int INPUT_TYPE_CDNA = 1;
	public static final int INPUT_TYPE_PROTEIN = 2;
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
		indexHomologyUnits.loadQualifiedSequences(unitSequences, null);
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
	public static HomologyCatalog loadFromFasta(String filename, int genomeId, int inputType) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		if(inputType==INPUT_TYPE_PROTEIN) handler.setSequenceType(StringBuilder.class);
		List<QualifiedSequence> sequences = handler.loadSequences(filename);
		ProteinTranslator translator = new ProteinTranslator();
		List<HomologyUnit> units = new ArrayList<>();
		for(QualifiedSequence seq : sequences) {
			String proteinSequence = seq.getCharacters().toString();
			if(inputType==INPUT_TYPE_CDNA) proteinSequence = translator.getProteinSequence(proteinSequence);
			if(proteinSequence.length()<10) System.out.println("Small sequence "+proteinSequence+" with name: "+seq.getName()+" length: "+proteinSequence.length());
			HomologyUnit unit = new HomologyUnit(genomeId, seq.getName(), proteinSequence);
			units.add(unit);
		}
		
		HomologyCatalog catalog = new HomologyCatalog(units);
		return catalog;
	}
}
