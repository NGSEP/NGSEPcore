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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.sequences.AminoacidSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transcriptome.ProteinTranslator;

/**
 * 
 * @author Jorge Gomez
 *
 */
public class HomologyCatalog {
	public static final int INPUT_TYPE_CDNA = 1;
	public static final int INPUT_TYPE_PROTEIN = 2;
	private Map<String, HomologyUnit> homologyUnitsMap= new HashMap<String, HomologyUnit>();
	private FMIndex indexHomologyUnits=null;
	
	public HomologyCatalog (List<HomologyUnit> units) {
		for(HomologyUnit unit: units) {
			homologyUnitsMap.put(unit.getUniqueKey(), unit);
		}
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
	
	public static HomologyCatalog loadFromFasta(String filename, int genomeId, int inputType) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		if(inputType==INPUT_TYPE_PROTEIN) handler.setSequenceType(AminoacidSequence.class);
		List<QualifiedSequence> sequences = handler.loadSequences(filename);
		ProteinTranslator translator = new ProteinTranslator();
		List<HomologyUnit> units = new ArrayList<>();
		for(QualifiedSequence seq : sequences) {
			CharSequence aminoacidSequence = seq.getCharacters();
			CharSequence cdsSequence = seq.getCharacters();
			if(inputType==INPUT_TYPE_CDNA) aminoacidSequence = new AminoacidSequence(translator.getProteinSequence(aminoacidSequence));
			if(aminoacidSequence.length()<10) System.out.println("Small sequence "+aminoacidSequence+" with name: "+seq.getName()+" length: "+aminoacidSequence.length());
			HomologyUnit unit = new HomologyUnit(genomeId, seq.getName(), aminoacidSequence);
			if(inputType==INPUT_TYPE_CDNA) unit.setCdsSequence(cdsSequence);
			units.add(unit);
		}
		
		HomologyCatalog catalog = new HomologyCatalog(units);
		return catalog;
	}
}
