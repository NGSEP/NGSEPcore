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
package ngsep.genome.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class SimpleGenomicRegionFileHandler {
	
	private Logger log = Logger.getLogger(SimpleGenomicRegionFileHandler.class.getName());
	private QualifiedSequenceList sequenceNames = new QualifiedSequenceList();
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	public QualifiedSequenceList getSequenceNames() {
		return sequenceNames;
	}
	public void setSequenceNames(QualifiedSequenceList sequenceNames) {
		this.sequenceNames = sequenceNames;
	}
	public List<GenomicRegion> loadRegions (String filename) throws IOException  {
		List<GenomicRegion> regions = new ArrayList<GenomicRegion>();
		
		try (FileReader fis = new FileReader(filename);
				BufferedReader in = new BufferedReader(fis)) {
			String line=in.readLine();
			while(line!=null) {
				String [] items = line.split(" |\t");
				QualifiedSequence seq;
				try {
					seq = sequenceNames.addOrLookupName(items[0]);
					regions.add(new GenomicRegionImpl(seq.getName(), Integer.parseInt(items[1]), Integer.parseInt(items[2])));
				} catch (NumberFormatException e) {
					log.warning("Can not load line: "+line+". "+e.getMessage());
				} catch (RuntimeException e) {
					log.warning("Can not load line: "+line+". Unrecognized sequence name. "+e.getMessage());
				}
				
				line=in.readLine();
			}
		}
		
		return regions;
	}
	
	public Map<String, List<GenomicRegion>> loadRegionsAsMap(String filename) throws IOException {
		Map<String, List<GenomicRegion>> regionsMap = new HashMap<>();
		List<GenomicRegion> regionsList = loadRegions(filename);
		for(GenomicRegion region:regionsList) {
			List<GenomicRegion> regionsSeq = regionsMap.get(region.getSequenceName());
			if(regionsSeq==null) {
				regionsSeq = new ArrayList<>();
				regionsMap.put(region.getSequenceName(), regionsSeq);
			}
			regionsSeq.add(region);
		}
		for(List<GenomicRegion> regionsSeq:regionsMap.values()) {
			Collections.sort(regionsSeq, GenomicRegionPositionComparator.getInstance());
		}
		return regionsMap;
	}
 	public void saveRegions (List<GenomicRegion> regions, PrintStream out) {
		for(GenomicRegion r:regions) {
			out.print(r.getSequenceName()+"\t");
			out.print(r.getFirst()+"\t");
			out.print(r.getLast());
			out.println();
		}
		out.flush();
	}
}
