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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.main.io.ParseUtils;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

/**
 * Generic loader for files in gff3 format
 * @author Jorge Duitama
 */
public class GFF3Loader {
	private Logger log = Logger.getAnonymousLogger();
	private QualifiedSequenceList sequenceNames;
	
	public GFF3Loader(QualifiedSequenceList sequenceNames) {
		super();
		this.sequenceNames = sequenceNames;
	}
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	public List<GFF3GenomicFeatureLine> loadFeatureLines (String filename) throws IOException {
		InputStream is = null;
		try {
			is = new FileInputStream(filename);
			if(filename.toLowerCase().endsWith(".gz")) {
				is = new ConcatGZIPInputStream(is);
			}
			return loadFeatureLines(is);
		} finally {
			if(is!=null) is.close();
		}
	}
	public List<GFF3GenomicFeatureLine> loadFeatureLines(InputStream is) throws IOException {
		List<GFF3GenomicFeatureLine> gffLines = new ArrayList<>();
		try (BufferedReader in = new BufferedReader(new InputStreamReader(is))) {
			String line=in.readLine();
			for(int i=1;line!=null;i++) {
				if("##FASTA".equals(line.trim())) break;
				if(line.length()>0 && line.charAt(0)!='#') {
					GFF3GenomicFeatureLine featureLine;
					try {
						featureLine = loadFeatureLine(line,i); 
					} catch (RuntimeException e) {
						e.printStackTrace();
						log.warning("Error loading line "+line+". "+e.getMessage());
						line=in.readLine();
						continue;
					}
					
					if(featureLine!=null) gffLines.add(featureLine);
				
				}
				line=in.readLine();
			}
		}
		return gffLines;
	}
	public GFF3GenomicFeatureLine loadFeatureLine (String line, int lineNumber) {
		String [] items = ParseUtils.parseString(line, '\t');
		//Memory saver for sequenceNames
		QualifiedSequence seq;
		try {
			seq = sequenceNames.addOrLookupName(items[0]);
		} catch (Exception e2) {
			log.warning("Can not load genomic feature at line: "+line+". Unrecognized sequence name. "+items[0]);
			return null;
		}
		String type = loadType(items[2]);
		int first;
		try {
			first = Integer.parseInt(items[3]);
		} catch (NumberFormatException e1) {
			log.warning("Error loading feature at line "+line+". Start coordinate must be a positive integer");
			return null;
		}
		int last;
		try {
			last = Integer.parseInt(items[4]);
		} catch (NumberFormatException e1) {
			log.warning("Error loading feature at line "+line+". End coordinate must be a positive integer");
			return null;
		}
		GFF3GenomicFeatureLine answer = new GFF3GenomicFeatureLine(seq.getName(), first, last,type);
		answer.setSource(items[1]);
		answer.setLineNumber(lineNumber);
		answer.setNegativeStrand(items[6].charAt(0)=='-');
		if(items[7].charAt(0)!='.') {
			try {
				answer.setPhase(Byte.parseByte(items[7]));
			} catch (NumberFormatException e) {
				log.warning("Error loading phase of feature at "+answer.getSequenceName()+":"+answer.getFirst()+"-"+answer.getLast()+". Value: "+items[7]);
			}
		} else if (GFF3GenomicFeatureLine.FEATURE_TYPE_CDS.equals(type)) {
			log.warning("Error loading phase of feature at "+answer.getSequenceName()+":"+answer.getFirst()+"-"+answer.getLast()+". The phase field is required for CDS features");
		}
		String errorMessage = "Error loading annotations of feature at "+answer.getSequenceName()+":"+answer.getFirst()+"-"+answer.getLast()+". ";
		String [] annItems = ParseUtils.parseStringWithText(items[8], ';', '\"');
		for(int i=0;i<annItems.length;i++) {
			if(annItems[i].trim().length()==0) continue;
			int idx = annItems[i].indexOf("=");
			if(idx <=0 || idx == annItems[i].length()-1) log.warning(errorMessage+"Annotation "+annItems[i]+" is not a key-value pair");
			else answer.addAnnotation(annItems[i].substring(0,idx),annItems[i].substring(idx+1));
		}
		return answer;
	}
	/**
	 * Returns the known feature type given the type string. This helps to save memory with feature types
	 * @param typeStr Type to search
	 * @return String If typeStr is a supported SOFA id, it returns the correspoonding name.
	 * otherwise it returns a supported type equals to typeStr or typeStr itself if it is not a known type
	 */
	public String loadType(String typeStr) {
		String type = typeStr;
		int idxType = -1;
		//Search type name
		for(int i=0;i<GFF3GenomicFeatureLine.supportedFeatureTypes.length && idxType == -1;i++) {
			if(GFF3GenomicFeatureLine.supportedFeatureTypes[i].equals(typeStr)) {
				idxType = i;
				type = GFF3GenomicFeatureLine.supportedFeatureTypes[idxType];
			}
		}
		//Search type id
		for(int i=0;i<GFF3GenomicFeatureLine.supportedFeatureTypesSOFAIDs.length && idxType == -1;i++) {
			if(GFF3GenomicFeatureLine.supportedFeatureTypesSOFAIDs[i].equals(typeStr)) {
				idxType = i;
				type = GFF3GenomicFeatureLine.supportedFeatureTypes[idxType];
			}
		}
		return type;
	}
}
