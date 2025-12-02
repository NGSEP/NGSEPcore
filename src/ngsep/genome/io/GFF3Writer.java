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

import java.io.PrintStream;
import java.util.List;
import java.util.Map;

/**
 * Generic writer for files in gff3 format
 * @author Jorge Duitama
 */
public class GFF3Writer {
	public void printHeader(PrintStream out) {
		out.println("##gff-version 3");
	}
	public void printGffFeatureLines(List<GFF3GenomicFeatureLine> lines, PrintStream out) {
		for(GFF3GenomicFeatureLine line: lines) {
			printGffFeatureLine(line, out);
		}
	}
	public void printGffFeatureLine(GFF3GenomicFeatureLine line, PrintStream out) {
		byte phase = line.getPhase();
		boolean positiveStrand = line.isPositiveStrand();
		String id = line.getId();
		out.print(line.getSequenceName());
		out.print("\t"+line.getSource());
		out.print("\t"+line.getType());
		out.print("\t"+line.getFirst());
		out.print("\t"+line.getLast());
		out.print("\t"+((line.getQuality()>=0)?""+line.getQuality():"."));
		out.print("\t"+(positiveStrand?'+':'-'));
		out.print("\t"+((phase>=0)?""+phase:"."));
		boolean attributesStarted = false;
		if(id!=null) {
			out.print("\t"+GFF3GenomicFeatureLine.ATTRIBUTE_ID+"="+id);
			attributesStarted = true;
		}
		Map<String,String> anns = line.getAnnotations();
		for(Map.Entry<String, String> entry:anns.entrySet()) {
			String attId = entry.getKey();
			if(!GFF3GenomicFeatureLine.ATTRIBUTE_ID.equals(attId)) {
				out.print(attributesStarted?';':'\t');
				out.print(attId+"="+entry.getValue());
				attributesStarted = true;
			}
		}
		out.println();
	}
	
	public String getOntologyId(String typeId) {
		for(int i=0;i<GFF3GenomicFeatureLine.supportedFeatureTypes.length;i++) {
			if(GFF3GenomicFeatureLine.supportedFeatureTypes[i].equals(typeId)) {
				return GFF3GenomicFeatureLine.supportedFeatureTypesSOFAIDs[i];
			}
		}
		return null;
	}
}
