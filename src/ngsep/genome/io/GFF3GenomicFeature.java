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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class GFF3GenomicFeature {
	private List<GFF3GenomicFeatureLine> lines = new ArrayList<>();
	private List<GFF3GenomicFeature> children = new ArrayList<GFF3GenomicFeature>();
	private String id;
	private String type;
	
	public GFF3GenomicFeature(GFF3GenomicFeatureLine line) {
		id = line.getId();
		type = line.getType();
		lines.add(line);
	}
	
	/**
	 * Adds the given line to the feature
	 * @param line New line to add
	 * @throws IOException If the ID of the new line is different from the ID of this feature
	 * @throws IOException If the type of the new line is different from the type of this feature
	 */
	public void addLine(GFF3GenomicFeatureLine line) throws IOException {
		if(id==null) throw new IOException("Error loading line at "+line.getSequenceName()+":"+line.getFirst()+"-"+line.getLast()+". Features with null ids can not have multiple lines");
		if(!id.equals(line.getId())) throw new IOException("Error loading line at "+line.getSequenceName()+":"+line.getFirst()+"-"+line.getLast()+" the type must be the same for all lines with the same id. Expected: "+type+" given: "+line.getType());;
		if(!type.equals(line.getType())) throw new IOException("Error loading line at "+line.getSequenceName()+":"+line.getFirst()+"-"+line.getLast()+" the type must be the same for all lines with the same id. Expected: "+type+" given: "+line.getType());;
		lines.add(line);
	}
	public List<GFF3GenomicFeatureLine> getLines() {
		return lines;
	}
	public void addChild(GFF3GenomicFeature child) {
		children.add(child);
	}
	public void addChildWithoutId(GFF3GenomicFeatureLine featureLine) {
		children.add(new GFF3GenomicFeature(featureLine)); 
		
	}
	public List<GFF3GenomicFeature> getChildren() {
		return children;
	}
	public String getType() {
		return type;
	}
	public String getId() {
		return id;
	}
	public String getName() {
		for(GFF3GenomicFeatureLine line:lines) {
			String name = line.getName(); 
			if(name!=null) return name; 
		}
		return null;
	}
	public String getAnnotation(String attribute) {
		for(GFF3GenomicFeatureLine line:lines) {
			String ann = line.getAnnotation(attribute); 
			if(ann!=null) return ann; 
		}
		return null;
	}
	public boolean hasParents() {
		for(GFF3GenomicFeatureLine line:lines) {
			if(line.hasParents()) return true; 
		}
		return false;
	}
	public Set<String> getParentIds() {
		Set<String> ids = new HashSet<>();
		for(GFF3GenomicFeatureLine line:lines) {
			ids.addAll(line.getParentIds());
		}
		return ids;
	}
	public List<String> getOntologyTerms() {
		return getValuesList(GFF3GenomicFeatureLine.ATTRIBUTE_ONTOLOGY);
	}
	
	public List<String> getDatabaseReferences() {
		return getValuesList(GFF3GenomicFeatureLine.ATTRIBUTE_DBXREF);
	}
	
	public List<String> getProducts() {
		return getValuesList(GFF3GenomicFeatureLine.ATTRIBUTE_PRODUCT);
	}

	public List<String> getValuesList(String attribute) {
		List<String> values = new ArrayList<>();
		for(GFF3GenomicFeatureLine line:lines) {
			values.addAll(line.getValuesList(attribute));
		}
		return values;
	}
}