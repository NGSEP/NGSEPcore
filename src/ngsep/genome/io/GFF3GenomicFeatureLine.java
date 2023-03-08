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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.genome.GenomicRegion;

public class GFF3GenomicFeatureLine implements GenomicRegion {
	public static final String FEATURE_TYPE_GENE = "gene";
	public static final String FEATURE_TYPE_PCGENE = "protein_coding_gene";
	public static final String FEATURE_TYPE_TRGENE = "transposable_element_gene";
	public static final String FEATURE_TYPE_PSEUDOGENE = "pseudogene";
	public static final String FEATURE_TYPE_EXON = "exon";
	public static final String FEATURE_TYPE_MRNA = "mRNA";
	public static final String FEATURE_TYPE_5PUTR = "five_prime_UTR";
	public static final String FEATURE_TYPE_3PUTR = "three_prime_UTR";
	public static final String FEATURE_TYPE_CDS = "CDS";
	public static final String FEATURE_TYPE_POLYPEPTIDE = "polypeptide";
	public static final String FEATURE_TYPE_SIMILARITY = "similarity";
	
	public static final String [] supportedFeatureTypes = {FEATURE_TYPE_GENE,FEATURE_TYPE_PCGENE, FEATURE_TYPE_PSEUDOGENE, FEATURE_TYPE_CDS,FEATURE_TYPE_MRNA,FEATURE_TYPE_TRGENE,FEATURE_TYPE_5PUTR,FEATURE_TYPE_3PUTR, FEATURE_TYPE_SIMILARITY};
	//TODO: Use a file resource
	public static final String [] supportedFeatureTypesSOFAIDs = {"SO:0000704","","SO:0000336","SO:0000316","SO:0000234","SO:0000111","SO:0000204","SO:0000205",""};
	//Predefined attributes according to the gff3 specification
	public static final String ATTRIBUTE_ID = "ID";
	public static final String ATTRIBUTE_NAME = "Name";
	public static final String ATTRIBUTE_ALIAS = "Alias";
	public static final String ATTRIBUTE_PARENT = "Parent";
	public static final String ATTRIBUTE_TARGET = "Target";
	public static final String ATTRIBUTE_GAP = "Gap";
	public static final String ATTRIBUTE_DERIVES_FROM = "Derives_from";
	public static final String ATTRIBUTE_PRODUCT = "product";
	public static final String ATTRIBUTE_NOTE = "Note";
	public static final String ATTRIBUTE_DBXREF = "Dbxref";
	public static final String ATTRIBUTE_ONTOLOGY = "Ontology_term";
	public static final String ATTRIBUTE_CIRCULAR = "Is_circular";
	private String sequenceName;
	private String source;
	private String type;
	private int first;
	private int last;
	private boolean negativeStrand=false;
	private byte phase = -1;
	private int lineNumber = 0;
	private Map<String,String> annotations = new HashMap<String, String>();
	
	public GFF3GenomicFeatureLine(String sequenceName, int first, int last, String type) {
		super();
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		this.type = type;
	}
	
	public String getSequenceName() {
		return sequenceName;
	}
	
	public int getFirst() {
		return first;
	}
	
	public int getLast() {
		return last;
	}
	
	public String getType() {
		return type;
	}
	
	public String getSource() {
		return source;
	}
	
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	public Map<String, String> getAnnotations() {
		return annotations;
	}
	public void addAnnotation(String key, String value) {
		annotations.put(key, value);
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return !negativeStrand;
	}
	public String getAnnotation(String tag) {
		return annotations.get(tag);
	}
	public void setSource(String source) {
		this.source = source;
	}
	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}
	public byte getPhase() {
		return phase;
	}
	public void setPhase(byte phase) {
		this.phase = phase;
	}
	public String getId() {
		return annotations.get(ATTRIBUTE_ID);
	}
	public String getName() {
		return annotations.get(ATTRIBUTE_NAME);
	}
	public Set<String> getParentIds() {
		String parentsStr = annotations.get(ATTRIBUTE_PARENT);
		if(parentsStr==null) return new HashSet<>();
		String [] parentIdsArr = parentsStr.split(",");
		Set<String> answer = new HashSet<>();
		for(int i = 0;i<parentIdsArr.length;i++) {
			answer.add(parentIdsArr[i]);
		}
		return answer;
	}
	public boolean hasParents() {
		return annotations.get(ATTRIBUTE_PARENT)!=null;
	}
	public List<String> getOntologyTerms() {
		return getValuesList(ATTRIBUTE_ONTOLOGY);
	}
	public List<String> getDatabaseReferences() {
		return getValuesList(ATTRIBUTE_DBXREF);
	}
	public List<String> getProducts() {
		return getValuesList(ATTRIBUTE_PRODUCT);
	}
	public List<String> getValuesList(String attribute) {
		String parentsStr = annotations.get(attribute);
		if(parentsStr==null) return new ArrayList<>();
		String [] parentIdsArr = parentsStr.split(",");
		return Arrays.asList(parentIdsArr);
	}

	public int getLineNumber() {
		return lineNumber;
	}

	public void setLineNumber(int lineNumber) {
		this.lineNumber = lineNumber;
	}
	
	
}