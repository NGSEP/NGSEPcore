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
package ngsep.transcriptome;

import java.util.List;

import ngsep.genome.GenomicRegion;

/**
 * Information for a Gene mapped to a reference genome
 * @author Jorge Duitama
 */
public class Gene implements GenomicRegion {
	private String id;
	private String name;
	private String sequenceName;
	private int first;
	private int last;
	private boolean negativeStrand;
	private List<String> ontologyTerms;
	private List<String> databaseReferences;
	
	/**
	 * Creates a new Gene with the given information
	 * @param id Id of the gene
	 * @param name Name of the gene
	 */
	public Gene(String id, String name, String sequenceName, int first, int last, boolean negativeStrand) {
		this.id = id;
		this.name = name;
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		this.negativeStrand = negativeStrand;
	}
	/**
	 * @return String id of the gene
	 */
	public String getId() {
		return id;
	}
	/**
	 * Changes the id of the gene
	 * @param id New id
	 */
	public void setId(String id) {
		this.id = id;
	}
	/**
	 * @return String Name of the gene
	 */
	public String getName() {
		return name;
	}
	/**
	 * Changes the name of the gene
	 * @param name New gene name
	 */
	public void setName(String name) {
		this.name = name;
	}
	/**
	 * @return the sequenceName
	 */
	public String getSequenceName() {
		return sequenceName;
	}
	/**
	 * @return the first
	 */
	public int getFirst() {
		return first;
	}
	/**
	 * @param first the first to set
	 */
	public void setFirst(int first) {
		this.first = first;
	}
	/**
	 * @return the last
	 */
	public int getLast() {
		return last;
	}
	/**
	 * @param last the last to set
	 */
	public void setLast(int last) {
		this.last = last;
	}
	
	/**
	 * @return the negativeStrand
	 */
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
	 * @return List<String> Ontology terms associated with this gene
	 */
	public List<String> getOntologyTerms() {
		return ontologyTerms;
	}
	/**
	 * Changes the ontology terms associated with this gene
	 * @param ontologyTerms New list
	 */
	public void setOntologyTerms(List<String> ontologyTerms) {
		this.ontologyTerms = ontologyTerms;
	}
	/**
	 * @return List<String> IDs of external databases associated with this gene
	 */
	public List<String> getDatabaseReferences() {
		return databaseReferences;
	}
	/**
	 * Changes the external database ids associated with this gene
	 * @param databaseIds New list
	 */
	public void setDatabaseReferences(List<String> databaseReferences) {
		this.databaseReferences = databaseReferences;
	}
	@Override
	public int length() {
		return last-first+1;
	}
	@Override
	public boolean isPositiveStrand() {
		return !negativeStrand;
	}
	
}
