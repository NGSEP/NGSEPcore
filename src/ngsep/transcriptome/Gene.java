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

/**
 * Information for a Gene mapped to a reference genome
 * @author Jorge Duitama
 */
public class Gene {
	private String id;
	private String name;
	/**
	 * Creates a new Gene with the given information
	 * @param id Id of the gene
	 * @param name Name of the gene
	 */
	public Gene(String id, String name) {
		super();
		this.id = id;
		this.name = name;
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
}
