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
package ngsep.sequences;

import java.io.Serializable;

/**
 * Class to store sequences identifiable by a name and with additional metedata
 * The sequence and quality information is optional so  
 * @author Jorge Duitama
 *
 */
public class QualifiedSequence implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String name;
	private int length;
	private String comments;
	private CharSequence characters;
	private String qualityScores;
	
	
	
	public QualifiedSequence(String name) {
		super();
		this.setName(name);
	}
	
	public QualifiedSequence(String name, CharSequence characters) {
		super();
		this.setName(name);
		this.setCharacters(characters);
	}

	public QualifiedSequence(String name, CharSequence characters, String qualityScores) {
		super();
		this.setName(name);
		this.setCharacters(characters);
		this.setQualityScores(qualityScores);
	}

	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public int getLength() {
		return length;
	}
	public void setLength(int length) {
		if(characters!=null) throw new RuntimeException("Length can not be modified if the sequence has characters");
		this.length = length;
	}
	public String getComments() {
		return comments;
	}
	public void setComments(String comments) {
		this.comments = comments;
	}
	public CharSequence getCharacters() {
		return characters;
	}
	public void setCharacters(CharSequence characters) {
		this.characters = characters;
		this.length = characters.length();
	}
	public String getQualityScores() {
		return qualityScores;
	}
	public void setQualityScores(String qualityScores) {
		this.qualityScores = qualityScores;
	}
}
