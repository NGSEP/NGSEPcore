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
package ngsep.vcf;

import java.util.LinkedHashMap;
import java.util.Map;

public class VCFHeaderLine {
	public static final String ATTRIBUTE_ID= "ID";
	public static final String ATTRIBUTE_DATATYPE= "Type";
	public static final String ATTRIBUTE_NUMBER= "Number";
	public static final String ATTRIBUTE_DESCRIPTION= "Description";
	
	private String headerType;
	private String id;
	private String description;
	private String number;
	private String dataType;
	private Map<String,String> attributes = new LinkedHashMap<String,String>();
	public VCFHeaderLine(String headerType, String id) {
		super();
		this.headerType = headerType;
		this.id = id;
	}
	
	public VCFHeaderLine(String headerType, String id, String description, String number, String dataType) {
		super();
		this.headerType = headerType;
		this.setId(id);
		this.setNumber(number);
		this.setDataType(dataType);
		this.setDescription(description);
	}

	public String getHeaderType() {
		return headerType;
	}
	public void setHeaderType(String headerType) {
		this.headerType = headerType;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
		attributes.put(ATTRIBUTE_ID, id);
	}
	public String getDescription() {
		return description;
	}
	public void setDescription(String description) {
		this.description = description;
		attributes.put(ATTRIBUTE_DESCRIPTION, description);
	}
	public String getNumber() {
		return number;
	}
	public void setNumber(String number) {
		this.number = number;
		attributes.put(ATTRIBUTE_NUMBER, number);
	}
	public String getDataType() {
		return dataType;
	}
	public void setDataType(String dataType) {
		this.dataType = dataType;
		attributes.put(ATTRIBUTE_DATATYPE, dataType);
	}
	public String getAttribute(String name) {
		return attributes.get(name);
	}
	public void setAttribute(String property, String value) {
		attributes.put(property, value);
	}
	public Map<String, String> getAttributes() {
		return attributes;
	}
	
	
}
