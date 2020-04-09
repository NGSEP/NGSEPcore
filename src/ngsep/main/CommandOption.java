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
package ngsep.main;

import java.lang.reflect.Method;

public class CommandOption {
	public static final String TYPE_INT = "INT";
	public static final String TYPE_LONG = "LONG";
	public static final String TYPE_FLOAT = "FLOAT";
	public static final String TYPE_DOUBLE = "DOUBLE";
	public static final String TYPE_GENOME = "GENOME";
	public static final String TYPE_STRING = "STRING";
	public static final String TYPE_FILE = "FILE";
	public static final String TYPE_DIR = "DIR";
	public static final String TYPE_BOOLEAN = "BOOLEAN";
	
	private String id;
	private String type;
	private String defaultValue=null;
	private String description;
	private String attribute;
	private boolean deprecated = false;
	
	public CommandOption(String id) {
		super();
		this.id = id;
	}
	public CommandOption(String id, String type, String defaultValue, String description) {
		super();
		this.id = id;
		this.type = type;
		this.defaultValue = defaultValue;
		this.description = description;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getType() {
		return type;
	}
	public void setType(String type) {
		this.type = type;
	}
	public String getDefaultValue() {
		return defaultValue;
	}
	public void setDefaultValue(String defaultValue) {
		this.defaultValue = defaultValue;
	}
	public String getDescription() {
		return description;
	}
	public void setDescription(String description) {
		this.description = description;
	}
	
	public String getAttribute() {
		return attribute;
	}
	public void setAttribute(String attribute) {
		this.attribute = attribute;
	}
	
	public boolean isDeprecated() {
		return deprecated;
	}
	public void setDeprecated(boolean deprecated) {
		this.deprecated = deprecated;
	}
	public boolean printType () {
		return type!=null && !TYPE_BOOLEAN.equals(type); 
	}
	public int getPrintLength() {
		int length = id.length()+9;
		if(printType()) length+=type.length()+1;
		return length;
	}
	public Method findGetMethod (Object instance) {
		if(attribute==null) throw new RuntimeException("Attribute not set for option: "+id);
		String methodName = "get"+Character.toUpperCase(attribute.charAt(0));
		if(attribute.length()>0) methodName+=attribute.substring(1);
		try {
			return instance.getClass().getMethod(methodName);
		} catch (NoSuchMethodException | SecurityException e) {
			throw new RuntimeException(e);
		}
	}
	public Method findSetMethod (Object instance) {
		if(attribute==null) throw new RuntimeException("Attribute not set for option: "+id);
		String methodName = "set"+Character.toUpperCase(attribute.charAt(0));
		if(attribute.length()>0) methodName+=attribute.substring(1);
		try {
			return instance.getClass().getMethod(methodName,getTypeClass());
		} catch (NoSuchMethodException | SecurityException e) {
			try {
				return instance.getClass().getMethod(methodName,String.class);
			} catch (NoSuchMethodException | SecurityException e1) {
				throw new RuntimeException(e);
			}
		}
	}
	public Method findStringSetMethod (Object instance) {
		if(attribute==null) throw new RuntimeException("Attribute not set for option: "+id);
		String methodName = "set"+Character.toUpperCase(attribute.charAt(0));
		if(attribute.length()>0) methodName+=attribute.substring(1);
		try {
			return instance.getClass().getMethod(methodName,String.class);
		} catch (NoSuchMethodException | SecurityException e) {
			throw new RuntimeException(e);
		}
	}
	private Class<?> getTypeClass() {
		if(TYPE_BOOLEAN.equals(type)) {
			return Boolean.class;
		}
		if(TYPE_INT.equals(type)) {
			return Integer.class;
		}
		if(TYPE_LONG.equals(type)) {
			return Long.class;
		}
		if(TYPE_FLOAT.equals(type)) {
			return Float.class;
		}
		if(TYPE_DOUBLE.equals(type)) {
			return Double.class;
		}
		if(TYPE_GENOME.equals(type)) {
			return String.class;
		}
		if(TYPE_STRING.equals(type)) {
			return String.class;
		}
		if(TYPE_FILE.equals(type)) {
			return String.class;
		}
		if(TYPE_DIR.equals(type)) {
			return String.class;
		}
		throw new RuntimeException("Can not decode option of unrecognized type: "+type);
	}
	public Object decodeValue (String value) {
		return OptionValuesDecoder.decode(value, getTypeClass());
	}
}
