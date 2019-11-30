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
package ngsep.assembly;

import java.io.Serializable;

/**
 * 
 * @author Jorge Duitama
 * @author Juan Camilo Bojaca
 * @author David Guevara
 *
 */
public class AssemblyVertex implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 5185570963647916821L;
	private CharSequence read;
	private boolean start;
	private int index;

	public AssemblyVertex(CharSequence read, boolean start, int index) {
		this.read = read;
		this.start = start;
		this.index = index;
	}

	/**
	 * @return the read
	 */
	public CharSequence getRead() {
		return read;
	}

	/**
	 * @return the start
	 */
	public boolean isStart() {
		return start;
	}
	public int getIndex () {
		return index;
	}
	
}

