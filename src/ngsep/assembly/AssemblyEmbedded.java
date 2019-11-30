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
 * @author Jorge Duitama
 */
public class AssemblyEmbedded implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1120528617112863153L;
	
	private int sequenceId;
	private CharSequence read;
	private int startPosition;
	private boolean isReverse;

	public AssemblyEmbedded(int sequenceId, CharSequence read, int startPosition, boolean isReverse) {
		this.sequenceId = sequenceId;
		this.read = read;
		this.startPosition = startPosition;
		this.isReverse = isReverse;
	}

	
	/**
	 * @return the sequenceId
	 */
	public int getSequenceId() {
		return sequenceId;
	}


	public CharSequence getRead() {
		return read;
	}

	public int getStartPosition() {
		return startPosition;
	}

	public void setStartPosition(int startPosition) {
		this.startPosition = startPosition;
	}

	public boolean isReverse() {
		return isReverse;
	}

	public void setReverse(boolean isReverse) {
		this.isReverse = isReverse;
	}
}
