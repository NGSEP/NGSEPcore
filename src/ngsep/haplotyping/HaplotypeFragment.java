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
package ngsep.haplotyping;

import ngsep.variants.CalledGenomicVariant;

public class HaplotypeFragment {
	
	/**
	 * Represents the first column of a fragment.
	 */
	private int firstColumn;
	
	/**
	 * Represents the list of calls of a fragment.
	 */
	private byte calls [];
	
	/**
	 * Constructor that initializes the attributes of a HaplotypeFragment with the given parameters.
	 * @param firstColumn.
	 * @param calls.
	 */
	public HaplotypeFragment(int firstColumn, byte[] calls) {
		
		this.firstColumn = firstColumn;
		this.calls = calls;
	}

	/**
	 * Returns the first column of a fragment.
	 * @return first column.
	 */
	public int getFirstColumn() {
		return firstColumn;
	}
	
	/**
	 * Returns the list of calls of fragment.
	 * @return calls.
	 */
	public byte[] getCalls() {
		return calls;
	}
	
	/**
	 * Returns the call in a given column of the fragment.
	 * @param column.
	 * @return call.
	 */
	public byte getCall(int column)
	{
		byte call = CalledGenomicVariant.ALLELE_UNDECIDED;
		int posIni = getFirstColumn();
		int length = calls.length;
		int posLas = posIni + length - 1;
		if (column >= posIni && column <= posLas)
		{
			int relativePosition = column - posIni;
			call = calls[relativePosition];
		}
		return call;	
		
	}
	/**
	 * Returns the last column of a fragment.
	 * @return last column 
	 */
	public int getLastColumn()
	{
		int lastColumn = getFirstColumn() + calls.length -1;
		return lastColumn;
	}
	
	

}
