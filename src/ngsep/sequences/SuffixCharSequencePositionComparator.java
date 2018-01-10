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
import java.util.Comparator;

/**
 * Comparator for integers interpreted as starts of suffixes of a single sequence
 * @author German Andrade
 * @author Jorge Duitama
 */
public class SuffixCharSequencePositionComparator implements Comparator<Integer> 
{
	private CharSequence word;
	private int length;
	public SuffixCharSequencePositionComparator(CharSequence word) 
	{
		super();
		this.word = word;
		this.length = word.length();
	}
	
	@Override
	public int compare(Integer o1, Integer o2) {
		while (o1<length && o2<length) {
			char c1 = word.charAt(o1);
			char c2 = word.charAt(o2); 
			if(c1!=c2) return c1-c2;
			o1++;
			o2++;
		}
		if( o2<length) return -1;
		if( o1<length) return 1;
		return 0;
	}

}
