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

/**
 * @author Nicolas Rozo
 */
public class MinimapShortKmerCodesHashFunction implements ShortKmerCodesHashFunction {

	@Override
	public int getHash(long dnaHash) {
		long range = (long)(Math.pow(2, Long.toBinaryString(dnaHash).length())-1);
		long answer = (~dnaHash +(dnaHash << 21)) & range;
		answer = (answer ^ answer >> 24);
		answer = (answer + (answer << 3) + (answer << 8)) & range;
		answer = (answer ^ answer >> 14);
		answer = (answer + (answer << 2) + (answer << 4)) & range;
		answer = (answer ^ answer >> 28);
		answer = (answer + (answer << 31)) & range;

		return (int) answer;
	}

}
