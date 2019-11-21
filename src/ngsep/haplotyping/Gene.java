/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2019 Jorge Duitama
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

import java.util.Random;

public class Gene 
{
	private int value;

	public Gene() 
	{
		generatePosition();
	}

	public void mutatePosition() 
	{
		if(value==0)
			value = 1;
		else
			value = 0;
	}

	public final int getValue()
	{
		return 0;

	}

	private void generatePosition()
	{
		Random r1 = new Random();
		int randomValue = r1.nextInt(Integer.MAX_VALUE);
		value = randomValue%2;		
	}
}
