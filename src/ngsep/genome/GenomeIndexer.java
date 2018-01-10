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
package ngsep.genome;

import java.io.IOException;

import ngsep.main.CommandsDescriptor;
import ngsep.sequences.FMIndex;

/**
 * Program that build the FM-index related to a genome
 * @author German Andrade
 * @author Jorge Duitama
 */
public class GenomeIndexer 
{

	private int tallyDistance;
	public static void main(String[] args) throws Exception  {
		GenomeIndexer instance = new GenomeIndexer();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String genomeFile=args[i++];
		String outputFile=args[i++];
		instance.createIndex (genomeFile,outputFile);
		
	}

	public void createIndex(String genomeFile, String outputFile) throws IOException {
		System.out.println("Building index for genome in file "+genomeFile);
		long time = System.currentTimeMillis();
		FMIndex fMIndex= new FMIndex();
		fMIndex.loadGenome(genomeFile);
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		System.out.println("Built index in "+seconds+" seconds. Saving in "+outputFile);
		fMIndex.save(outputFile);
		System.out.println("Process completed");
		
	}

	public int getTallyDistance() {
		return tallyDistance;
	}

	public void setTallyDistance(int tallyDistance) {
		this.tallyDistance = tallyDistance;
	}
	public void setTallyDistance(Integer tallyDistance) {
		this.setTallyDistance(tallyDistance.intValue());
	}
	

}
