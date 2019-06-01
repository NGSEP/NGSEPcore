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
package ngsep.transcriptome.io;

import ngsep.genome.ReferenceGenome;
import ngsep.transcriptome.Transcriptome;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class GTF2GFF3Converter {

	public static void main(String[] args) throws Exception {
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		System.err.println("Loaded genome");
		GTF2TranscriptomeHandler gtf2Handler = new GTF2TranscriptomeHandler(genome);
		Transcriptome t = gtf2Handler.loadMap(args[1]);
		System.err.println("Loaded transcriptome with "+t.getAllTranscripts().size()+" transcripts");
		GFF3TranscriptomeWriter writer = new GFF3TranscriptomeWriter();
		writer.printTranscriptome(t, System.out);

	}
}
