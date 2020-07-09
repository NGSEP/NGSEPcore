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
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;

/**
 * Program that build the FM-index related to a genome
 * @author German Andrade
 * @author Jorge Duitama
 */
public class GenomeIndexer {
	// Constants for default values
	
	
	// Logging and progress
	private Logger log = Logger.getLogger(GenomeIndexer.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String outputFile = null;
	
	
	// Get and set methods
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}

	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}

	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}

	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}

	public static void main(String[] args) throws Exception  {
		GenomeIndexer instance = new GenomeIndexer();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run () throws IOException {
		if (inputFile==null) throw new IOException("The reference genome is a required parameter");
		if (outputFile==null) throw new IOException("The path of the output file is a required parameter");
		createIndex (inputFile,outputFile);
	}

	public void createIndex(String genomeFile, String outputFile) throws IOException {
		log.info("Loading genome from file "+genomeFile);
		ReferenceGenome genome = new ReferenceGenome(genomeFile);
		log.info("Building index for genome in file "+genomeFile);
		long time = System.currentTimeMillis();
		ReferenceGenomeFMIndex fMIndex= new ReferenceGenomeFMIndex(genome);
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		log.info("Built index in "+seconds+" seconds. Saving in "+outputFile);
		fMIndex.save(outputFile);
		log.info("Process completed");
	}
}
