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
package ngsep.discovery;

import java.io.IOException;
import java.io.PrintStream;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;


public class CoverageStatisticsCalculator implements PileupListener {
	
	private Logger log = Logger.getLogger(CoverageStatisticsCalculator.class.getName());
	//Pileup generator
	private AlignmentsPileupGenerator generator;
	//Progress tracking for external control
	private ProgressNotifier progressNotifier = null;
	private long coveredGenomeSize = 0;
	private long genomeSizeBAMFile = 0;
	private String outFilename = null;
	private int minMQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		CoverageStatisticsCalculator calculator = new CoverageStatisticsCalculator();
		int i = CommandsDescriptor.getInstance().loadOptions(calculator, args);
		String inFile = args[i++];
		calculator.outFilename = args[i++];
		calculator.processFile(inFile);
	}
	
	/**
	 * @return the minMQ
	 */
	public int getMinMQ() {
		return minMQ;
	}

	/**
	 * @param minMQ the minMQ to set
	 */
	public void setMinMQ(int minMQ) {
		this.minMQ = minMQ;
	}
	
	/**
	 * @param minMQ the minMQ to set
	 */
	public void setMinMQ(Integer minMQ) {
		this.setMinMQ(minMQ.intValue());
	}
	
	public void processFile(String filename) throws IOException {
		CoverageStatsPileupListener listener = new CoverageStatsPileupListener();
		coveredGenomeSize = 0;
		loadGenomeSize(filename);
		if(genomeSizeBAMFile == 0) genomeSizeBAMFile = 1000000000;
		generator = new AlignmentsPileupGenerator();
		generator.setLog(log);
		generator.setProcessSecondaryAlignments(true);
		generator.setMaxAlnsPerStartPos(100);
		generator.setMinMQ(minMQ);
		generator.addListener(listener);
		generator.addListener(this);
		generator.processFile(filename);
		if(outFilename!=null) {
			try (PrintStream outFile = new PrintStream(outFilename)){
				listener.printCoverageStats(outFile);
			}
		}
		
		 
	}
	@Override
	public void onPileup(PileupRecord pileup) {
		coveredGenomeSize++;
		if(progressNotifier!=null && coveredGenomeSize%10000==0) {
			int progress = 10+(int)Math.round(85.0*coveredGenomeSize/genomeSizeBAMFile);
			generator.setKeepRunning(progressNotifier.keepRunning(progress));
		}
		
	}
	@Override
	public void onSequenceStart(QualifiedSequence sequenceName) {
		
	}
	@Override
	public void onSequenceEnd(QualifiedSequence sequenceName) {
		
	}
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
	
	
	/**
	 * @return the log
	 */
	public Logger getLog() {
		return log;
	}

	/**
	 * @param log the log to set
	 */
	public void setLog(Logger log) {
		this.log = log;
	}

	/**
	 * @return the outFilename
	 */
	public String getOutFilename() {
		return outFilename;
	}

	/**
	 * @param outFilename the outFilename to set
	 */
	public void setOutFilename(String outFilename) {
		this.outFilename = outFilename;
	}

	private void loadGenomeSize(String alignmentsFile) throws IOException {
		//Load from the alignments file
		ReadAlignmentFileReader reader = null;
		try {
			reader = new ReadAlignmentFileReader(alignmentsFile);
			genomeSizeBAMFile = reader.getSequences().getTotalLength();
		} finally {
			if(reader!=null)reader.close();
		}
	}
	

}
