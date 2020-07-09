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
import java.util.Arrays;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;


public class CoverageStatisticsCalculator implements PileupListener {
	
	// Constants for default values
	public static final int DEF_MIN_MQ_UNIQUE_ALIGNMENT = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
		
	// Logging and progress
	private Logger log = Logger.getLogger(CoverageStatisticsCalculator.class.getName());
	private ProgressNotifier progressNotifier = null;
	private long coveredGenomeSize = 0;
	private long genomeSizeBAMFile = 0;
	
	// Parameters
	private String inputFile = null;
	private String outputFile = null;
	private int minMQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	// Model attributes
	private AlignmentsPileupGenerator generator;
	private int maxCoverage = 300;
	private int [] coverageCounts;
	private int highCoverageCount = 0;
	private int [] coverageCountUniqueAlignments;
	private int highCoverageCountUniqueAlignments = 0;
	
	
	
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
	
	public int getMinMQ() {
		return minMQ;
	}
	public void setMinMQ(int minMQ) {
		this.minMQ = minMQ;
	}
	public void setMinMQ(Integer minMQ) {
		this.setMinMQ(minMQ.intValue());
	}
	
	public static void main(String[] args) throws Exception {
		CoverageStatisticsCalculator instance = new CoverageStatisticsCalculator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
		
	}
	
	public void run() throws IOException {
		if (inputFile == null) throw new IOException("The alignments input file is a required parameter");
		processFile(inputFile,outputFile);
		log.info("Process finished");
	}
	public void processFile(String inputFile, String outputFile) throws IOException {
		reset();
		coveredGenomeSize = 0;
		try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(inputFile)){ 
			genomeSizeBAMFile = reader.getSequences().getTotalLength();
		}
		if(genomeSizeBAMFile == 0) genomeSizeBAMFile = 1000000000;
		generator = new AlignmentsPileupGenerator();
		generator.setLog(log);
		generator.setProcessSecondaryAlignments(true);
		generator.setMaxAlnsPerStartPos(100);
		generator.setMinMQ(minMQ);
		generator.addListener(this);
		generator.processFile(inputFile);
		if(outputFile!=null) {
			try (PrintStream out = new PrintStream(outputFile)){
				printCoverageStats(out);
			}
		} else printCoverageStats(System.out);	 
	}
	@Override
	public void onPileup(PileupRecord pileup) {
		processPileup(pileup);
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
	public void reset() {
		this.coverageCounts = new int[maxCoverage];
		Arrays.fill(coverageCounts, 0);
		highCoverageCount = 0;
		this.coverageCountUniqueAlignments = new int[maxCoverage];
		Arrays.fill(coverageCountUniqueAlignments, 0);
		highCoverageCountUniqueAlignments = 0;
	}
	public int[] getCoverageCounts() {
		return coverageCounts;
	}
	public int getMaxCoverage() {
		return maxCoverage;
	}
	public void setMaxCoverage(int maxCoverage) {
		this.maxCoverage = maxCoverage;
		reset();
	}


	public int getHighCoverageCount() {
		return highCoverageCount;
	}

	public void processPileup (PileupRecord record) {
		int calls = record.getNumAlignments();
		if(calls < coverageCounts.length) {
			coverageCounts[calls]++;
		} else {
			highCoverageCount++;
		}
		int callsUniqueAlns = record.getNumUniqueAlns();
		if(callsUniqueAlns < coverageCountUniqueAlignments.length) {
			coverageCountUniqueAlignments[callsUniqueAlns]++;
		} else {
			highCoverageCountUniqueAlignments++;
		}
	}
	
	public int getCoverageMaxCount() {
		int maxIndex = 1;
		for(int i=1;i<coverageCounts.length;i++) {
			if(coverageCounts[maxIndex]<coverageCounts[i]) {
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	public void printCoverageStats(PrintStream out) {
		for(int i=1;i<coverageCounts.length;i++) {
			out.println(""+i+"\t"+coverageCounts[i]+"\t"+coverageCountUniqueAlignments[i]);
		}
		out.println("More\t"+highCoverageCount+"\t"+highCoverageCountUniqueAlignments);
		out.flush();
	}
}
