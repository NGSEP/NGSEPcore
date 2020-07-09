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
package ngsep.alignments;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNASequence;


/**
 * Program that takes a set of alignments and a reference genome and counts the
 * number of mismatches with the reference for each read position from 5' to 3'
 * end
 * @author Jorge Duitama
 */
public class BasePairQualityStatisticsCalculator {
	
	// Constants for default values
	public static final int DEF_MIN_MQ_UNIQUE_ALIGNMENT = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	// Logging and progress
	private Logger log = Logger.getLogger(BasePairQualityStatisticsCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private List<String> inputFiles = new ArrayList<String>();
	private ReferenceGenome genome = null;
	private String outputFile;
	private int minMQ = DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	// Model attributes
	private List<Long> mismatches=new ArrayList<Long>();
	private List<Long> alnsByLength=new ArrayList<Long>();
	private long totalAlignments = 0;
	private long totalBases = 0;
	private List<Long> mismatchesReadsUniqueMapping=new ArrayList<Long>();
	private List<Long> readsUniqueMappingByLength=new ArrayList<Long>();
	private long readsUniqueMapping = 0;
	private long basesUniqueMapping = 0;
	
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

	public List<String> getInputFiles() {
		return inputFiles;
	}
	public void setInputFiles(List<String> inputFiles) {
		this.inputFiles = inputFiles;
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	public void setGenome(String genomeFile) throws IOException {
		setGenome(OptionValuesDecoder.loadGenome(genomeFile,log));
	}
	
	public int getMinMQ() {
		return minMQ;
	}
	public void setMinMQ(int minMQ) {
		this.minMQ = minMQ;
	}
	public void setMinMQ(String value) {
		this.setMinMQ((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public long getTotalAlignments() {
		return totalAlignments;
	}

	public long getTotalBases() {
		return totalBases;
	}

	public long getReadsUniqueMapping() {
		return readsUniqueMapping;
	}

	public long getBasesUniqueMapping() {
		return basesUniqueMapping;
	}
	
	public static void main(String[] args) throws Exception {
		BasePairQualityStatisticsCalculator instance = new BasePairQualityStatisticsCalculator();
		int i=CommandsDescriptor.getInstance().loadOptions(instance, args);
		for(;i<args.length;i++) {
			instance.inputFiles.add(args[i]);
		}
		instance.run();
	}
	
	public void run() throws IOException {
		if (genome == null) throw new IOException("The file with the reference genome is a required parameter");
		init();
		for(String alnsFile:inputFiles) {
			log.info("Processing alignments file "+alnsFile);
			processFile(alnsFile);
		}
		if(outputFile==null) printStatistics(System.out);
		else {
			try (PrintStream out = new PrintStream(outputFile)) {
				printStatistics(out);
			}
		}
		log.info("Process finished");
	}
	
	public void init() {
		mismatches.clear();
		mismatchesReadsUniqueMapping.clear();
		alnsByLength.clear();
		readsUniqueMappingByLength.clear();
		totalAlignments=0;
		totalBases=0;
		readsUniqueMapping=0;
		basesUniqueMapping=0;
	}
	/**
	 * Updates the counts of reads, alignments and mismatches processing the alignments within the given file
	 * @param filename SAM/BAM file with alignments
	 * @throws IOException If the file can not be read
	 */
	public void processFile(String filename) throws IOException {
		try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(filename)) {
			reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_SEQUENCE);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			reader.setFilterFlags(filterFlags);
			reader.setMinMQ(minMQ);
			Iterator<ReadAlignment> it = reader.iterator();
			while (it.hasNext()) {
				ReadAlignment aln = it.next();
				CharSequence read = aln.getReadCharacters();
				boolean reverse = aln.isNegativeStrand();
				
				boolean uniqueAln = aln.isUnique();
				int basesInRegion = 0;
				int readLength = read.length();
				while(readLength>mismatches.size()) {
					long z = 0;
					mismatches.add(z);
					mismatchesReadsUniqueMapping.add(z);
					alnsByLength.add(z);
					readsUniqueMappingByLength.add(z);
				}
				for (int j = 0; j < read.length(); j++) {
					char baseRead = Character.toUpperCase(read.charAt(j));
					if (DNASequence.isInAlphabeth(baseRead)) {
						int pos = aln.getReferencePositionAlignedRead(j);
						if (pos > 0) {
							char baseRef = Character.toUpperCase(genome.getReferenceBase(aln.getSequenceName(),pos));
							if (baseRef != 0 && DNASequence.isInAlphabeth(baseRef)) {
								basesInRegion++;
								if (baseRead != baseRef) {
									int posStats = j;
									if (reverse) posStats = read.length() - 1 - j;
									addOne(mismatches,posStats);
									
									if (uniqueAln) addOne(mismatchesReadsUniqueMapping,posStats);
								}
							}
						}
					}
				}
				if (basesInRegion > 0) {
					
					totalAlignments++;
					totalBases += basesInRegion;
					addOne(alnsByLength,readLength-1);
					if (uniqueAln) {
						readsUniqueMapping++;
						basesUniqueMapping += basesInRegion;
						addOne(readsUniqueMappingByLength,readLength-1);
					}
				}
				
				if (progressNotifier!=null && totalAlignments%10000==0) {
					int progrees = (int) (totalAlignments/10000);
					if (!progressNotifier.keepRunning(progrees)) break;
				}
			}
		}
	}
	
	private void addOne(List<Long> list, int i) {
		long value = list.get(i);
		list.set(i, value+1);
	}

	/**
	 * Prints counts to the given stream
	 * @param out Stream to print the statistics
	 */
	public void printStatistics(PrintStream out) {
		int n = mismatches.size();
		long  [] cumulativeAlnsByLength = BasePairQualityStatisticsCalculator.makeCumulative(alnsByLength);
		long  [] cumulativeRUM = BasePairQualityStatisticsCalculator.makeCumulative(readsUniqueMappingByLength);
		for (int i = 0; i < n; i++) {
			out.println("" + (i + 1) + "\t" + mismatches.get(i) + "\t"+ mismatchesReadsUniqueMapping.get(i)+ "\t" + cumulativeAlnsByLength[i] + "\t"+ cumulativeRUM[i]);
		}
		out.println();
		out.println("Alignments\t" + totalAlignments + "\t"+ readsUniqueMapping);
		out.println("Bases\t" + totalBases + "\t"+ basesUniqueMapping);
	}
	
	/**
	 * Calculates the percentage of reads/alignments having a base different than the reference
	 * @param uniqueAlignments true if only unique alignments need to be taken into account
	 * @return double [] Percentage [0 to 100] of reads/alignments having a base different than the reference for each read position from 5' to 3'
	 * 
	 */
	public double [] calculatePercentages(boolean uniqueAlignments) {
		int n = mismatches.size();
		double [] percentages = new double [n];
		long  [] cumulative;
		if(uniqueAlignments) cumulative = BasePairQualityStatisticsCalculator.makeCumulative(readsUniqueMappingByLength); 
		else cumulative = BasePairQualityStatisticsCalculator.makeCumulative(alnsByLength); 
		for (int i = 0; i < n; i++) {
			if(uniqueAlignments) percentages[i]= 100.0*mismatchesReadsUniqueMapping.get(i);
			else percentages[i]= 100.0*mismatches.get(i);
			percentages[i]/=(double)cumulative[i];
		}
		return percentages;
	}
	
	/**
	 * Calculates the percentages of mismatches reading the counts from a statistics file
	 * @param statsFile File from which statistics are calculated
	 * @param uniqueAlignments true if only unique alignments need to be taken into account
	 * @return double [] Percentage [0 to 100] of reads/alignments having a base different than the reference for each read position from 5' to 3'
	 * @throws IOException If the file can not be read 
	 */
	public static double [] calculatePercentages(String statsFile, boolean uniqueAlignments) throws IOException {
		List<Double> percentages = new ArrayList<Double>();
		FileInputStream fis = null;
		BufferedReader in = null;
		try {
			fis = new FileInputStream(statsFile);
			in = new BufferedReader(new InputStreamReader(fis));
			String line = in.readLine();
			while (line != null) {
				String[] items = line.split("\t");
				if(items.length>=5) {
					double mismatches;
					double count;
					if(uniqueAlignments) {
						mismatches = Integer.parseInt(items[2]);
						count = Integer.parseInt(items[4]);
					} else {
						mismatches = Integer.parseInt(items[1]);
						count = Integer.parseInt(items[3]);
					}
					percentages.add(100.0*mismatches/count);
				}
				line = in.readLine();
			}
		} finally {
			if (in != null)
				in.close();
			if (fis != null)
				fis.close();
		}
		double [] answer = new double [percentages.size()];
		for(int i=0;i<answer.length;i++) answer[i] = percentages.get(i);
		return answer;
	}

	private static long[] makeCumulative(List<Long> data) {
		int n = data.size();
		long [] answer = new long[n];
		Arrays.fill(answer, 0);
		for(int i=0;i<n;i++) {
			for(int j=i;j<n;j++) {
				answer[i]+=data.get(j);
			}
		}
		return answer;
	}

}
