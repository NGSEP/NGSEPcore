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
package ngsep.discovery.rd;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.ReferenceGenome;
import ngsep.math.Distribution;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class ReadDepthDistribution {
	public static final int DEFAULT_BIN_SIZE=100;
	private Logger log = Logger.getLogger(ReadDepthDistribution.class.getName());
	//Parameters set before starting
	private int binSize = DEFAULT_BIN_SIZE;
	private Map<String, List<ReadDepthBin>> bins = new TreeMap<String, List<ReadDepthBin>>();
	private QualifiedSequenceList sequences;
	private long genomeSize = 0;
	private long totalReads = 0;
	//Parameters of the distribution of read depth
	private double meanReadDepth=0;
	private double sigmaReadDepth = 1;
	private int minMQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	public long getGenomeSize() {
		return genomeSize;
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

	public QualifiedSequenceList getSequences() {
		return sequences;
	}

	public double getMeanReadDepth() {
		return meanReadDepth;
	}

	public double getSigmaReadDepth() {
		return sigmaReadDepth;
	}
	
	public int getBinSize() {
		return binSize;
	}
	
	public long getTotalReads() {
		return totalReads;
	}

	public ReadDepthDistribution(String partitionFile, boolean includeLevels) throws IOException {
		long partitionGenomeSize = 0;
		sequences = new QualifiedSequenceList();
		FileInputStream fis = new FileInputStream(partitionFile);
		BufferedReader in = new BufferedReader(new InputStreamReader(fis));
		String line = in.readLine();
		while (line != null) {
			String[] items = line.split("\t| ");
			String seqName = sequences.addOrLookupName(items[0]).getName(); 
			List<ReadDepthBin> binsSeq = bins.get(seqName);
			if(binsSeq==null) {
				binsSeq = new ArrayList<ReadDepthBin>();
				bins.put(seqName, binsSeq);
			}
			int first = Integer.parseInt(items[1]);
			int last = Integer.parseInt(items[2]);
			ReadDepthBin bin = new ReadDepthBin(seqName, first, last, Double.parseDouble(items[3])/100.0);
			bin.setRawReadDepth(Double.parseDouble(items[4]));
			bin.setCorrectedReadDepth(Double.parseDouble(items[5]));
			if(includeLevels) bin.setReadDepthLevel(Double.parseDouble(items[6]));
			binsSeq.add(bin);
			if(partitionGenomeSize==0) {
				//Obtain bin size from the first bin
				binSize = last - first + 1;
			}
			partitionGenomeSize+=binSize;
			line = in.readLine();
		}
		fis.close();
		genomeSize = partitionGenomeSize;
	}
	
	public ReadDepthDistribution(ReferenceGenome genome, int binSize) {
		
		this.binSize = binSize;
		genomeSize = genome.getTotalLength();
		sequences = genome.getSequencesMetadata();
		System.out.println("Number of sequences: "+sequences.size());
		int n = genome.getNumSequences();
		for(int h=0;h<n;h++) {
			QualifiedSequence sequence = genome.getSequenceByIndex(h);
			String seqName = sequence.getName();
			List<ReadDepthBin> seqBins = new ArrayList<ReadDepthBin>(); 
			bins.put(seqName, seqBins);
			CharSequence sequenceChars = sequence.getCharacters();
			int l = sequenceChars.length();
			int end = l - l%binSize;
			//Ignore the last basepairs to avoid going over the end of the chromosome
			for(int i=0;i<end;i+=binSize) {
				String binSequence = sequenceChars.subSequence(i, i+binSize).toString();
				double gcContent = 0;
				int nBases = 0;
				for(int j=0;j<binSize;j++) {
					char base = Character.toUpperCase(binSequence.charAt(j)); 
					if(DNASequence.isInAlphabeth(base)) {
						nBases++;
						if(base == 'G' || base == 'C') {
							gcContent++;
						}
					}
				}
				if(nBases>0) {
					gcContent/=nBases;
				} else {
					gcContent = -1;
				}
				seqBins.add(new ReadDepthBin(seqName, i+1, i+binSize,gcContent));
			}
			//System.out.println("Sequence name: "+seqName+" Sequence length "+sequence.length+" end: "+end+" bins: "+seqBins.size());
		}
	}
	public void processAlignments (String filename) throws IOException {
		
		
		try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(filename)) {
			reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_MINIMAL);
			reader.setLog(log);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			reader.setFilterFlags(filterFlags);
			reader.setMinMQ(minMQ);
			Iterator<ReadAlignment> it = reader.iterator();
			//Sequence under processing
			while(it.hasNext()) {
				ReadAlignment aln = it.next();
				boolean uniqueRead = aln.isUnique();
				
				int middle = aln.getFirst()+aln.getReadLength()/2;
				List<ReadDepthBin> seqBins = bins.get(aln.getSequenceName());
				if(seqBins==null) continue;
				int binPos = middle/binSize;
				if(seqBins!=null && seqBins.size()>binPos) {
					ReadDepthBin bin = seqBins.get(binPos);
					if(!uniqueRead) bin.setInRepetitiveRegion(true);
					bin.addRead();
				}
				totalReads++;
				if(totalReads%1000000 == 0) log.info("Processed "+totalReads+" alignments");
				//if(totalReads%100000 == 0) log.info("Processing read: "+aln.getReadName()+". Location: "+aln.getSequenceName()+":"+aln.getFirst()+" flags: "+aln.getFlags()+". Unique: "+aln.isUnique()+". Bins size: "+seqBins.size()+" bin pos: "+binPos);
			}
		}
		//Set corrected depth back to raw depth
		for(List<ReadDepthBin> binsSeq:bins.values()) {
			for(ReadDepthBin bin:binsSeq) {
				bin.setCorrectedReadDepth(bin.getRawReadDepth());
			}
		}
		double sum = 0;
		int n=0;
		for(List<ReadDepthBin> binsSeq:bins.values()) {
			for(ReadDepthBin bin:binsSeq) {
				if(!bin.isInRepetitiveRegion()) {
					sum+=bin.getRawReadDepth();
					n++;
				}
			}
		}
		if(n==0 || sum/n <1) throw new IOException("The average coverage in unique regions ("+(sum/n)+") is too low for reliable CNV detection. "
				+ "Check if the XS field is present for all alignments in the bam file and if so, use the option -ignoreXS. If the average genome-wide coverage is low, then skip detection of CNVs"); 
	}
	public void correctDepthByGCContent () {
		int gcContentBins = 100;
		double [] readDepthGC = new double [gcContentBins];
		int [] gcNBins = new int [gcContentBins];
		
		Arrays.fill(readDepthGC, 0);
		Arrays.fill(gcNBins, 0);
		double globalReadDepth = 0;
		int globalNBins = 0;
		log.info("Calculating average read depth in unique bins");
		for(String seqName:bins.keySet()) {
			List<ReadDepthBin> seqBins = bins.get(seqName);
			double seqReadDepth = 0;
			int seqNBins = 0;
			for(ReadDepthBin bin:seqBins) {
				if(bin.isGoodForAverage()) {
					seqReadDepth += bin.getRawReadDepth();
					seqNBins++;
					globalReadDepth+= bin.getRawReadDepth();
					globalNBins++;
					int gcBinPos = (int)(gcContentBins*bin.getGcContent());
					if(gcBinPos == gcContentBins) gcBinPos--;
					//if(gcBinPos==0) log.info("Bin "+seqName+":"+bin.getStart()+"-"+bin.getEnd()+" has low gc Content: "+bin.getGcContent()+". Read depth: "+bin.getRawReadDepth());
					readDepthGC[gcBinPos]+=bin.getRawReadDepth();
					gcNBins[gcBinPos]++;
				}
			}
			if(seqNBins>0) {
				seqReadDepth/=seqNBins;
				if(bins.size()<100)log.info("Read depth average for "+seqName+": "+seqReadDepth+". Bins unique regions: "+seqNBins);
			} else {
				log.info("No bins in unique regions for sequence: "+seqName);
			}
		}
		if(globalNBins>0) {
			globalReadDepth/=globalNBins;
			log.info("Global read depth average: "+globalReadDepth+". Bins unique regions: "+globalNBins);
		}
		for(int i=0;i<gcNBins.length;i++) {
			if(gcNBins[i]>0) {
				readDepthGC[i]/=gcNBins[i];
			} else {
				readDepthGC[i]=globalReadDepth;
			}
			log.info("Read depth average gcContent bin "+i+": "+readDepthGC[i]+". Bins: "+gcNBins[i]);
		}
		for(String seqName:bins.keySet()) {
			if(bins.size()<100) log.info("Correcting GC for bins in sequence "+seqName);
			List<ReadDepthBin> seqBins = bins.get(seqName);
			for(ReadDepthBin bin:seqBins) {
				double gcContentBin = bin.getGcContent();
				if(gcContentBin>=0) {
					int gcBinPos = (int)(gcContentBins*gcContentBin);
					if(gcBinPos == gcContentBins) gcBinPos--;
					if(readDepthGC[gcBinPos]>0) {
						bin.setCorrectedReadDepth(bin.getRawReadDepth()*globalReadDepth/readDepthGC[gcBinPos]);
					} else {
						bin.setCorrectedReadDepth(bin.getRawReadDepth());
					}
				}
			}
		}
	}
	public void calculateReadDepthDistParameters() {
		//Calculate sample average
		double maxReadDepth = 0;
		double sum=0;
		int n=0;
		for(List<ReadDepthBin> seqBins:bins.values()) {
			for(ReadDepthBin bin:seqBins) {
				if(bin.isGoodForAverage()) {
					sum+=bin.getCorrectedReadDepth();
					if(bin.getCorrectedReadDepth()>maxReadDepth) {
						maxReadDepth = bin.getCorrectedReadDepth();
					}
					n++;
				}
			}
		}
		double average = sum/n;
		log.info("Average read depth: "+average+". Maximum: "+maxReadDepth+". Number of bins: "+n);
		int maxValueDistribution = (int) Math.round(maxReadDepth);
		int maxValueDist2 = (int) Math.round(average*4.0); 
		if(maxValueDist2> 10 && maxValueDist2<maxValueDistribution) maxValueDistribution = maxValueDist2;
		if(maxValueDistribution<2) maxValueDistribution = 2;
		//Calculate distribution
		Distribution distCalc = new Distribution(1, maxValueDistribution , 1);
		Map<String,Distribution> seqDistCalc=null;
		if(bins.size()<100) seqDistCalc =  new TreeMap<String, Distribution>(); 
		for(String seqName:bins.keySet()) {		
			List<ReadDepthBin> seqBins = bins.get(seqName);
			Distribution seqDist = null;
			if(seqDistCalc!=null) {
				seqDist = new Distribution(1, maxValueDistribution, 1);
				seqDistCalc.put(seqName, seqDist);
			}
			for(ReadDepthBin bin:seqBins) {
				if(bin.isGoodForAverage()) {
					distCalc.processDatapoint(bin.getCorrectedReadDepth());
					if(seqDist!=null) seqDist.processDatapoint(bin.getCorrectedReadDepth());
				}
			}
		}
		if(seqDistCalc!=null) printDistributions(seqDistCalc,System.out);
		
		System.out.println("Global read depth distribution");
		distCalc.printDistribution(System.out,average*2);
		double [] frequencies = distCalc.getDistribution();
		double [] midPoints = new double[frequencies.length];
		for(int i=0;i<frequencies.length;i++) {
			midPoints[i] = 1.5+i;
		}
		//Use rough estimates while I get something better
		double mode = distCalc.getMaximumBinStart()+0.5;
		sigmaReadDepth = distCalc.getEstimatedStandardDeviationPeak(mode);
		double avg = distCalc.getAverage();
		if(Math.abs(mode-avg) < sigmaReadDepth) meanReadDepth = (mode+avg)/2; 
		else meanReadDepth = average;
		//Regression regression = new Regression(midPoints,frequencies);
		/*
        GaussFunction gf = new GaussFunction();
        double[] start = new double[3];
        double[] step = new double[3];
        start[0] = distCalc.getMaximumBinStart()+0.5;
        double roughEstimateSD = distCalc.getEstimatedStandardDeviationPeak(start[0]);
        step[0] = start[0]/10;
        start[1] = roughEstimateSD;
        step[1] = roughEstimateSD/10.0;
        start[2] = distCalc.getMaximumFrequency()*roughEstimateSD*Math.sqrt(2.0*Math.PI);
        step[2] = start[2]/10.0;

        regression.simplex(gf, start, step);

        double[] bestEstimates = regression.getBestEstimates();
        meanReadDepth = bestEstimates[0];
        sigmaReadDepth = bestEstimates[1];
		 */
		
	}
	
	public List<ReadDepthBin> getBins(String seqName) {
		return Collections.unmodifiableList(bins.get(seqName));
	}
	
	public List<ReadDepthBin> getAllBins() {
		List<ReadDepthBin> answer = new ArrayList<ReadDepthBin>();
		for(String seqName:sequences.getNamesStringList()){
			answer.addAll(getBins(seqName));
		}
		return answer;
	}
	
	private void printDistributions(Map<String, Distribution> seqDistCalc,PrintStream out) {
		out.println("Read depth distribution per chromosome");
		int maxValue = 0;
		for(String seqName:seqDistCalc.keySet()) {
			out.print("\t"+seqName);
			int maxDist = (int) seqDistCalc.get(seqName).getMaxValueDistribution();
			if(maxValue <maxDist) maxValue = maxDist;
		}
		out.println();
		for(int i=0;i<maxValue;i++) {
			out.print(""+(i+1));
			for(String seqName:seqDistCalc.keySet()) {
				Distribution d = seqDistCalc.get(seqName);
				out.print("\t"+d.getDistribution()[i]);
			}
			out.println();
		}
	}
	
}
