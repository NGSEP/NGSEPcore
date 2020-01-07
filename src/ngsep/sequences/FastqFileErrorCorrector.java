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

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Stack;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import ngsep.main.CommandsDescriptor;
import ngsep.math.Distribution;
import ngsep.sequences.io.FastqFileReader;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class FastqFileErrorCorrector {
	private Logger log = Logger.getLogger(FastqFileErrorCorrector.class.getName());
	public static final int DEF_MIN_KMER_ABUNDANCE = 5;
	
	private KmersMap kmersMap;
	private int kmerSize = KmersCounter.DEFAULT_KMER_SIZE;
	private int minAbundance = DEF_MIN_KMER_ABUNDANCE;
	private int correctedErrors = 0;
	/**
	 * @return the minCount
	 */
	public int getMinAbundance() {
		return minAbundance;
	}

	/**
	 * @param minAbundance the minAbundance to set
	 */
	public void setMinAbundance(int minAbundance) {
		this.minAbundance = minAbundance;
	}
	
	/**
	 * @param minAbundance the minAbundance to set
	 */
	public void setMinAbundance(Integer minAbundance) {
		this.setMinAbundance(minAbundance.intValue());
	}

	public static void main(String[] args) throws Exception {
		FastqFileErrorCorrector instance = new FastqFileErrorCorrector();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String inFilename = args[i++];
		String outFilename = args[i++];
		instance.process(inFilename,outFilename);

	}

	public void process(String inFilename, String outFilename) throws IOException {
		correctedErrors = 0;
		buildKmersMap(inFilename);
		System.out.println("Processing file: "+inFilename);
		try (FastqFileReader reader = new FastqFileReader(inFilename);
			 OutputStream os = new GZIPOutputStream(new FileOutputStream(outFilename));
			 PrintStream out = new PrintStream(os)) {
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				processRead (read);
				read.save(out);
			}
		}
		System.out.println("Corrected "+correctedErrors+" potential errors. Output written to "+outFilename);
	}

	private void buildKmersMap(String inFilename) throws IOException {
		System.out.println("Calculating k-mers map from: "+inFilename);
		KmersCounter counter = new KmersCounter();
		counter.setLog(log);
		counter.setIgnoreLowComplexity(true);
		counter.setKmerSize(kmerSize);
		counter.setBothStrands(true);
		counter.processFile(inFilename);
		kmersMap = counter.getKmersMap();
		System.out.println("Extracted "+kmersMap.size()+" k-mers from: "+inFilename+ " (ignoring low complexity)");
		Distribution kmersDist = kmersMap.calculateAbundancesDistribution();
		int mode = (int)kmersDist.getLocalMode(5, 125);
		kmersDist.printDistributionInt(System.out);
		log.info("Distribution mode: "+mode);
		if(kmerSize>15) {
			kmersMap.filterKmers(minAbundance);
			log.info("The Map now has "+kmersMap.size()+" k-mers");
		}
	}

	public void processRead(RawRead read) {
		//processReadBestSNPChange(read);
		processReadDeBruijnExploration (read);
	}
	private void processReadDeBruijnExploration(RawRead read) {
		String readStr = read.getCharacters().toString();
		String rq = read.getQualityScores();
		StringBuilder correctedRead = new StringBuilder();
		StringBuilder correctedQualities = new StringBuilder();
		CharSequence [] readKmers = KmersCounter.extractKmers(readStr, kmerSize , true, false);
		int [] readKmerCounts = new int [readKmers.length];
		Arrays.fill(readKmerCounts, 0);
		
		for(int i=0;i<readKmers.length;i++) {
			CharSequence kmer = readKmers[i];
			if (kmer!=null) readKmerCounts[i] = kmersMap.getCount(kmer);
		}
		boolean corrected = false;
		int lastRepresented= -1;
		for(int i=0;i<readKmers.length;i++) {
			if(readKmerCounts[i] >= minAbundance) {
				if(i-1!=lastRepresented) {
					//TODO: Try to correct sequence starts
					String correctedSegment = null;
					if(lastRepresented>=0) correctedSegment = buildCorrectedSegment(lastRepresented, readKmers[lastRepresented].toString(), i, readKmers[i].toString());
					if(correctedSegment!=null  && correctedSegment.length()>=kmerSize-1 ) {
						corrected = true;
						correctedRead.append(correctedSegment);
						correctedQualities.append(rq.substring(lastRepresented+1,lastRepresented+kmerSize));
						correctedQualities.append(RawRead.generateFixedQSString('+', correctedSegment.length()-(kmerSize-1)));
					}
					else {
						correctedRead.append(readStr.substring(lastRepresented+1,i));
						correctedQualities.append(rq.substring(lastRepresented+1,i));
						
					}
				}
				correctedRead.append(readStr.charAt(i));
				correctedQualities.append(rq.charAt(i));
				lastRepresented = i;
			}
		}
		if(lastRepresented==readKmerCounts.length-1) {
			correctedRead.append(readStr.substring(lastRepresented+1));
			correctedQualities.append(rq.substring(lastRepresented+1));
		}
		else {
			String correctedSegment = buildCorrectedSegment(lastRepresented, readKmers[lastRepresented].toString(), readStr.length(), null);
			if(correctedSegment!=null && correctedSegment.length()>=kmerSize-1 ) {
				corrected = true;
				correctedRead.append(correctedSegment);
				correctedQualities.append(rq.substring(lastRepresented+1,lastRepresented+kmerSize));
				correctedQualities.append(RawRead.generateFixedQSString('+', correctedSegment.length()-(kmerSize-1)));
				
			}
			else {
				correctedRead.append(readStr.substring(lastRepresented+1));
				correctedQualities.append(rq.substring(lastRepresented+1));
			}
		}
		if(corrected) {
			read.setCharacters(correctedRead);
			read.setQualityScores(correctedQualities.toString());
		}
		
	}

	private String buildCorrectedSegment(int sourceKmerIdx, String sourceKmer, int destKmerIdx, String destKmer) {
		int kmerLength = sourceKmer.length();
		if(destKmerIdx-sourceKmerIdx<kmerLength) return null;
		int expectedAssemblyLength = destKmerIdx-sourceKmerIdx;
		if(destKmer!=null) expectedAssemblyLength+=destKmer.length();
		if(expectedAssemblyLength>3*kmerLength) return null;
		Stack<String> agenda = new Stack<>();
		agenda.push(sourceKmer);
		while (agenda.size()>0) {
			String nextState = agenda.pop();
			//Satisfability
			if((destKmer==null && nextState.length()==expectedAssemblyLength)) {
				correctedErrors++;
				return nextState.substring(1);
			}
			if(destKmer!=null && nextState.length()>destKmer.length()+2 && nextState.endsWith(destKmer)) {
				correctedErrors++;
				return nextState.substring(1, nextState.length()-destKmer.length());
			}
			//Viability
			if(nextState.length()>expectedAssemblyLength+10) continue;
			//Next states
			String kMinus1Mer = nextState.substring(nextState.length()-kmerLength+1);
			String dna = DNASequence.BASES_STRING;
			for(int i=0;i<dna.length();i++) {
				char bp = dna.charAt(i);
				String nextKmer = kMinus1Mer+bp;
				if(kmersMap.getCount(nextKmer)>=minAbundance) {
					agenda.push(nextState+bp); 
				}
			}
		}
		return null;
	}

	public void processReadBestSNPChange(RawRead read) {
		for(int h=0;h<3;h++) {
			String readStr = read.getCharacters().toString();
			char [] readChars = readStr.toCharArray();
			CharSequence [] readKmers = KmersCounter.extractKmers(readStr, kmerSize , true, false);
			int [] readKmerCounts = new int [readKmers.length];
			Arrays.fill(readKmerCounts, 0);
			
			for(int i=0;i<readKmers.length;i++) {
				CharSequence kmer = readKmers[i];
				if (kmer!=null) readKmerCounts[i] = kmersMap.getCount(kmer);
			}
			int lastRepresented= -1;
			boolean corrected = false;
			for(int i=0;i<readKmers.length;i++) {
				if(readKmerCounts[i] >= minAbundance) {
					if(i-1!=lastRepresented) {
						corrected = corrected || correctErrors (readChars,lastRepresented,i);
					}
					lastRepresented = i;
				}
			}
			corrected = corrected || correctErrors (readChars,lastRepresented,readChars.length);
			if (corrected) {
				read.setCharacters(new String(readChars));
			} else break;
		}
		
	}

	private boolean correctErrors(char [] readChars, int lastRepresented, int nextRepresented) {
		int first = 0;
		if(lastRepresented>=0) first = lastRepresented+kmerSize;
		int last = nextRepresented-1;
		if(last-first < 0) return false;
		double bestScore = getScore(readChars,lastRepresented+1,last);
		int bestI=-1;
		char bestBP = 0;
		for(int i=first;i<=last;i++) {
			char origBP = readChars[i];
			for(int j=0;j<4;j++) {
				char changeBP = DNASequence.BASES_STRING.charAt(j);
				if(changeBP != origBP) {
					readChars[i] = changeBP;
					double score = getScore(readChars,lastRepresented+1,last);
					if(score > bestScore) {
						bestScore = score;
						bestI = i;
						bestBP = changeBP;
					}
				}
			}
			readChars[i] = origBP;
		}
		if(bestI>=0) {
			readChars[bestI] = bestBP;
			correctedErrors++;
			return true;
		}
		return false;
	}

	private double getScore(char[] readChars, int first, int last) {
		CharSequence [] readKmers = KmersCounter.extractKmers(new String (readChars), kmerSize , first, last, true,false);
		double score = 0;
		for(int i=first;i<=last&& i<readKmers.length;i++) {
			CharSequence kmer = readKmers[i];
			if (kmer!=null) {
				int count = kmersMap.getCount(kmer);
				score+=count;
			}
		}
		return score;
	}

}
