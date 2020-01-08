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

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Stack;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.math.Distribution;
import ngsep.sequences.io.FastqFileReader;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class ReadsFileErrorsCorrector {
	private Logger log = Logger.getLogger(ReadsFileErrorsCorrector.class.getName());
	public static final int DEF_KMER_LENGTH = KmersCounter.DEF_KMER_LENGTH;
	public static final byte INPUT_FORMAT_FASTQ=KmersCounter.INPUT_FORMAT_FASTQ;
	public static final byte INPUT_FORMAT_FASTA=KmersCounter.INPUT_FORMAT_FASTA;
	public static final int DEF_MIN_KMER_ABUNDANCE = 5;
	
	private KmersMap kmersMap;
	private int kmerLength = DEF_KMER_LENGTH;
	private boolean countOnlyForwardStrand=false;
	private byte inputFormat = INPUT_FORMAT_FASTQ;
	private int minAbundance = DEF_MIN_KMER_ABUNDANCE;
	private int correctedErrors = 0;
	
	
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
		if(kmerLength<=15) kmersMap = new ByteArrayKmersMapImpl((byte) kmerLength);
		else kmersMap = new DefaultKmersMapImpl();
	}
	public void setKmerLength(String value) {
		setKmerLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public boolean isCountOnlyForwardStrand() {
		return countOnlyForwardStrand;
	}
	public void setCountOnlyForwardStrand(boolean countOnlyForwardStrand) {
		this.countOnlyForwardStrand = countOnlyForwardStrand;
	}
	public void setCountOnlyForwardStrand(Boolean countOnlyForwardStrand) {
		this.setCountOnlyForwardStrand(countOnlyForwardStrand.booleanValue());
	}
	/**
	 * @return the inputFormat
	 */
	public byte getInputFormat() {
		return inputFormat;
	}

	/**
	 * @param inputFormat the inputFormat to set
	 */
	public void setInputFormat(byte inputFormat) {
		if(inputFormat!=INPUT_FORMAT_FASTA && inputFormat!=INPUT_FORMAT_FASTQ) throw new IllegalArgumentException("Invalid input format "+inputFormat);
		this.inputFormat = inputFormat;
	}

	public void setInputFormat(String value) {
		this.setInputFormat((byte) OptionValuesDecoder.decode(value, Byte.class));
	}
	
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
		ReadsFileErrorsCorrector instance = new ReadsFileErrorsCorrector();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String inFilename = args[i++];
		String outFilename = args[i++];
		instance.process(inFilename,outFilename);

	}

	public void process(String inFilename, String outFilename) throws IOException {
		correctedErrors = 0;
		buildKmersMap(inFilename);
		System.out.println("Processing file: "+inFilename);
		int numReads=0;
		long numBp = 0;
		long mbp = 0;
		if(inputFormat==INPUT_FORMAT_FASTQ) {
			try (FastqFileReader reader = new FastqFileReader(inFilename);
				 OutputStream os = new GZIPOutputStream(new FileOutputStream(outFilename));
				 PrintStream out = new PrintStream(os)) {
				Iterator<RawRead> it = reader.iterator();
				while (it.hasNext()) {
					RawRead read = it.next();
					processRead (read);
					read.save(out);
					numReads++;
					numBp+=read.getLength();
					if(mbp<numBp/1000000) {
						mbp = numBp/1000000;
						log.info("Processed "+numReads+" reads and "+mbp+" Mbp. Corrected "+correctedErrors+" potential errors");
					}
				}
			}
		} else if (inputFormat==INPUT_FORMAT_FASTA) {
			try (FileReader reader = new FileReader(inFilename);
				 BufferedReader in = new BufferedReader(reader);
				 OutputStream os = new GZIPOutputStream(new FileOutputStream(outFilename));
				 PrintStream out = new PrintStream(os)) {
				 String line = in.readLine();
				 while (line!=null) {
					String readName = line.substring(1);
					String readSeq = in.readLine();
					RawRead read = new RawRead(readName, readSeq, RawRead.generateFixedQSString('5', readSeq.length()));
					processRead (read);
					read.save(out);
					numReads++;
					numBp+=read.getLength();
					if(mbp<numBp/1000000) {
						mbp = numBp/1000000;
						log.info("Processed "+numReads+" reads and "+mbp+" Mbp. Corrected "+correctedErrors+" potential errors");
					}
					line = in.readLine();	
				}
			}
		}
		log.info("Processed "+numReads+" reads and "+mbp+" Mbp. Corrected "+correctedErrors+" potential errors. Output written to "+outFilename);
	}

	private void buildKmersMap(String inFilename) throws IOException {
		System.out.println("Calculating k-mers map from: "+inFilename);
		KmersCounter counter = new KmersCounter();
		counter.setLog(log);
		counter.setIgnoreLowComplexity(true);
		counter.setKmerLength(kmerLength);
		counter.setCountOnlyForwardStrand(countOnlyForwardStrand);
		counter.processFile(inFilename);
		kmersMap = counter.getKmersMap();
		System.out.println("Extracted "+kmersMap.size()+" k-mers from: "+inFilename+ " (ignoring low complexity)");
		Distribution kmersDist = kmersMap.calculateAbundancesDistribution();
		int mode = (int)kmersDist.getLocalMode(5, 125);
		kmersDist.printDistributionInt(System.out);
		log.info("Distribution mode: "+mode);
		if(kmerLength>15) {
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
		CharSequence [] readKmers = KmersCounter.extractKmers(readStr, kmerLength , 1, false, true, false);
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
					if(correctedSegment!=null  && correctedSegment.length()>=kmerLength-1 ) {
						corrected = true;
						correctedRead.append(correctedSegment);
						correctedQualities.append(rq.substring(lastRepresented+1,lastRepresented+kmerLength));
						correctedQualities.append(RawRead.generateFixedQSString('+', correctedSegment.length()-(kmerLength-1)));
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
		if(lastRepresented==-1 || lastRepresented==readKmerCounts.length-1) {
			correctedRead.append(readStr.substring(lastRepresented+1));
			correctedQualities.append(rq.substring(lastRepresented+1));
		}
		else {
			String correctedSegment = buildCorrectedSegment(lastRepresented, readKmers[lastRepresented].toString(), readStr.length(), null);
			if(correctedSegment!=null && correctedSegment.length()>=kmerLength-1 ) {
				corrected = true;
				correctedRead.append(correctedSegment);
				correctedQualities.append(rq.substring(lastRepresented+1,lastRepresented+kmerLength));
				correctedQualities.append(RawRead.generateFixedQSString('+', correctedSegment.length()-(kmerLength-1)));
				
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
			CharSequence [] readKmers = KmersCounter.extractKmers(readStr, kmerLength , 1, false, true, false);
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
		if(lastRepresented>=0) first = lastRepresented+kmerLength;
		int last = nextRepresented-1;
		if(last-first < 0) return false;
		double bestScore = getScore(readChars,lastRepresented+1,last);
		int bestI=-1;
		char bestBP = 0;
		for(int i=first;i<=last;i++) {
			char origBP = readChars[i];
			//Simulate all single bp changes
			for(int j=0;j<DNASequence.BASES_STRING.length();j++) {
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
		String segment = (new String(readChars)).substring(first,last+1);
		CharSequence [] segmentKmers = KmersCounter.extractKmers(segment, kmerLength , 1, false, true, false);
		double score = 0;
		for(int i=0;i<segmentKmers.length;i++) {
			CharSequence kmer = segmentKmers[i];
			if (kmer!=null) {
				int count = kmersMap.getCount(kmer);
				score+=count;
			}
		}
		return score;
	}

}
