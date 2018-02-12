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
import java.util.Iterator;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import ngsep.main.CommandsDescriptor;
import ngsep.sequences.io.FastqFileReader;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class FastqFileErrorCorrector {
	private Logger log = Logger.getLogger(FastqFileErrorCorrector.class.getName());
	private KmersMap kmersMap;
	private int kmerSize = KmersCounter.DEFAULT_KMER_SIZE;
	private int minAbundance = 5;
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
		System.out.println("Calculating k-mers map from: "+inFilename);
		KmersCounter counter = new KmersCounter();
		counter.setLog(log);
		counter.setKmerSize(kmerSize);
		counter.processFile(inFilename);
		kmersMap = counter.getKmersMap();
		log.info("Filtering from "+kmersMap.size()+" k-mers by minimum abundance: "+minAbundance);
		kmersMap.filterKmers(minAbundance);
		log.info("The Map now has "+kmersMap.size()+" k-mers");
		kmersMap = counter.getKmersMap();
		kmerSize = counter.getKmerSize();
		System.out.println("Extracted "+kmersMap.size()+" filtered k-mers from: "+inFilename);
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

	public void processRead(RawRead read) {
		for(int h=0;h<3;h++) {
			String readStr = read.getCharacters().toString();
			char [] readChars = readStr.toCharArray();
			CharSequence [] readKmers = KmersCounter.extractKmers(readStr, kmerSize , true);
			int [] readKmerCounts = new int [readKmers.length];
			
			for(int i=0;i<readKmers.length;i++) {
				CharSequence kmer = readKmers[i];
				if (kmer==null) readKmerCounts[i] = 0;
				else {
					readKmerCounts[i] = kmersMap.getCount(kmer);
				}
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
		CharSequence [] readKmers = KmersCounter.extractKmers(new String (readChars), kmerSize , first, last, true);
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
