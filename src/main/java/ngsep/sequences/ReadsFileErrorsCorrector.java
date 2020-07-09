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
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.math.Distribution;
import ngsep.sequences.io.FastqFileReader;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class ReadsFileErrorsCorrector {
	
	// Constants for default values
	public static final int DEF_KMER_LENGTH = KmersExtractor.DEF_KMER_LENGTH;
	public static final int DEF_MIN_KMER_COUNT = KmersExtractor.DEF_MIN_KMER_COUNT;
	public static final byte INPUT_FORMAT_FASTQ=KmersExtractor.INPUT_FORMAT_FASTQ;
	public static final byte INPUT_FORMAT_FASTA=KmersExtractor.INPUT_FORMAT_FASTA;
	
	// Logging and progress
	private Logger log = Logger.getLogger(ReadsFileErrorsCorrector.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String outputFile = null;
	private String kmersMapFile = null;
	private int kmerLength = DEF_KMER_LENGTH;
	private int minKmerCount = DEF_MIN_KMER_COUNT;
	private boolean onlyForwardStrand=false;
	private byte inputFormat = INPUT_FORMAT_FASTQ;
	
	
	// Model attributes
	private KmersMap kmersMap;
	private int correctedErrors = 0;
	
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
	public String getKmersMapFile() {
		return kmersMapFile;
	}
	public void setKmersMapFile(String kmersMapFile) {
		this.kmersMapFile = kmersMapFile;
	}
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
	
	public int getMinKmerCount() {
		return minKmerCount;
	}
	public void setMinKmerCount(int minKmerCount) {
		this.minKmerCount = minKmerCount;
	}
	public void setMinKmerCount(String value) {
		setMinKmerCount((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public boolean isOnlyForwardStrand() {
		return onlyForwardStrand;
	}
	public void setOnlyForwardStrand(boolean onlyForwardStrand) {
		this.onlyForwardStrand = onlyForwardStrand;
	}
	public void setOnlyForwardStrand(Boolean onlyForwardStrand) {
		this.setOnlyForwardStrand(onlyForwardStrand.booleanValue());
	}
	
	public byte getInputFormat() {
		return inputFormat;
	}
	public void setInputFormat(byte inputFormat) {
		if(inputFormat!=INPUT_FORMAT_FASTA && inputFormat!=INPUT_FORMAT_FASTQ) throw new IllegalArgumentException("Invalid input format "+inputFormat);
		this.inputFormat = inputFormat;
	}
	public void setInputFormat(String value) {
		this.setInputFormat((byte) OptionValuesDecoder.decode(value, Byte.class));
	}

	public static void main(String[] args) throws Exception {
		ReadsFileErrorsCorrector instance = new ReadsFileErrorsCorrector();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
		

	}

	public void run() throws IOException {
		logParameters();
		if(inputFile==null) {
			log.severe("The input file is mandatory");
			return;
		}
		if(outputFile==null) {
			log.severe("The output file is mandatory");
			return;
		}
		process(inputFile,outputFile);
		log.info("Process finished");
	}
	
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		if(kmersMapFile!=null) out.println("K-mers map file: "+ kmersMapFile);
		out.println("K-mer length: "+ kmerLength);
		out.println("Minimum count to save k-mer: "+ minKmerCount);
		if (onlyForwardStrand) out.println("Extract k-mers only from the forward strand");
		if (inputFormat == INPUT_FORMAT_FASTQ)  out.println("Fastq format");
		if (inputFormat == INPUT_FORMAT_FASTA)  out.println("Fasta format");
		log.info(os.toString());
		
	}
	public void process(String inFilename, String outFilename) throws IOException {
		correctedErrors = 0;
		if (kmersMapFile!=null) loadKmersMap();
		else buildKmersMap(inFilename);
		Distribution kmersDist = kmersMap.calculateAbundancesDistribution();
		int mode = (int)kmersDist.getLocalMode(5, 125);
		//kmersDist.printDistributionInt(System.out);
		log.info("Distribution mode: "+mode);
		kmersMap.filterKmers(minKmerCount);
		log.info("The Map now has "+kmersMap.size()+" k-mers");
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

	private void loadKmersMap() throws IOException {
		log.info("Loading k-mers map from : "+kmersMapFile);
		if(kmerLength<=15) kmersMap = new ByteArrayKmersMapImpl((byte) kmerLength);
		else kmersMap = new DefaultKmersMapImpl();
		try (FileInputStream fis = new FileInputStream(kmersMapFile)) {
			InputStream is=fis;
			if(kmersMapFile.toLowerCase().endsWith(".gz")) {
				is = new ConcatGZIPInputStream(is);
			}
			try (BufferedReader in = new BufferedReader(new InputStreamReader(is))) {
				String line = in.readLine();
				while(line!=null) {
					String [] items = line.split("\t| ");
					String kmer = items[0];
					int count = Integer.parseInt(items[1]);
					kmersMap.setCount(kmer,count);
					line = in.readLine();
				}
			}
		}
		System.out.println("Extracted "+kmersMap.size()+" k-mers from: " + kmersMapFile);
		
	}
	private void buildKmersMap(String inFilename) throws IOException {
		log.info("Calculating k-mers map from reads in : "+inFilename);
		KmersExtractor counter = new KmersExtractor();
		counter.setLog(log);
		counter.setIgnoreLowComplexity(false);
		counter.setKmerLength(kmerLength);
		counter.setOnlyForwardStrand(onlyForwardStrand);
		counter.processFile(inFilename);
		kmersMap = counter.getKmersMap();
		System.out.println("Extracted "+kmersMap.size()+" k-mers from: " + inFilename);
		
	}

	public void processRead(RawRead read) {
		// TODO: Option to choose algorithm
		//processReadBestSNPChange(read);
		processReadDeBruijnExploration (read);
	}
	public void processReadDeBruijnExploration(RawRead read) {
		String readStr = read.getCharacters().toString();
		String rq = read.getQualityScores();
		StringBuilder correctedRead = new StringBuilder();
		StringBuilder correctedQualities = new StringBuilder();
		CharSequence [] readKmers = KmersExtractor.extractKmers(readStr, kmerLength , 1, 0, readStr.length(), false, true, false);
		int [] readKmerCounts = new int [readKmers.length];
		Arrays.fill(readKmerCounts, 0);
		
		for(int i=0;i<readKmers.length;i++) {
			CharSequence kmer = readKmers[i];
			if (kmer!=null) readKmerCounts[i] = kmersMap.getCount(kmer);
		}
		boolean corrected = false;
		int lastRepresented= -1;
		for(int i=0;i<readKmers.length;i++) {
			if(readKmerCounts[i] >= minKmerCount) {
				if(i-1!=lastRepresented) {
					//TODO: Try to correct sequence starts
					int regionLength = i-lastRepresented-1;
					//if(lastRepresented>=0) System.out.println("Trying to correct from "+lastRepresented+" to "+i+" Length: "+regionLength+" last kmer: "+readKmers[lastRepresented]+" count: "+readKmerCounts[lastRepresented]+" next kmer: "+readKmers[i]+" count: "+readKmerCounts[i]);
					
					String correctedSegment = null;
					if(lastRepresented>=0) correctedSegment = buildCorrectedSegment(lastRepresented, readKmers[lastRepresented].toString(), i, readKmers[i].toString());
					//System.out.println("Corrected segment "+correctedSegment);
					if(correctedSegment!=null ) {
						int segmentLength = correctedSegment.length();
						//System.out.println("Corrected segment length "+segmentLength);
						corrected = true;
						correctedRead.append(correctedSegment);
						if(segmentLength>=regionLength) {
							correctedQualities.append(rq.substring(lastRepresented+1,i));
							if (segmentLength>regionLength) correctedQualities.append(RawRead.generateFixedQSString('+', segmentLength-regionLength));
						} else {
							int endCorrected = lastRepresented+1+segmentLength;
							correctedQualities.append(rq.substring(lastRepresented+1,endCorrected));
						}
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
		if(expectedAssemblyLength>4*kmerLength) return null;
		//Stack<String> agenda = new Stack<>();
		PriorityQueue<String> agenda = new PriorityQueue<String>(new Comparator<String>() {
			@Override
			public int compare(String state1, String state2) {
				int score1 = getScore (state1, destKmer);
				int score2 = getScore (state2, destKmer);
				return score2-score1;
			}	
		});
		
		agenda.add(sourceKmer);
		for (int candidates = 0;agenda.size()>0 && candidates < 5000;candidates++) {
			String nextState = agenda.remove();
			//if(sourceKmerIdx==2528) System.out.println("Next state: "+nextState+" length: "+nextState.length()+" score: "+getScore(nextState, destKmer)+" agenda size: "+agenda.size());
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
			if(nextState.length()>expectedAssemblyLength+5) continue;
			//Next states
			String kMinus1Mer = nextState.substring(nextState.length()-kmerLength+1);
			String dna = DNASequence.BASES_STRING;
			for(int i=0;i<dna.length();i++) {
				char bp = dna.charAt(i);
				String nextKmer = kMinus1Mer+bp;
				if(kmersMap.getCount(nextKmer)>=minKmerCount) {
					agenda.add(nextState+bp); 
				}
			}
		}
		return null;
	}
	
	private static int getScore(String state, String destKmer) {
		if(destKmer==null) return 0;
		for(int i=destKmer.length();i>0;i--) {
			if(state.endsWith(destKmer.substring(0,i))) return i;
		}
		return 0;
	}

	public void processReadBestSNPChange(RawRead read) {
		for(int h=0;h<3;h++) {
			String readStr = read.getCharacters().toString();
			char [] readChars = readStr.toCharArray();
			CharSequence [] readKmers = KmersExtractor.extractKmers(readStr, kmerLength , 1, 0, readStr.length(), false, true, false);
			int [] readKmerCounts = new int [readKmers.length];
			Arrays.fill(readKmerCounts, 0);
			
			for(int i=0;i<readKmers.length;i++) {
				CharSequence kmer = readKmers[i];
				if (kmer!=null) readKmerCounts[i] = kmersMap.getCount(kmer);
			}
			int lastRepresented= -1;
			boolean corrected = false;
			for(int i=0;i<readKmers.length;i++) {
				if(readKmerCounts[i] >= minKmerCount) {
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
		CharSequence [] segmentKmers = KmersExtractor.extractKmers(segment, kmerLength , 1, 0, segment.length(), false, true, false);
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


