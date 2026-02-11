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

import java.io.ByteArrayOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.main.ThreadPoolManager;
import ngsep.math.Distribution;
import ngsep.math.ShannonEntropyCalculator;
import ngsep.sequences.io.FastaFileReader;
import ngsep.sequences.io.FastqFileReader;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class KmersExtractor {
	
	// Constants for default values
	public static final byte DEF_KMER_LENGTH = 15;
	public static final int DEF_MIN_KMER_COUNT = 5;
	public static final int DEF_NUM_THREADS = 1;
	public static final byte INPUT_FORMAT_FASTQ=0;
	public static final byte INPUT_FORMAT_FASTA=1;
	public static final byte INPUT_FORMAT_BAM=2;
	
	private static final int MAX_LENGTH_SINGLE_TASK = 100000;
	
	// Logging and progress
	private Logger log = Logger.getLogger(KmersExtractor.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String outputPrefix = null;
	private int kmerLength = DEF_KMER_LENGTH;
	private int minKmerCount = DEF_MIN_KMER_COUNT;
	private boolean onlyForwardStrand=false;
	private byte inputFormat = INPUT_FORMAT_FASTQ;
	private boolean freeText = false;
	private boolean ignoreLowComplexity = false;
	private int numThreads = DEF_NUM_THREADS;
	private int minReadLength = 0;
	private int minReadAverageQuality = 0;
	private boolean readNCharacters = true;
	
	// Model attributes
	private KmersMap kmersMap = null;
	private boolean loadSequences = false;
	private List<QualifiedSequence> loadedSequences = null;
	
	
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
	
	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
		kmersMap=null;
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
	
	public boolean isFreeText() {
		return freeText;
	}
	public void setFreeText(boolean freeText) {
		this.freeText = freeText;
		if(freeText) this.setOnlyForwardStrand(true);
	}
	public void setFreeText(Boolean freeText) {
		this.setFreeText(freeText.booleanValue());
	}
	
	public boolean isIgnoreLowComplexity() {
		return ignoreLowComplexity;
	}
	public void setIgnoreLowComplexity(boolean ignoreLowComplexity) {
		this.ignoreLowComplexity = ignoreLowComplexity;
	}
	public void setIgnoreLowComplexity(Boolean ignoreLowComplexity) {
		this.setIgnoreLowComplexity(ignoreLowComplexity.booleanValue());
	}
	
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	public void setNumThreads(String value) {
		this.setNumThreads((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	
	public int getMinReadLength() {
		return minReadLength;
	}
	public void setMinReadLength(int minReadLength) {
		this.minReadLength = minReadLength;
	}
	
	
	public int getMinReadAverageQuality() {
		return minReadAverageQuality;
	}
	public void setMinReadAverageQuality(int minReadAverageQuality) {
		this.minReadAverageQuality = minReadAverageQuality;
	}
	/**
	 * @return the hashKmers
	 */
	public KmersMap getKmersMap() {
		return kmersMap;
	}
	
	public boolean isLoadSequences() {
		return loadSequences;
	}
	public void setLoadSequences(boolean loadSequences) {
		this.loadSequences = loadSequences;
	}
	public List<QualifiedSequence> getLoadedSequences() {
		return loadedSequences;
	}
	
	
	public boolean isReadNCharacters() {
		return readNCharacters;
	}
	public void setReadNCharacters(boolean readNCharacters) {
		this.readNCharacters = readNCharacters;
	}
	/**
	 * Receives the parameters from the command line interface and distributes the duties
	 * @param args
	 * @throws Exception 
	 */
	public static void main (String [ ] args) throws Exception {

		KmersExtractor instance = new KmersExtractor();
		//Parameters
		int k=CommandsDescriptor.getInstance().loadOptions(instance, args);
		List<String> files = new ArrayList<>();
		for(;k<args.length;k++) {
			files.add(args[k]);
		} 
		instance.processFiles(files);
		instance.saveResults();
	}
	
	/**
	 * Processes a list of input files as fasta or fastq and updates the kmers table
	 * @param files List of names of the files to process.
	 * @throws IOException If a file can not be read
	 * @throws InterruptedException 
	 */
	public void processFiles(List<String> files) throws IOException, InterruptedException {
		logParameters ();
		if(files.size()==1 && "-".equals(files.get(0))) processFastqFile(System.in);
		for(String filename:files) processFile(filename);
	}
	
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Prefix for output files:"+ outputPrefix);
		out.println("K-mer length: "+ kmerLength);
		out.println("Minimum count to save k-mer: "+ minKmerCount);
		if (onlyForwardStrand) out.println("Extract k-mers only from the forward strand");
		if (inputFormat == INPUT_FORMAT_FASTQ)  out.println("Fastq format");
		if (inputFormat == INPUT_FORMAT_FASTA)  out.println("Fasta format");
		if (ignoreLowComplexity) out.println("Ignore low complexity k-mers");
		log.info(os.toString());
		
	}
	/**
	 * Processes the file with the given name as fasta or fastq and updates the kmers table
	 * @param sequenceFileName Name of the file with the sequences to process.
	 * @throws IOException If the file can not be read
	 * @throws InterruptedException 
	 */
	public void processFile(String filename) throws IOException, InterruptedException {
		initialize();
		//Is fasta or fastq? and read it
		if(inputFormat==INPUT_FORMAT_FASTA){
			processFastaFile(filename);
		} else {
			processFastqFile(filename);	
		}
	}
	
	private void initialize() {
		if(kmersMap==null) {
			if(!isFreeText() && kmerLength<=15) kmersMap = new ShortArrayDNAKmersMapImpl((byte)kmerLength);
			else kmersMap = new DefaultKmersMapImpl();
			if(loadSequences) loadedSequences=new ArrayList<QualifiedSequence>();
		}
	}
	/**
	 * Processes the file with the given name as fastq and updates the kmers table
	 * @param filename Name of the file with the sequences to process.
	 * @throws IOException If the file can not be read
	 * @throws InterruptedException 
	 */
    public void processFastqFile(String filename) throws IOException, InterruptedException {
    	initialize();
    	ThreadPoolManager poolKmers = new ThreadPoolManager(numThreads, 100);
    	long totalLength = 0;
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			if(freeText) reader.setSequenceType(StringBuilder.class);
			else if(readNCharacters) reader.setSequenceType(DNAMaskedSequence.class);
			else reader.setSequenceType(DNASequence.class);
			if(minReadAverageQuality==0) reader.setLoadMode(FastqFileReader.LOAD_MODE_WITH_NAME);
			else reader.setMinAverageQuality(minReadAverageQuality);
			Iterator<RawRead> it = reader.iterator();
			for (int i=0;it.hasNext();i++) {
				RawRead read = it.next();
				if(read.getLength()<minReadLength) continue;
				countSequenceKmers (read, poolKmers);
				if(loadSequences) loadedSequences.add(read);
				totalLength+=read.getLength();
				if((i+1)%1000==0) log.info("Processed "+(i+1)+" sequences. Total length: "+totalLength);
			}
		}
		poolKmers.terminatePool();
	 }
    
	/**
     * Process the given input stream to get kmers count sequence by sequence
     * @param fis Input stream in fastq format
     * @throws IOException if there is an error reading the stream
	 * @throws InterruptedException 
     */
	public void processFastqFile(InputStream fis) throws IOException, InterruptedException {
		initialize();
		ThreadPoolManager poolKmers = new ThreadPoolManager(numThreads, 1000);
		try (FastqFileReader reader = new FastqFileReader(fis)) {
			if(freeText) reader.setSequenceType(StringBuilder.class);
			else if(readNCharacters) reader.setSequenceType(DNAMaskedSequence.class);
			else reader.setSequenceType(DNASequence.class);
			Iterator<RawRead> it = reader.iterator();
			for (int i=0;it.hasNext();i++) {
				RawRead read = it.next();
				if(read.getLength()<minReadLength) continue;
				countSequenceKmers (read, poolKmers);
				if(loadSequences) loadedSequences.add(read);
				if((i+1)%1000==0) log.info("Processed "+(i+1)+" sequences");
			}
		}
		poolKmers.terminatePool();
	}
	
	/**
	 * Processes the file with the given name as fasta and updates the kmers table
	 * @param filename Name of the file with the sequences to process.
	 * @throws IOException If the file can not be read
	 * @throws InterruptedException 
	 */
    public void processFastaFile(String filename) throws IOException, InterruptedException {
    	initialize();
    	ThreadPoolManager poolKmers = new ThreadPoolManager(numThreads, 1000);
    	try (FastaFileReader reader = new FastaFileReader(filename)) {
    		if(freeText) reader.setSequenceType(StringBuilder.class);
			else if(readNCharacters) reader.setSequenceType(DNAMaskedSequence.class);
			else reader.setSequenceType(DNASequence.class);
			Iterator<QualifiedSequence> it = reader.iterator();
			for (int i=0;it.hasNext();i++) {
				QualifiedSequence seq = it.next();
				if(seq.getLength()<minReadLength) continue;
				if(seq.getLength()>1000000) log.info("Processing sequence "+seq.getName());
				countSequenceKmers (seq, poolKmers);
				if(loadSequences) loadedSequences.add(seq);
				if(seq.getLength()>1000000) log.info("Processed sequence "+seq.getName()+" total k-mers: "+kmersMap.size());
				if((i+1)%1000==0) log.info("Processed "+(i+1)+" sequences");
			}
    	}
    	poolKmers.terminatePool();
	}
    public void processQualifiedSequences(List<QualifiedSequence> sequences) {
    	initialize();
    	ThreadPoolManager poolKmers = new ThreadPoolManager(numThreads, 1000);
    	int i = 0;
    	for(QualifiedSequence qseq:sequences) {
    		if(qseq.getLength()<minReadLength) continue;
    		if(qseq.getLength()>1000000) log.info("Processing sequence "+qseq.getName());
    		try {
				countSequenceKmers (qseq, poolKmers);
			} catch (InterruptedException e) {
				e.printStackTrace();
				//throw new RuntimeException("Concurrence error extracting k-mers",e);
			}
    		if(qseq.getLength()>1000000) log.info("Processed sequence "+qseq.getName()+" total k-mers: "+kmersMap.size());
    		i++;
    		if(i%100==0) log.info("Processed "+i+" sequences");
    	}
    	try {
			poolKmers.terminatePool();
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException("Concurrence error extracting k-mers",e);
		}
    }
   
    public void countSequenceKmers(QualifiedSequence qseq, ThreadPoolManager manager) throws InterruptedException {
    	if(qseq.getLength() <= MAX_LENGTH_SINGLE_TASK) {
    		manager.queueTask(()->countSequenceKmers(qseq));
    		return;
    	}
    	CharSequence seq = qseq.getCharacters();
    	int n = seq.length();
    	for(int i=0;i<n;i+=MAX_LENGTH_SINGLE_TASK) {
    		int end = i+MAX_LENGTH_SINGLE_TASK+kmerLength-1;
    		end = Math.min(end, n);
    		String subsequence = seq.subSequence(i, end).toString();
    		manager.queueTask(()->countSequenceKmers(new QualifiedSequence("",subsequence)));
    	}
    }
	public void countSequenceKmers(QualifiedSequence qseq) {
		//Forward		
		CharSequence sequence = qseq.getCharacters();
		countSequenceKmers(sequence.toString());
		//Reverse complement
		if(!onlyForwardStrand){
			CharSequence reverseSequence = DNAMaskedSequence.getReverseComplement(sequence);
			countSequenceKmers(reverseSequence.toString());
		}
	}
	/**
	 * Updates the k-mers table using the information of the given sequence
	 * @param seq CharSequence object to extract the k-mers
	 */
	public void countSequenceKmers(String seq)
	{
		initialize();
		int seqLength = seq.length();
		
		if(seqLength < kmerLength) {
			log.warning("Sequence "+seq+" smaller than k-mer length");
			return;
		}
		if(!freeText && !ignoreLowComplexity && kmerLength<=15) {
			//Faster alternative
			long [] codes = extractDNAKmerCodes(seq, kmerLength, 0, seq.length());
			synchronized (kmersMap) {
				ShortArrayDNAKmersMapImpl skmersMap = (ShortArrayDNAKmersMapImpl) kmersMap;
				for(int i=0;i<codes.length;i++) {	
					if(codes[i]>=0)skmersMap.addCodeOccurance(codes[i]);
				}
			}
			return;
		}
		String [] kmers = extractKmers(seq, kmerLength, 1, 0, seq.length(), false, freeText, ignoreLowComplexity);
		//synchronized (kmersMap) {
			for(String kmer:kmers) {
				if(kmer==null) continue;
				if(kmer.length()<=15) kmersMap.addOcurrance(kmer);
				else kmersMap.addOcurrance(pack(kmer));
			}
		//}
	}
	
	/**
	 * Extracts the k-mers present in the given sequence
	 * @param source Sequence to process. The sequence is processed as is (no uppercase or other transformation).
	 * @param kmerLength Length of the output sequences. It must be less or equal than the length of the sequence
	 * @param offset Distance between kmers
	 * @param forceLast Includes last k-mer even if it does not meet the offset requirement
	 * @param freeText Tells if the sequence should be treated as free text and then non-DNA kmers should be considered
	 * @param ignoreLowComplexity If true, ignores kmers having low complexity sequences 
	 * @return Map<Integer,String> Map of kmers indexed and sorted by (zero based) start position in the source sequence
	 */
	public static Map<Integer,String> extractKmersAsMap (String source, int kmerLength, int offset, boolean forceLast, boolean freeText, boolean ignoreLowComplexity) {
		return extractKmersAsMap(source, kmerLength, offset, 0, source.length(), forceLast, freeText, ignoreLowComplexity);
				
	}
	/**
	 * Extracts the k-mers present in the given sequence
	 * @param source Sequence to process. The sequence is processed as is (no uppercase or other transformation).
	 * @param kmerLength Length of the output sequences. It must be less or equal than the length of the sequence
	 * @param offset Distance between kmers
	 * @param start index to extract kmers
	 * @param end index to extract kmers
	 * @param forceLast Includes last k-mer even if it does not meet the offset requirement
	 * @param freeText Tells if the sequence should be treated as free text and then non-DNA kmers should be considered
	 * @param ignoreLowComplexity If true, ignores kmers having low complexity sequences 
	 * @return Map<Integer,String> Map of kmers indexed and sorted by (zero based) start position in the source sequence
	 */
	private static Map<Integer,String> extractKmersAsMap (String source, int kmerLength, int offset, int start, int end, boolean forceLast, boolean freeText, boolean ignoreLowComplexity) {
		validateLimits(source, start, end);
		int n = source.length();
		Map<Integer,String> kmersMap = new LinkedHashMap<Integer, String>();
		if (n < kmerLength) return kmersMap;
		int lastKmerStart = end - kmerLength;
		int lastPos = 0;
		for(int i = start; i <=lastKmerStart; i+=offset) {
			String kmer = source.substring(i, i+kmerLength);
			if(passFilters(kmer, freeText, ignoreLowComplexity)) kmersMap.put(i, kmer);
			lastPos = i;
		}
		if(forceLast && lastPos < lastKmerStart) {
			String kmer = source.substring(lastKmerStart, end );
			if(passFilters(kmer, freeText, ignoreLowComplexity)) kmersMap.put(lastKmerStart, kmer);
		}

		return kmersMap;
	}
	/**
	 * Extracts the k-mers present in the given sequence
	 * @param source Sequence to process. The sequence is processed as is (no uppercase or other transformation).
	 * @param kmerLength Length of the output sequences. It must be less or equal than the length of the sequence
	 * @param offset Distance between kmers
	 * @param forceLast Includes last k-mer even if it does not meet the offset requirement
	 * @param freeText Tells if the sequence should be treated as free text and then non-DNA kmers should be considered
	 * @param ignoreLowComplexity If true, ignores kmers having low complexity sequences
	 * @return CharSequence [] Array of k-mers within the source sequence. The index in the array corresponds
	 * to the index in the sequence of the start of the k-mer
	 */
	public static String [] extractKmers (String source, int kmerLength, int offset, int start, int end, boolean forceLast, boolean freeText, boolean ignoreLowComplexity) {
		validateLimits(source, start, end);
		int n = source.length();
		if(n<kmerLength) return new String[0];
		int lastKmerStart = end - kmerLength; 
		String [] kmers = new String [lastKmerStart+1];
		Arrays.fill(kmers, null);
		for(int i = start; i <=lastKmerStart; i+=offset) {
			String kmer = source.substring(i, i+kmerLength);
			if (passFilters(kmer, freeText, ignoreLowComplexity)) kmers[i]=kmer;
		}
		if(forceLast && kmers[lastKmerStart]==null) {
			String kmer = source.substring(lastKmerStart, end );
			if (passFilters(kmer, freeText, ignoreLowComplexity)) kmers[lastKmerStart]=kmer;
		}
		return kmers;
	}
	private static void validateLimits(CharSequence source, int start, int end) {
		if(start >=end) throw new IllegalArgumentException("Start index must be smaller than end index. Given start: "+start+" given end: "+end);
		if(start < 0) throw new IllegalArgumentException("Start index must be positive. Given start: "+start);
		if(end > source.length()) throw new IllegalArgumentException("End index must be at most equal to source length. Given end: "+end+" length: "+source.length());
	}
	/**
	 * Extracts the codes representing DNA kmers from the given sequence
	 * @param source Sequence to extract kmers. Usually a String but it works with StringBuilder or other types of sequences
	 * @param kmerLength must be at most 31 to allow unique encoding of DNA kmers
	 * @param start of the source sequence
	 * @param end of the source sequence
	 * @return long[] Array of kmer codes. The index corresponds to the sequence origin minus start
	 */
	public static long[] extractDNAKmerCodes (CharSequence source, int kmerLength, int start, int end) {
		if(kmerLength>31) throw new IllegalArgumentException("This method only works with kmer lengths up to 31");
		return extractKmerCodes(source, kmerLength, start, end, DNASequence.EMPTY_DNA_SEQUENCE, false);
	}
	/**
	 * Extracts the codes representing DNA kmers from the given sequence
	 * @param source Sequence to extract kmers. Usually a String but it works with StringBuilder or other types of sequences
	 * @param kmerLength must be at most 31 to allow unique encoding of DNA kmers
	 * @param start of the source sequence
	 * @param end of the source sequence
	 * @return long[] Array of kmer codes. The index corresponds to the sequence origin minus start
	 */
	public static long[] extractKmerCodes (CharSequence source, int kmerLength, int start, int end, LimitedSequence targetSeq, boolean ignoreLowComplexity) {
		validateLimits(source, start, end);
		int n = source.length();
		int nDebug = -1;
		boolean freeText = !(targetSeq instanceof DNASequence);
		//WARN: Big Maps have concurrency issues
		//List<Long> kmerCodes = new ArrayList<Long>();
		//Map<Integer,Long> kmerCodesMap = new LinkedHashMap<Integer, Long>();
		if (n < kmerLength) return new long[0];
		int rangeLength = end-start+1;
		if(rangeLength<=kmerLength) return new long[0];
		int lastKmerStart = end - kmerLength;
		long [] kmerCodesArray = new long[lastKmerStart-start+1];
		Arrays.fill(kmerCodesArray, -1);
		long lastCode = -1;
		for(int i = start; i <=lastKmerStart; i++) {
			long code;
			if(lastCode==-1) {
				CharSequence kmer = source.subSequence(i, i+kmerLength);
				if (!passFilters(kmer, freeText, ignoreLowComplexity)) continue;
				code = targetSeq.getLongCode(kmer, 0, kmerLength);
				if(n==nDebug) System.err.println("Kmer: "+kmer+" code: "+code);
			} else {
				char lastCharNextKmer = source.charAt(i+kmerLength-1);
				if(!targetSeq.isInAlphabet(lastCharNextKmer)) {
					lastCode = -1;
					//Ignore all kmers spanning the non DNA character
					i+=kmerLength-1;
					continue;
				}
				code = targetSeq.getNextCode(lastCode,kmerLength,lastCharNextKmer);
			}
			kmerCodesArray[i-start] = code;
			//kmerCodesMap.put(i, code);
			//kmerCodes.add(code);
			//lastCode = code;
		}
		return kmerCodesArray;
	}
	public static Map<Integer,Long>  extractDNAKmerCodesAsMap (CharSequence source, int kmerLength, int start, int end) {
		return extractDNAKmerCodesAsMap(source, kmerLength, start, end, false);
	}
	public static Map<Integer,Long>  extractDNAKmerCodesAsMap (CharSequence source, int kmerLength, int start, int end, boolean ignoreLowComplexity) {
		//Not good for concurrency
		Map<Integer,Long> kmerCodesMap = new LinkedHashMap<Integer, Long>();
		long [] codes = extractKmerCodes(source, kmerLength, start, end, DNASequence.EMPTY_DNA_SEQUENCE, ignoreLowComplexity);
		for(int i=0;i<codes.length;i++) {
			if(codes[i]>=0) kmerCodesMap.put(start+i, codes[i]);
		}
		return kmerCodesMap;
	}
	/*public static Map<Long, Integer> extractLocallyUniqueKmerCodes(CharSequence sequence, int kmerLength, int start, int end) {
		Map<Integer,Long> rawCodes = KmersExtractor.extractDNAKmerCodes(sequence, kmerLength, start, end);
		Map<Long, Integer> answer = new HashMap<Long, Integer>();
		Map<Long, Integer> reverseMap = new HashMap<Long,Integer>();
		Set<Integer> multiple = new HashSet<>();
		for(int i=start;i<end;i++) {
			Long code = rawCodes.get(i);
			if(code == null) {
				multiple.add(i);
				continue;
			}
			Integer previousStart = reverseMap.get(code);
			if(previousStart!=null) {
				multiple.add(i);
				continue;
			}
			reverseMap.put(code,i);
		}
		for(int i=start;i<end;i++) {
			Long code = rawCodes.get(i);
			if(code!=null && !multiple.contains(i)) {
				answer.put(code,i);
			}
		}
		return answer;
	}*/
	
	
	private static boolean passFilters(CharSequence kmer, boolean freeText, boolean ignoreLowComplexity) {
		if(!freeText && !DNASequence.isDNA(kmer)) return false;
		if(ignoreLowComplexity && isLowComplexity(kmer)) return false;
		return true;
	}
	
	public static CharSequence pack (String initialKmer) {
		CharSequence kmer = initialKmer;
		try {
			if(initialKmer.length()<=31) {
				kmer = new DNAShortKmer(kmer);
			} else {
				kmer = new DNASequence(kmer);
			}
		} catch (IllegalArgumentException e) {
		}
		return kmer;
	}
	public static final boolean isLowComplexity(CharSequence kmer) {
		ShannonEntropyCalculator calculator = new ShannonEntropyCalculator(4);
		double entropy = calculator.calculateEntropy(kmer);
		return entropy < 1.2;
	}
	public static Map<Long, List<Integer>> getReverseMap(Map<Integer, Long> kmersMap) {
		Map<Long,List<Integer>> reverseMapF = new HashMap<Long, List<Integer>>();
		for(Map.Entry<Integer, Long> entry:kmersMap.entrySet()) {
			List<Integer> posKmer = reverseMapF.computeIfAbsent(entry.getValue(), v-> new ArrayList<Integer>());
			posKmer.add(entry.getKey());
		}
		return reverseMapF;
	}
	public void saveResults () throws IOException {
		log.info("Calculating distribution of abundances from "+kmersMap.size()+" k-mers");
		Distribution kmerSpectrum = kmersMap.calculateAbundancesDistribution();
		try (PrintStream out=new PrintStream(outputPrefix+"_kmers_distribution.txt")) {
			out.println("Kmer_frequency\tNumber_of_distinct_kmers");
			kmerSpectrum.printDistributionInt(out);
		}
		kmersMap.filterKmers(minKmerCount);
		log.info("Saving "+kmersMap.size()+" filtered k-mers with minimum count "+minKmerCount);
		try (OutputStream os = new GZIPOutputStream(new FileOutputStream(outputPrefix+"_kmers.txt.gz"));
			 PrintStream out = new PrintStream(os)) {
			kmersMap.save(out);
		}
		log.info("Process finished");
	}
	public void dispose() {
		kmersMap = null;
		loadedSequences = null;
	}
}
