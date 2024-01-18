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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment.Platform;
import ngsep.alignments.io.ReadAlignmentFileWriter;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.ReferenceGenomeFMIndex;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.main.ThreadPoolManager;
import ngsep.main.io.ParseUtils;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaFileReader;
import ngsep.sequences.io.FastqFileReader;

/**
 * Program to align reads to a reference genome
 * @author German Andrade
 * @author Jorge Duitama
 */
public class ReadsAligner {

	// Constants for default values
	
	public static final String DEF_SAMPLE_ID = "Sample";
	public static final byte INPUT_FORMAT_FASTQ=KmersExtractor.INPUT_FORMAT_FASTQ;
	public static final byte INPUT_FORMAT_FASTA=KmersExtractor.INPUT_FORMAT_FASTA;
	public static final int DEF_MAX_ALNS_PER_READ=1;
	public static final int DEF_KMER_LENGTH = 25;
	public static final int DEF_WINDOW_LENGTH = 40;
	public static final int DEF_MIN_INSERT_LENGTH=0;
	public static final int DEF_MAX_INSERT_LENGTH=1000;
	public static final int DEF_NUM_THREADS=1;
	

	public static final int MAX_SPACE_BETWEEN_KMERS = 50;
	
	// Logging and progress
	private Logger log = Logger.getLogger(ReadsAligner.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	// Parameters
	private String inputFile = null;
	private String inputFile2 = null;
	private String outputFile = null;
	private String fmIndexFile = null;
	private String knownSTRsFile = null;
	private String sampleId = DEF_SAMPLE_ID;
	private ReadAlignment.Platform platform = null;
	private byte inputFormat = INPUT_FORMAT_FASTQ;
	private int maxAlnsPerRead = DEF_MAX_ALNS_PER_READ;
	private int kmerLength = DEF_KMER_LENGTH;
	private int minInsertLength = DEF_MIN_INSERT_LENGTH;
	private int maxInsertLength = DEF_MAX_INSERT_LENGTH;
	private int windowLength = DEF_WINDOW_LENGTH;
	
	
	private int numThreads = DEF_NUM_THREADS;
	
	// Model attributes
	private ReferenceGenome genome;
	private ReferenceGenomeFMIndex fMIndex=null;
	private ReadAlignmentObjectsFactory factory;
	
	
	private ThreadPoolManager pool;
	
	// Statistics
	private int totalReads = 0;
	private int readsAligned = 0;
	private int numProperPairs = 0;
	private int numNonProperPairs = 0;
	private int numAlignedSingle = 0;
	private int uniqueAlignments=0;
	
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
	
	public String getInputFile2() {
		return inputFile2;
	}
	public void setInputFile2(String inputFile2) {
		this.inputFile2 = inputFile2;
		this.setPlatform(ReadAlignment.Platform.ILLUMINA);
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

	public ReferenceGenomeFMIndex getFmIndex() {
		return fMIndex;
	}
	public void setFmIndex(ReferenceGenomeFMIndex fMIndex) {
		this.fMIndex = fMIndex;
	}
	public String getFmIndexFile() {
		return fmIndexFile;
	}
	public void setFmIndexFile(String fmIndexFile) {
		this.fmIndexFile = fmIndexFile;
	}

	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public String getKnownSTRsFile() {
		return knownSTRsFile;
	}
	public void setKnownSTRsFile(String knownSTRsFile) {
		this.knownSTRsFile = knownSTRsFile;
		//TODO: Find this method
		//if(knownSTRsFile!=null && !knownSTRsFile.isEmpty()) loadSTRsFile(knownSTRsFile);
	}
	
	public String getSampleId() {
		return sampleId;
	}
	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}
	
	public ReadAlignment.Platform getPlatform() {
		return platform;
	}
	public void setPlatform(ReadAlignment.Platform platform) {
		this.platform = platform;
	}
	public void setPlatform(String value) {
		ReadAlignment.Platform platform = ReadAlignment.Platform.valueOf(value);
		this.setPlatform(platform);
	}
	
	
	public byte getInputFormat() {
		return inputFormat;
	}
	public void setInputFormat(byte inputFormat) {
		this.inputFormat = inputFormat;
	}
	public void setInputFormat(String value) {
		setInputFormat((byte)OptionValuesDecoder.decode(value, Byte.class));
	}
	
	public int getMaxAlnsPerRead() {
		return maxAlnsPerRead;
	}
	public void setMaxAlnsPerRead(int maxAlnsPerRead) {
		this.maxAlnsPerRead = maxAlnsPerRead;
	}
	public void setMaxAlnsPerRead(String value) {
		setMaxAlnsPerRead((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
	}
	public void setKmerLength(String value) {
		setKmerLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getMinInsertLength() {
		return minInsertLength;
	}
	public void setMinInsertLength(int minInsertLength) {
		this.minInsertLength = minInsertLength;
	}
	public void setMinInsertLength(String value) {
		setMinInsertLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getMaxInsertLength() {
		return maxInsertLength;
	}
	public void setMaxInsertLength(int maxInsertLength) {
		this.maxInsertLength = maxInsertLength;
	}
	public void setMaxInsertLength(String value) {
		setMaxInsertLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getWindowLength() {
		return windowLength;
	}
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}
	public void setWindowLength(String value) {
		setWindowLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	public void setNumThreads(String value) {
		setNumThreads((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public static void main(String[] args) throws Exception 
	{
		ReadsAligner instance = new ReadsAligner();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run () throws IOException {
		long time = System.currentTimeMillis();
		logParameters ();
		if(genome==null) throw new IOException("The reference genome is a required parameter");
		if(fmIndexFile!=null) {
			if (fmIndexFile!=null) {
				log.info("Loading reference index from file: "+fmIndexFile);
				fMIndex = ReferenceGenomeFMIndex.load(genome, fmIndexFile);
			}
		}
		initializeFactory();
		QualifiedSequenceList sequences = genome.getSequencesMetadata();
		
		pool = new ThreadPoolManager(numThreads, platform.isLongReads()?100:10000);
		boolean paired = false;
		PrintStream out = System.out;
		if(outputFile!=null) out = new PrintStream(outputFile); 
		try (ReadAlignmentFileWriter writer = new ReadAlignmentFileWriter(sequences, out)){
			writer.setSampleInfo(sampleId, platform);
			if(inputFile!=null && inputFile2!=null) {
				log.info("Aligning paired end reads from files: "+inputFile + " and "+inputFile2);
				paired = true;
				alignReads(inputFile,inputFile2, writer);
			} else if (inputFile!=null) {
				//if(longReads && inputFile2!=null) log.warning("Paired end alignment not supported for long reads. Ignoring file "+ inputFile2);
				log.info("Aligning single reads from file: "+inputFile);
				alignReads(inputFile, writer);
			} else if (inputFile2!=null ) {
				throw new IOException("The first input file is required for paired end alignment");
			} else {
				log.info("Aligning single reads from standard input");
				alignReads(System.in, writer);
			}
			pool.terminatePool();
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		log.info("Time: "+seconds+" seconds");
		printStatistics(paired);
		
	}
	
	//Constructors
	public ReadsAligner() {
	}
	
	public ReadsAligner(ReferenceGenome genome, Platform platform) {
		this.genome = genome;
		this.platform = platform;
		//if (platform.isLongReads()) initializeLongReadsFactory();
	}
	
	private void initializeFactory() {
		if(factory!=null) return;
		if(platform==null) {
			if(inputFile!=null)
				try {
					guessPlatform(inputFile);
				} catch (IOException e) {
					throw new RuntimeException("IO error guessing platform",e);
				}
			else throw new RuntimeException("Null platform. The platform must be initialized before mapping reads");
		}
		factory = new ReadAlignmentObjectsFactory(genome);
		factory.setLog(log);
		factory.setKmerLength(kmerLength);
		factory.setWindowLength(windowLength);
		factory.setNumThreads(numThreads);
		factory.setPlatform(platform);
		factory.setFmIndex(fMIndex);
		if(platform.isLongReads()) factory.setAlignmentAlgorithm(UngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_DYNAMIC_KMERS);
		//if(platform.isLongReads()) factory.setAlignmentAlgorithm(UngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_STATIC_BAND);
		//if(platform.isLongReads()) factory.setAlignmentAlgorithm(UngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_SIMPLE_GAP);
		else factory.setAlignmentAlgorithm(UngappedSearchHitsClusterAligner.ALIGNMENT_ALGORITHM_SHORT_READS);
		factory.requestClustersFinder();
		//System.out.println("Long reads: "+platform.isLongReads());
		factory.requestAligner();
		log.info("Initialized aligner");
	}
	
	private void guessPlatform(String inputFile) throws IOException {
		log.info("Guessing platform");
		int sampleSize = 1000;
		int maxLength = 0;
		int sumQS = 0;
		int numQSSampled = 0;
		if(inputFormat == INPUT_FORMAT_FASTQ) {
			try (FastqFileReader reader = new FastqFileReader(inputFile)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<RawRead> it = reader.iterator();
				for(int i=1;it.hasNext() && i<sampleSize;i++) {
					RawRead read = it.next();
					String qs = read.getQualityScores();
					maxLength = Math.max(maxLength, read.getLength());
					for(int j=10;j<qs.length() && j<200;j+=20) {
						int nextVal = qs.charAt(j)-33;
						sumQS+=nextVal;
						numQSSampled++;
					}
				}
			}
		}  else if(inputFormat== INPUT_FORMAT_FASTA) {
			try (FastaFileReader reader = new FastaFileReader(inputFile)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<QualifiedSequence> it = reader.iterator();
				for(int i=1;it.hasNext() && i<sampleSize;i++) {
					QualifiedSequence seq = it.next();
					RawRead read = new RawRead(seq.getName(), seq.getCharacters(),null);
					maxLength = Math.max(maxLength, read.getLength());
				}
			}
		}
		platform = ReadAlignment.Platform.PACBIO;
		int avgQS = 100;
		if(numQSSampled>0) avgQS = sumQS/numQSSampled;
		if(maxLength<400 && avgQS > 25 ) platform = ReadAlignment.Platform.ILLUMINA;
		else if (maxLength>30000 && avgQS < 25) platform = ReadAlignment.Platform.ONT;
		log.info("Max length: "+maxLength+" "+" average quality score: "+avgQS+" Guessed platform: "+platform);
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		if (inputFile2!=null) out.println("Second file with paired-end reads  "+inputFile2);
		out.println("Output file:"+ outputFile);
		if (genome!=null) out.println("Reference genome loaded from file: "+genome.getFilename());
		if (fmIndexFile!=null) out.println("FM index file "+fmIndexFile);
		out.println("Sample id: "+ sampleId);
		if(platform!=null) out.println("Platform: "+ platform);
		out.println("K-mer length: "+ kmerLength);
		if (inputFormat == INPUT_FORMAT_FASTQ)  out.println("Fastq format");
		if (inputFormat == INPUT_FORMAT_FASTA)  out.println("Fasta format");
		if (knownSTRsFile!=null) out.println("Fie with known short tandem repeats "+knownSTRsFile);
		out.println("Maximum alignments per read: "+ maxAlnsPerRead);
		out.println("Window length to calculate minimizers: "+ windowLength);
		if (inputFile2!=null) {
			out.println("Proper limits for paired-end alignment. Minimum: "+ minInsertLength+" maximum: "+maxInsertLength);
		}
		out.println("Number of threads: "+ numThreads);
		log.info(os.toString());
	}
	
	/**
	 * Aligns readsFile to the reference genome
	 * @param readsFile Fastq file with the reads to align
	 * @param writer
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void alignReads( String readsFile, ReadAlignmentFileWriter writer) throws IOException, InterruptedException {
		if(platform==null) guessPlatform(readsFile);
		if(inputFormat == INPUT_FORMAT_FASTQ) {
			try (FastqFileReader reader = new FastqFileReader(readsFile)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<RawRead> it = reader.iterator();
				for(int i=1;it.hasNext();i++) {
					RawRead read = it.next();
					final int readNumber = i;
					pool.queueTask( () -> processSingleRead(readNumber, read, writer));
				}
			}
		} else if(inputFormat== INPUT_FORMAT_FASTA) {
			try (FastaFileReader reader = new FastaFileReader(readsFile)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<QualifiedSequence> it = reader.iterator();
				for(int i=1;it.hasNext();i++) {
					QualifiedSequence seq = it.next();
					RawRead read = new RawRead(seq.getName(), seq.getCharacters(),null);
					final int readNumber = i;
					pool.queueTask( () -> processSingleRead(readNumber, read, writer));
				}
			}
		}
	}
	private void alignReads(InputStream in, ReadAlignmentFileWriter writer) throws IOException, InterruptedException {
		if(inputFormat == INPUT_FORMAT_FASTQ) {
			try (FastqFileReader reader = new FastqFileReader(in)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<RawRead> it = reader.iterator();
				for(int i=1;it.hasNext();i++) {
					RawRead read = it.next();
					final int readNumber = i;
					pool.queueTask( () -> processSingleRead(readNumber, read, writer));
				}
			}
		}  else if(inputFormat== INPUT_FORMAT_FASTA) {
			try (FastaFileReader reader = new FastaFileReader(in)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<QualifiedSequence> it = reader.iterator();
				for(int i=1;it.hasNext();i++) {
					QualifiedSequence seq = it.next();
					RawRead read = new RawRead(seq.getName(), seq.getCharacters(),null);
					final int readNumber = i;
					pool.queueTask( () -> processSingleRead(readNumber, read, writer));
				}
			}
		}
	}
	

	/**
	 * Aligns readsFile with the fMIndexFile for pairend
	 * @param fMIndexFile Binary file with the serialization of an FMIndex
	 * @param readsFile1 Fastq file with the reads to align
	 * @param readsFile2 Fastq file with the reads to align
	 * @param out
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void alignReads( String readsFile1, String readsFile2, ReadAlignmentFileWriter writer) throws IOException, InterruptedException {
		if(platform==null) platform = ReadAlignment.Platform.ILLUMINA;
		try (FastqFileReader reader1 = new FastqFileReader(readsFile1);
			 FastqFileReader reader2 = new FastqFileReader(readsFile2)) {
			reader1.setSequenceType(DNAMaskedSequence.class);
			reader2.setSequenceType(DNAMaskedSequence.class);
			Iterator<RawRead> it1 = reader1.iterator();
			Iterator<RawRead> it2 = reader2.iterator();
			for(int i=1;it1.hasNext() && it2.hasNext();i++) {
				RawRead read1 = it1.next();
				RawRead read2 = it2.next();
				final int readNumber = i;
				pool.queueTask(()->processPairedEndRead(readNumber, read1, read2, writer));
			}
		}
	}
	
	private void processSingleRead(int readNumber, RawRead read, ReadAlignmentFileWriter writer) {
		List<ReadAlignment> alns = alignRead(read);
		//System.out.println("Alignments for: "+read.getName()+" "+alns.size());
		int numAlns = alns.size();
		if (alns.size()==0) {
			alns.add(ReadAlignment.createMockAlignmentUnmappedRead(read, false, false));
		}
		synchronized (writer) {
			for(ReadAlignment aln:alns) writer.write(aln);
			totalReads++;
			if(numAlns>0) readsAligned++;
			if(numAlns==1 && alns.get(0).getAlignmentQuality()>20) uniqueAlignments++;
		}
		checkProgress(readNumber);
	}
	public List<ReadAlignment> alignRead (QualifiedSequence read) {
		initializeFactory();
		SingleReadsAligner aligner = new SingleReadsAligner(genome, factory.requestClustersFinder(), factory.requestAligner());
		return aligner.alignRead(read);
	}
	
	private void processPairedEndRead (int readNumber, RawRead read1, RawRead read2, ReadAlignmentFileWriter writer) {
		List<ReadAlignment> alns = alignPairedEndReads(read1, read2,true);
		synchronized (writer) {
			for(ReadAlignment aln:alns) writer.write(aln);
		}
		checkProgress(readNumber);
	}
	public List<ReadAlignment> alignPairedEndReads(RawRead read1, RawRead read2, boolean createUnmappedReadRecords) {
		initializeFactory();
		PairedReadsAligner aligner = new PairedReadsAligner(genome, factory.requestClustersFinder(), factory.requestAligner());
		aligner.setCreateUnmappedReadRecords(createUnmappedReadRecords);
		aligner.setMaxAlnsPerRead(maxAlnsPerRead);
		aligner.setMinInsertLength(minInsertLength);
		aligner.setMaxInsertLength(maxInsertLength);
		List<ReadAlignment> alns = aligner.alignReads(read1, read2);
		synchronized (this) {
			totalReads+=2;
			readsAligned+=aligner.getNumReadsAligned();
			if(aligner.isProperPair()) numProperPairs+=2;
			else if (aligner.isPair()) numNonProperPairs+=2;
			else numAlignedSingle+=aligner.getNumReadsAligned();
			uniqueAlignments+=aligner.getNumUniqueAlignments();
		}
		return alns;
	}
	
	private void checkProgress (int readNumber) {
		if(readNumber%1000>0) return;
		if(!platform.isLongReads() && readNumber%100000>0) return;
		
		int progress = readNumber/100000;
		if(platform.isLongReads() ) {
			progress = readNumber/1000;
		}
		
		log.info("Processed "+readNumber+" fragments. Aligned "+readsAligned+" reads");
		if (progressNotifier!=null && !progressNotifier.keepRunning(progress)) {
			log.info("Process cancelled by user");
			pool.setCancelled(true);
		}
	}
	
	private void printStatistics(boolean paired) {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
		//TODO: Report again
		//if(shortReadsAligner!=null) {
			//out.println("Reads with less than 2 mismatches: "+shortReadsAligner.getFewMismatchesAlns());
			//out.println("Complete alignments tried: "+shortReadsAligner.getCompleteAlns());
		//}
		
		out.println("Total reads: "+totalReads);
		out.println("Reads aligned: "+readsAligned);
		if(paired) {
			out.println("    Reads aligned as proper pairs: "+numProperPairs+ " Percentage: "+fmt.format(100.0*numProperPairs/(double)totalReads)+"%");
			out.println("    Reads aligned as non proper pairs: "+numNonProperPairs+ " Percentage: "+fmt.format(100.0*numNonProperPairs/(double)totalReads)+"%");
			out.println("    Reads aligned single: "+numAlignedSingle+ " Percentage: "+fmt.format(100.0*numAlignedSingle/(double)totalReads)+"%");
		}
		out.println("Unique alignments: "+uniqueAlignments+ " Percentage: "+fmt.format(100.0*uniqueAlignments/(double)totalReads)+"%");
		out.println("Overall alignment rate: "+fmt.format(100.0*readsAligned/(double)totalReads)+"%");
		log.info(os.toString());
	}
	
}