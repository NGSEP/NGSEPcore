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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

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
	public static final ReadAlignment.Platform DEF_PLATFORM = ReadAlignment.Platform.ILLUMINA;
	public static final byte INPUT_FORMAT_FASTQ=KmersExtractor.INPUT_FORMAT_FASTQ;
	public static final byte INPUT_FORMAT_FASTA=KmersExtractor.INPUT_FORMAT_FASTA;
	public static final int DEF_MAX_ALNS_PER_READ=3;
	public static final int DEF_KMER_LENGTH = 15;
	public static final int DEF_WINDOW_LENGTH = 5;
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
	private ReadAlignment.Platform platform = DEF_PLATFORM;
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
	private FMIndexReadAlignmentAlgorithm shortReadsAligner;
	private MinimizersTableReadAlignmentAlgorithm longReadsAligner;
	
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
		
		QualifiedSequenceList sequences = genome.getSequencesMetadata();
		if (platform.isLongReads()) {
			longReadsAligner = new MinimizersTableReadAlignmentAlgorithm();
			longReadsAligner.setLog(log);
			longReadsAligner.setMaxAlnsPerRead(maxAlnsPerRead);
			longReadsAligner.loadGenome (genome, kmerLength, windowLength);
		} else {
			if (fMIndex!=null) {
				log.info("Aligning reads using built index with "+fMIndex.getSequencesMetadata().size()+" sequences");
			} else if (fmIndexFile!=null) {
				log.info("Loading reference index from file: "+fmIndexFile);
				fMIndex = ReferenceGenomeFMIndex.load(genome, fmIndexFile);
			} else {
				log.info("Calculating FM-index from genome file: "+genome.getFilename());
				fMIndex = new ReferenceGenomeFMIndex(genome);
			}
			shortReadsAligner = new FMIndexReadAlignmentAlgorithm(fMIndex,kmerLength,maxAlnsPerRead);
		}
		
		boolean longReads = platform.isLongReads();
		boolean paired = false;
		PrintStream out = System.out;
		if(outputFile!=null) out = new PrintStream(outputFile); 
		try (ReadAlignmentFileWriter writer = new ReadAlignmentFileWriter(sequences, out)){
			// TODO: Enable threads for long reads
			if(numThreads>1 && !longReads) pool = new ThreadPoolManager(numThreads, platform.isLongReads()?1000:100000);
			writer.setSampleInfo(sampleId, platform);
			if(!longReads && inputFile!=null && inputFile2!=null) {
				log.info("Aligning paired end reads from files: "+inputFile + " and "+inputFile2);
				alignReads(inputFile,inputFile2, writer);
				paired = true;
			} else if (inputFile!=null) {
				if(longReads && inputFile2!=null) log.warning("Paired end alignment not supported for long reads. Ignoring file "+ inputFile2);
				log.info("Aligning single reads from file: "+inputFile);
				alignReads(inputFile, writer);
			} else if (inputFile2!=null ) {
				throw new IOException("The first input file is required for paired end alignment");
			} else {
				log.info("Aligning single reads from standard input");
				alignReads(System.in, writer);
			}
			if(pool!=null) {
				pool.terminatePool();
			}
			
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		log.info("Time: "+seconds+" seconds");
		printStatistics(paired);
		
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
		out.println("Platform: "+ platform);
		out.println("K-mer length: "+ kmerLength);
		if (inputFormat == INPUT_FORMAT_FASTQ)  out.println("Fastq format");
		if (inputFormat == INPUT_FORMAT_FASTA)  out.println("Fasta format");
		if (knownSTRsFile!=null) out.println("Fie with known short tandem repeats "+knownSTRsFile);
		out.println("Maximum alignments per read: "+ maxAlnsPerRead);
		if(platform.isLongReads()) {
			out.println("Window length to calculate minimizers: "+ windowLength);
		} else if (inputFile2!=null) {
			out.println("Proper limits for paired-end alignment. Minimum: "+ minInsertLength+" maximum: "+maxInsertLength);
		}
		out.println("Number of threads: "+ numThreads);
		log.info(os.toString());
	}
	
	/**
	 * Aligns readsFile with the fMIndexFile
	 * @param fMIndexFile Binary file with the serialization of an FMIndex
	 * @param readsFile Fastq file with the reads to align
	 * @param out
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void alignReads( String readsFile, ReadAlignmentFileWriter writer) throws IOException, InterruptedException {
		if(knownSTRsFile!=null && !knownSTRsFile.isEmpty()) shortReadsAligner.loadSTRsFile(knownSTRsFile);
		if(inputFormat == INPUT_FORMAT_FASTQ) {
			try (FastqFileReader reader = new FastqFileReader(readsFile)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<RawRead> it = reader.iterator();
				while(it.hasNext()) {
					RawRead read = it.next();
					if (pool == null) processSingleRead(read, writer);
					else pool.queueTask( () -> processSingleRead(read, writer));
					if (!checkProgress()) break;
				}
			}
		} else if(inputFormat== INPUT_FORMAT_FASTA) {
			try (FastaFileReader reader = new FastaFileReader(readsFile)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<QualifiedSequence> it = reader.iterator();
				while(it.hasNext()) {
					QualifiedSequence seq = it.next();
					RawRead read = new RawRead(seq.getName(), seq.getCharacters(),null);
					if (pool == null) processSingleRead(read, writer);
					else pool.queueTask( () -> processSingleRead(read, writer));
					if (!checkProgress()) break;
				}
			}
		}
	}
	private void printStatistics(boolean paired) {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		DecimalFormat fmt = ParseUtils.ENGLISHFMT;
		if(shortReadsAligner!=null) {
			out.println("Reads with less than 2 mismatches: "+shortReadsAligner.getFewMismatchesAlns());
			out.println("Complete alignments tried: "+shortReadsAligner.getCompleteAlns());
		}
		
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
	private void alignReads(InputStream in, ReadAlignmentFileWriter writer) throws IOException, InterruptedException {
		if(inputFormat == INPUT_FORMAT_FASTQ) {
			try (FastqFileReader reader = new FastqFileReader(in)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<RawRead> it = reader.iterator();
				while(it.hasNext()) {
					RawRead read = it.next();
					if (pool == null) processSingleRead(read, writer);
					else pool.queueTask( () -> processSingleRead(read, writer));
					if (!checkProgress()) break;
				}
			}
		}  else if(inputFormat== INPUT_FORMAT_FASTA) {
			try (FastaFileReader reader = new FastaFileReader(in)) {
				reader.setSequenceType(DNAMaskedSequence.class);
				Iterator<QualifiedSequence> it = reader.iterator();
				while(it.hasNext()) {
					QualifiedSequence seq = it.next();
					RawRead read = new RawRead(seq.getName(), seq.getCharacters(),null);
					if (pool == null) processSingleRead(read, writer);
					else pool.queueTask( () -> processSingleRead(read, writer));
					if (!checkProgress()) break;
				}
			}
		}
	}
	private void processSingleRead(RawRead read, ReadAlignmentFileWriter writer) {
		List<ReadAlignment> alns = alignRead(read, true);
		//System.out.println("Alignments for: "+read.getName()+" "+alns.size());
		if(alns.size()>1) {
			for(ReadAlignment aln:alns) aln.setAlignmentQuality((byte) Math.round(0.2*aln.getAlignmentQuality()/(double)alns.size()));
		}
		synchronized (writer) {
			for(ReadAlignment aln:alns) writer.write(aln);
			if(alns.size()==0) {
				ReadAlignment alnNoMap = createUnmappedAlignment(read, false, false);
				writer.write(alnNoMap);
			}
			totalReads++;
			int numAlns = alns.size();
			if(numAlns>0) readsAligned++;
			if(numAlns==1) uniqueAlignments++;
		}
	}
	
	private int lastProgress = 0;
	private boolean checkProgress () {
		int totalReadsFixed = totalReads;
		int progress = totalReadsFixed/100;
		if (!platform.isLongReads()) progress = progress/100;
		if(lastProgress==progress) return true;
		
		log.info("Processed "+totalReadsFixed+" reads. Aligned: "+readsAligned);
		if (progressNotifier!=null && !progressNotifier.keepRunning(progress)) {
			log.info("Process cancelled by user");
			return false;
		}
		lastProgress = progress;
		return true;
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
		if(knownSTRsFile!=null && !knownSTRsFile.isEmpty())shortReadsAligner.loadSTRsFile(knownSTRsFile);
		try (FastqFileReader reader1 = new FastqFileReader(readsFile1);
			 FastqFileReader reader2 = new FastqFileReader(readsFile2)) {
			reader1.setSequenceType(DNAMaskedSequence.class);
			reader2.setSequenceType(DNAMaskedSequence.class);
			Iterator<RawRead> it1 = reader1.iterator();
			Iterator<RawRead> it2 = reader2.iterator();
			while(it1.hasNext() && it2.hasNext()) {
				RawRead read1 = it1.next();
				RawRead read2 = it2.next();
				if (pool==null) processPairedEndRead(read1, read2, writer);
				else pool.queueTask(()->processPairedEndRead(read1, read2, writer));
				if(!checkProgress()) break;
			}
		}
	}
	
	private void processPairedEndRead (RawRead read1, RawRead read2, ReadAlignmentFileWriter writer) {
		List<ReadAlignment> alns1 = alignRead(read1,false);
		for(ReadAlignment aln:alns1) aln.setFirstOfPair(true);
		List<ReadAlignment> alns2 = alignRead(read2,false);
		for(ReadAlignment aln:alns2) aln.setSecondOfPair(true);
		List<ReadAlignment> alns = new ArrayList<ReadAlignment>();
		//System.out.println("Alignments found: "+alns1.size()+" "+alns2.size());
		int numMapped = 0;
		boolean proper = false;
		boolean asPair = false;
		int numUnique = 0;
		if(alns1.size()==0 || alns2.size()==0) {
			ArrayList<ReadAlignment> unMapped = processUnMapped(read1, alns1,read2,alns2);
			int n = unMapped.size();
			boolean mappedFound = false;
			for (int i = 0; i < Math.min(n,maxAlnsPerRead+1); i++) {
				ReadAlignment aln = unMapped.get(i);
				if(!aln.isReadUnmapped()) {
					mappedFound=true;
					if (n>2) aln.setAlignmentQuality((byte) Math.round(0.2*aln.getAlignmentQuality()/(double)n));
				}
				alns.add(aln);
			}
			if(mappedFound) {
				numMapped = 1;
				if(n==2) numUnique = 1;
			}
		} else {
			numMapped = 2;
			//System.out.println("Alignments read 1: "+alns1.size()+" read 2: "+alns2.size());
			List<ReadAlignmentPair> pairAlns = findPairs(alns1, alns2, true);
			//System.out.println("Pairs proper: "+pairAlns.size());
			if(pairAlns.isEmpty()) {
				pairAlns = findPairs(alns1, alns2,false);
				//System.out.println("Pairs no proper: "+pairAlns.size());
				if(pairAlns.isEmpty()) {
					int n = alns1.size();
					for(int i=0;i<n;i++) {
						ReadAlignment aln1 = alns1.get(i);
						aln1.setPaired(true);
						aln1.setMateDifferentSequence(true);
						setMateInfo(aln1, alns2.get(0));
						if(i>0) aln1.setSecondary(true);
						if (n>1) aln1.setAlignmentQuality((byte) Math.round(0.2*aln1.getAlignmentQuality()/(double)n));
						else aln1.setAlignmentQuality((byte) Math.round(0.5*aln1.getAlignmentQuality()));
						alns.add(aln1);
					}
					if(n==1) numUnique++;
					n = alns2.size();
					for(int i=0;i<n;i++) {
						ReadAlignment aln2 = alns2.get(i);
						aln2.setPaired(true);
						aln2.setMateDifferentSequence(true);
						setMateInfo(aln2, alns1.get(0));
						if(i>0) aln2.setSecondary(true);
						if (n>1) aln2.setAlignmentQuality((byte) Math.round(0.2*aln2.getAlignmentQuality()/(double)n));
						else aln2.setAlignmentQuality((byte) Math.round(0.5*aln2.getAlignmentQuality()));
						alns.add(aln2);
					}
					if(alns2.size()==1) numUnique++;
				}
				else {
					addPairAlignments(alns, pairAlns, alns1.size(),alns2.size());
					asPair=true;
					if(pairAlns.size()==1) numUnique=2;
				}

			} else {
				addPairAlignments(alns, pairAlns, alns1.size(),alns2.size());
				asPair=true;
				proper=true;
				if(pairAlns.size()==1) numUnique=2;
			}
		}
		synchronized (writer) {
			for(ReadAlignment aln:alns) writer.write(aln);
			totalReads+=2;
			readsAligned+=numMapped;
			if(proper) numProperPairs+=2;
			else if (asPair) numNonProperPairs+=2;
			else numAlignedSingle+=numMapped;
			uniqueAlignments+=numUnique;
		}
	}

	private void addPairAlignments(List<ReadAlignment> alns, List<ReadAlignmentPair> pairAlns, int numAlnsUnpaired1, int numAlnsUnpaired2) {
		int n = Math.min(pairAlns.size(),maxAlnsPerRead);
		for (int i = 0; i < n; i++) {
			ReadAlignmentPair current = pairAlns.get(i);
			ReadAlignment aln1 = current.getAln1();
			ReadAlignment aln2 = current.getAln2();
			if(i>0) {
				aln1.setSecondary(true);
				aln2.setSecondary(true);
			}
			if (n>1) {
				aln1.setAlignmentQuality((byte) Math.round(0.2*aln1.getAlignmentQuality()/(double)n));
				aln2.setAlignmentQuality((byte) Math.round(0.2*aln2.getAlignmentQuality()/(double)n));
			} else if (!aln1.isProperPair() || !aln2.isProperPair() || (numAlnsUnpaired1>1 && numAlnsUnpaired2>1)) {
				double div = Math.max(numAlnsUnpaired1+numAlnsUnpaired2-1,1);
				aln1.setAlignmentQuality((byte) Math.round(0.5*aln1.getAlignmentQuality()/div));
				aln2.setAlignmentQuality((byte) Math.round(0.5*aln2.getAlignmentQuality()/div));
			}
			alns.add(aln1);
			alns.add(aln2);
		}
	}

	public List<ReadAlignmentPair> findPairs(List<ReadAlignment> alns1, List<ReadAlignment> alns2,boolean onlyProper){
		List<ReadAlignmentPair> pairEndAlns = new ArrayList<ReadAlignmentPair>();
		int n = Math.min(alns1.size(),maxAlnsPerRead);
		for (int i = 0; i < n; i++) {
			ReadAlignment aln1 = alns1.get(i);
			//System.out.println("Pairing "+aln1.getFirst()+" is paired: "+aln1.isPaired());
			if(aln1.isPaired()) continue;
			ReadAlignmentPair alnPair = findPairForAlignment(aln1,alns2,onlyProper);
			if(alnPair!=null) {
				pairEndAlns.add(alnPair);
			}
		}
		return pairEndAlns;
	}

	public ReadAlignmentPair findPairForAlignment(ReadAlignment aln1, List<ReadAlignment> alns2,boolean onlyProper) {
		List<ReadAlignment> candidates = new ArrayList<ReadAlignment>();
		int n = Math.min(alns2.size(),maxAlnsPerRead);
		for (int i = 0; i < n; i++) {
			ReadAlignment current =alns2.get(i);
			if(!current.isPaired()) {
				//System.out.println("Next candidate "+current.getFirst()+" is paired: "+current.isPaired());
				if(isValidPair(aln1,current,onlyProper)) {
					//System.out.println("Adding candidate");
					candidates.add(current);
				}
			}
		}
		if(candidates.size()==1) {
			return buildPair(aln1, candidates.get(0), onlyProper);
		}
		else if (candidates.size()>1){
			return buildPair(aln1, pickBestPair(aln1, candidates), onlyProper);
		}
		return null;
	}

	public boolean isValidPair(ReadAlignment aln1, ReadAlignment aln2,boolean onlyProper) {
		if(!aln1.getSequenceName().equals(aln2.getSequenceName())) return false;
		int start1 = aln1.getFirst();
		int start2 = aln2.getFirst();
		int end1 = aln1.getLast();
		int end2 = aln2.getLast();

		int endMax = Math.max(end1, end2);
		int startMinimum = Math.min(start1, start2);
		if(endMax==end1 && aln1.isNegativeStrand() && startMinimum==start2&& aln2.isPositiveStrand()
				|| endMax==end2 && aln2.isNegativeStrand() && startMinimum==start1&& aln1.isPositiveStrand())
		{
			if(onlyProper) {
				int insertLength = endMax-startMinimum+1;
				return insertLength>=minInsertLength && insertLength<=maxInsertLength;
			}
			else{
				return true;
			}
		}
		return false;
	}

	private ReadAlignmentPair buildPair(ReadAlignment aln1, ReadAlignment aln2, boolean proper) {
		aln1.setPaired(true);
		
		aln2.setPaired(true);
		
		if(proper) {
			aln1.setProperPair(true);
			aln2.setProperPair(true);
		}
		int insertLength1 = aln1.getLast()-aln2.getFirst()+1;
		int insertLength2 = aln2.getLast()-aln1.getFirst()+1;
		if(insertLength1>insertLength2) {
			aln1.setInferredInsertSize(-insertLength1);
			aln2.setInferredInsertSize(insertLength1);
		} else {
			aln1.setInferredInsertSize(insertLength2);
			aln2.setInferredInsertSize(-insertLength2);
		}
		setMateInfo(aln1,aln2);
		setMateInfo(aln2,aln1);
		return new ReadAlignmentPair(aln1,aln2);
	}

	private ArrayList<ReadAlignment> processUnMapped(RawRead read1, List<ReadAlignment> alns1, RawRead read2, List<ReadAlignment> alns2) {
		ArrayList<ReadAlignment> unMappedAlignments= new ArrayList<ReadAlignment>();
		ReadAlignment unmappedAln1 = null;
		ReadAlignment unmappedAln2 = null;
		if(alns1.size()==0) unmappedAln1 = createUnmappedAlignment(read1, true, true);
		if(alns2.size()==0) unmappedAln2 = createUnmappedAlignment(read2, true, false);
		if(unmappedAln1!= null && unmappedAln2!=null) {
			unmappedAln1.setMateUnmapped(true);
			unmappedAln2.setMateUnmapped(true);
			unMappedAlignments.add(unmappedAln1);
			unMappedAlignments.add(unmappedAln2);
		}
		else if(unmappedAln2!=null) readMappedmateUnMapped(unmappedAln2, alns1, unMappedAlignments);
		else if(unmappedAln1!=null) readMappedmateUnMapped(unmappedAln1, alns2, unMappedAlignments);
		return unMappedAlignments;
	}

	private void readMappedmateUnMapped(ReadAlignment unmappedAln, List<ReadAlignment> alns, List<ReadAlignment> unMappedAlignments) {
		unMappedAlignments.add(unmappedAln);
		
		for(int i=0;i<alns.size();i++) {
			ReadAlignment aln =  alns.get(i);
			aln.setPaired(true);
			aln.setMateUnmapped(true);
			if(i==0) {
				setMateInfo(unmappedAln, aln);
			} else {
				aln.setSecondary(true);
			}
			unMappedAlignments.add(aln);
		}
		
	}

	private ReadAlignment createUnmappedAlignment(RawRead read, boolean paired, boolean first) {
		ReadAlignment alnNoMap = new ReadAlignment(null, 0, 0, read.getLength(), ReadAlignment.FLAG_READ_UNMAPPED);
		alnNoMap.setReadName(read.getName());
		alnNoMap.setReadCharacters(read.getCharacters());
		alnNoMap.setQualityScores(read.getQualityScores());
		alnNoMap.setPaired(paired);
		if(paired) {
			alnNoMap.setFirstOfPair(first);
			alnNoMap.setSecondOfPair(!first);
		}
		return alnNoMap;
	}

	private ReadAlignment setMateInfo(ReadAlignment alignment,ReadAlignment mate) {
		alignment.setMateFirst(mate.getFirst());
		alignment.setMateSequenceName(mate.getSequenceName());
		alignment.setMateNegativeStrand(mate.isNegativeStrand());
		return alignment;
	}

	private ReadAlignment pickBestPair(ReadAlignment aln1, List<ReadAlignment> alns) {
		Collections.sort(alns,new Comparator<ReadAlignment>() {

			@Override
			public int compare(ReadAlignment r1, ReadAlignment r2) { 
				int middle = (maxInsertLength+minInsertLength)/2;
				
				int d1 = Math.abs(middle - (Math.abs(aln1.getFirst()-r1.getFirst())+r1.getReadLength()));
				int d2 = Math.abs(middle - (Math.abs(aln1.getFirst()-r2.getFirst())+r2.getReadLength()));
				//System.out.println("aln1: "+aln1.getFirst()+"r1 first: "+r1.getFirst()+" r2 first: "+r2.getFirst()+" distances: "+d1+" "+d2);
				return d1-d2;
			}
		});
		return alns.get(0);
	}

	public List<ReadAlignment> alignRead(RawRead read) {
		return alignRead(read, true);
	}
	public List<ReadAlignment> alignRead(RawRead read, boolean assignSecondaryStatus) {
		List<ReadAlignment> alignments;
		if(platform.isLongReads()) {
			if(longReadsAligner ==null) longReadsAligner = new MinimizersTableReadAlignmentAlgorithm();
			alignments = longReadsAligner.alignRead(read);
		} else {
			if(shortReadsAligner==null) shortReadsAligner=new FMIndexReadAlignmentAlgorithm(fMIndex, kmerLength, maxAlnsPerRead);
			alignments = shortReadsAligner.alignRead(read);
		}
		return filterAlignments(alignments, assignSecondaryStatus);
	}
	private List<ReadAlignment> filterAlignments(List<ReadAlignment> alignments, boolean assignSecondary) {
		if (alignments.size()==0) return alignments;
		Collections.sort(alignments, (aln1,aln2) -> aln2.getAlignmentQuality() - aln1.getAlignmentQuality());
		short bestQual = alignments.get(0).getAlignmentQuality();
		//TODO. Investigate alignment score
		int threshold = (int) (0.9*bestQual);
		List<ReadAlignment> filteredAlignments = new ArrayList<>();
		int n = Math.min(alignments.size(),maxAlnsPerRead);
		for (int i=0;i<n;i++) {
			ReadAlignment aln = alignments.get(i);	
			//System.out.println("read: "+aln.getReadCharacters()+" First: "+aln.getFirst()+" flags: "+aln.getFlags()+" qual: "+aln.getAlignmentQuality()+" threshold "+threshold);
			if(aln.getAlignmentQuality()<=threshold) break;
			if(assignSecondary && i>0) aln.setSecondary(true);
			filteredAlignments.add(aln);
		}
		//System.out.println("Initial alignments: "+alignments.size()+" final: "+filteredAlignments.size());
		return filteredAlignments;
	}
	
}