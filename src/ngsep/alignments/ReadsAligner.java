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
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.alignments.io.ReadAlignmentFileWriter;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.ReferenceGenomeFMIndex;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.main.ThreadPoolManager;
import ngsep.main.io.ParseUtils;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.LimitedSequence;
import ngsep.sequences.PairwiseAlignmentAffineGap;
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
	private int windowLength = DEF_WINDOW_LENGTH;
	private int minInsertLength = DEF_MIN_INSERT_LENGTH;
	private int maxInsertLength = DEF_MAX_INSERT_LENGTH;
	
	private int numThreads = DEF_NUM_THREADS;
	
	// Model attributes
	private Map<String, List<GenomicRegion>> knownSTRs;

	private boolean onlyPositiveStrand = false;
	
	private boolean runFullAlignment = true;
	
	private ReferenceGenome genome;

	private ReferenceGenomeFMIndex fMIndex;
	
	private Set<String> repetitiveKmers = new HashSet<String>();
	
	private LongReadsAligner longReadsAligner;
	
	private ThreadPoolManager pool;
	
	// Statistics
	private int totalReads = 0;
	private int readsAligned = 0;
	private int numProperPairs = 0;
	private int numNonProperPairs = 0;
	private int numAlignedSingle = 0;
	private int uniqueAlignments=0;
	
	private int fewMismatchesAlns = 0;
	private int completeAlns = 0;
	
	
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
	
	public ReferenceGenomeFMIndex getFMIndex() {
		return fMIndex;
	}
	public void setFMIndex(ReferenceGenomeFMIndex fMIndex) {
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
	public Map<String, List<GenomicRegion>> getKnownSTRs() {
		return knownSTRs;
	}
	public void setKnownSTRs(Map<String, List<GenomicRegion>> knownSTRs) {
		this.knownSTRs = knownSTRs;
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
		QualifiedSequenceList sequences;
		if(genome==null) {
			if (fmIndexFile!=null) {
				log.info("Loading reference index from file: "+fmIndexFile);
				fMIndex = ReferenceGenomeFMIndex.loadFromBinaries(fmIndexFile);
			} else if (fMIndex!=null) {
				log.info("Aligning reads using built index with "+fMIndex.getSequencesMetadata().size()+" sequences");
			} else {
				throw new IOException("The genome index file is a required parameter");
			}
			sequences = fMIndex.getSequencesMetadata();
		} else {
			sequences = genome.getSequencesMetadata();
			if (platform.isLongReads()) {
				longReadsAligner = new LongReadsAligner();
				longReadsAligner.setLog(log);
				longReadsAligner.loadGenome (genome, kmerLength, windowLength);
			} else {
				log.info("Calculating FM-index from genome file: "+genome.getFilename());
				fMIndex = new ReferenceGenomeFMIndex(genome);
			}
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
		else if (fmIndexFile!=null) out.println("FM index file "+fmIndexFile);
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
	public void loadSTRsFile(String strsFile) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		knownSTRs=handler.loadRegionsAsMap(strsFile);
		flat();
	}

	private void flat() {
		Set<String>keys=knownSTRs.keySet();
		Iterator<String>it=keys.iterator();
		while(it.hasNext()) {
			String key = it.next();
			List<GenomicRegion> l =knownSTRs.get(key);
			List<GenomicRegion> newList = flat(l);
			while(isOverlappging(newList)) {
				newList=flat(newList);
			}
			knownSTRs.put(key, newList);
		}
	}

	private boolean isOverlappging(List<GenomicRegion>l) {
		for (int i = 0; i < l.size()-1; i++) {
			GenomicRegion current = l.get(i);
			GenomicRegion next = l.get(i+1);
			if(isOverlappgingSorted(current, next)) {
				return true;
			}
		}
		return false;
	}

	private boolean isOverlappgingSorted(GenomicRegion a,GenomicRegion b) {
		return b.getFirst()>=a.getFirst()&&b.getFirst()<=a.getLast();
	}

	private List<GenomicRegion> flat(List<GenomicRegion> l) {
		List<GenomicRegion> newList = new ArrayList<GenomicRegion>();
		for (int i = 0; i < l.size(); i++) {
			if(i+1==l.size())newList.add(l.get(i));
			else {
				GenomicRegion current =l.get(i);
				GenomicRegion next =l.get(i+1);
				if(isOverlappgingSorted(current, next)) {
					newList.add(new GenomicRegionImpl(current.getSequenceName(), current.getFirst(), Math.max(current.getLast(), next.getLast())));
					i++;
				}
				else {
					newList.add(current);
				}
			}	
		}

		return newList;
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
		if(knownSTRsFile!=null && !knownSTRsFile.isEmpty()) loadSTRsFile(knownSTRsFile);
		if(inputFormat == INPUT_FORMAT_FASTQ) {
			try (FastqFileReader reader = new FastqFileReader(readsFile)) {
				//Load as DNAMaskedSequence to allow reverse complement
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
					//System.out.println("Aligning read "+seq.getName()+" of length: "+seq.getLength());
				
					RawRead read = new RawRead(seq.getName(), seq.getCharacters(),RawRead.generateFixedQSString('5', seq.getLength()));
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
		out.println("Reads with less than 2 mismatches: "+fewMismatchesAlns);
		out.println("Complete alignments tried: "+completeAlns);
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
				//Load as DNAMaskedSequence to allow reverse complement
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
					RawRead read = new RawRead(seq.getName(), seq.getCharacters(),RawRead.generateFixedQSString('5', seq.getLength()));
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
		if(knownSTRsFile!=null && !knownSTRsFile.isEmpty())loadSTRsFile(knownSTRsFile);
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
			List<ReadAlignmentPair> pairAlns = findPairs(alns1, alns2, true);
			if(pairAlns.isEmpty()) {
				pairAlns = findPairs(alns1, alns2,false);
				if(pairAlns.isEmpty()) {
					int n = alns1.size();
					for(int i=0;i<n;i++) {
						ReadAlignment aln1 = alns1.get(i);
						aln1.setPaired(true);
						aln1.setMateDifferentSequence(true);
						setMateInfo(aln1, alns2.get(0));
						if(i>0) aln1.setSecondary(true);
						if (n>1) aln1.setAlignmentQuality((byte) Math.round(0.2*aln1.getAlignmentQuality()/(double)n));
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
						alns.add(aln2);
					}
					if(alns2.size()==1) numUnique++;
				}
				else {
					addPairAlignments(alns, pairAlns);
					asPair=true;
					if(pairAlns.size()==1) numUnique=2;
				}

			} else {
				addPairAlignments(alns, pairAlns);
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

	private void addPairAlignments(List<ReadAlignment> alns, List<ReadAlignmentPair> pairAlns) {
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
				if(isValidPair(aln1,current,onlyProper)) {
					candidates.add(current);
				}
			}
		}
		if(candidates.size()==1) {
			return buildPair(aln1, candidates.get(0), onlyProper);
		}
		else if (candidates.size()>1){
			return buildPair(aln1, getRandomReadAlignment(candidates), onlyProper);
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

	private static ReadAlignment getRandomReadAlignment(List<ReadAlignment> alns) {
		Random r = new Random();
		return alns.get(r.nextInt(alns.size())) ;
	}

	private ReadAlignment buildAln(CharSequence query, String sequenceName, int first, int last, String cigar) {
		if(first <=0) return null;
		ReadAlignment aln = new ReadAlignment(sequenceName, first, last, query.length(), 0);
		aln.setReadCharacters(query);
		if(cigar!=null)aln.setCigarString(cigar);
		//verify last exists
		if(!fMIndex.isValidPosition(sequenceName,aln.getLast())) return null;
		return aln;
	}

	public List<ReadAlignment> alignRead(RawRead read, boolean assignSecondaryStatus) {
		List<ReadAlignment> alignments = new ArrayList<>();
		String readSeq = read.getSequenceString();
		String qual = read.getQualityScores();
		String reverseQS = null;
		if(qual == null || qual.length()!=readSeq.length()) {
			qual = RawRead.generateFixedQSString('5', readSeq.length());
			reverseQS = qual;
		} else if (!onlyPositiveStrand) {
			reverseQS = new StringBuilder(qual).reverse().toString();
		}
		String reverseComplement = null;
		if(!onlyPositiveStrand) {
			reverseComplement = DNAMaskedSequence.getReverseComplement(readSeq).toString();
			
		}
		if(fMIndex!=null && readSeq.length()<500) {
			int maxMismatches = 2;
			alignments.addAll(fewMismatchesSingleStrandSearch(readSeq,maxMismatches));
			//System.out.println("Read: "+read.getName()+" Forward exact alignments: "+alignments.size());
			if(reverseComplement!=null) {
				List<ReadAlignment> alnsR = fewMismatchesSingleStrandSearch(reverseComplement,maxMismatches);
				//System.out.println("Read: "+read.getName()+" Reverse exact alignments: "+alnsR.size());
				for (ReadAlignment aln:alnsR) aln.setNegativeStrand(true);
				alignments.addAll(alnsR);
			}
		}
		if(alignments.size()==0) {
			alignments.addAll(inexactSearchAlgorithm(readSeq));
			//System.out.println("Read: "+read.getName()+" Forward inexact alignments: "+alignments.size());
			if(reverseComplement!=null) {
				List<ReadAlignment> alnsR = inexactSearchAlgorithm(reverseComplement);
				//System.out.println("Read: "+read.getName()+" Reverse inexact alignments: "+alnsR.size());
				for (ReadAlignment aln:alnsR) aln.setNegativeStrand(true);
				alignments.addAll(alnsR);
			}
		} else fewMismatchesAlns++;
		
		//System.out.println("Read: "+read.getName()+" total alignments: "+alignments.size());
		for(ReadAlignment aln:alignments) {
			aln.setReadName(read.getName());
			if(!aln.isNegativeStrand()) aln.setQualityScores(qual);
			else aln.setQualityScores(reverseQS);
		}
		
		return filterAlignments(alignments, assignSecondaryStatus);
	}
	private List<ReadAlignment> fewMismatchesSingleStrandSearch(String query, int maxMismatches) {
		List<ReadAlignment> alns = new ArrayList<ReadAlignment>();
		String cigar = ""+query.length();
		cigar += ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		// Whole read exact search
		List<UngappedSearchHit> readHits=fMIndex.exactSearch(query);
		for(int i=0;i<readHits.size()&& i<maxAlnsPerRead;i++) {
			UngappedSearchHit hit = readHits.get(i);
			ReadAlignment aln = buildAln (query, hit.getSequenceName(), hit.getStart()+1, hit.getStart()+query.length(), cigar);
			aln.setAlignmentQuality((byte) 100);
			if(aln!=null) alns.add(aln);
		}
		if (alns.size()>0) return alns;
		//One mismatch search
		int middle = query.length()/2;
		if(middle < 50) return alns;
		String firstPart = query.substring(0,middle);
		readHits=fMIndex.exactSearch(firstPart);
		for(UngappedSearchHit hit: readHits) {
			ReadAlignment aln = buildAln (query, hit.getSequenceName(), hit.getStart()+1, hit.getStart()+query.length(), cigar);
			if (aln==null) continue;
			int[] mismatches = countMismatches(query, aln);
			if(mismatches==null) continue;
			if(mismatches[0]<=maxMismatches) {
				if (mismatches[1]+mismatches[2]>0) {
					aln = buildAln(query, hit.getSequenceName(), hit.getStart()+1+mismatches[1], hit.getStart()+query.length()-mismatches[2], makeCigar(query.length(),mismatches));
				}
				aln.setAlignmentQuality((byte) (100-5*mismatches[0]));
				aln.setNumMismatches((short) mismatches[0]);
				
				alns.add(aln);
				//Best alignments selected during the filtering step
				if(alns.size()>=3*maxAlnsPerRead) break;
			}
		}
		if (alns.size()>0) return alns;
		String secondPart = query.substring(middle);
		readHits=fMIndex.exactSearch(secondPart);
		for(UngappedSearchHit hit: readHits) {
			int start = hit.getStart()-middle;
			if(start<0) continue;
			ReadAlignment aln = buildAln (query, hit.getSequenceName(), start+1, start+query.length(), cigar);
			if (aln==null) continue;
			int[] mismatches = countMismatches(query, aln);
			if(mismatches==null) continue;
			if(mismatches[0]<=maxMismatches) {
				if (mismatches[1]+mismatches[2]>0) {
					aln = buildAln(query, hit.getSequenceName(), hit.getStart()+1+mismatches[1], hit.getStart()+query.length()-mismatches[2], makeCigar(query.length(),mismatches));
				}
				aln.setAlignmentQuality((byte) (100-5*mismatches[0]));
				aln.setNumMismatches((short) mismatches[0]);
				alns.add(aln);
				//Best alignments selected during the filtering step
				if(alns.size()>=3*maxAlnsPerRead) break;
			}
		}
		return alns;
	}
	private String makeCigar(int length, int[] mismatches) {
		String cigar = "";
		int l2 = length-mismatches[1]-mismatches[2];
		if(mismatches[1]>0) {
			cigar+=mismatches[1];
			cigar += ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_SKIPFROMREAD);
		}
		cigar+=l2;
		cigar += ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		if(mismatches[2]>0) {
			cigar+=mismatches[2];
			cigar += ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_SKIPFROMREAD);
		}
		return cigar;
	}
	private List<ReadAlignment> inexactSearchAlgorithm(String readSeq) {
		List<ReadAlignment> alignments;
		if(platform.isLongReads()) {
			alignments = longReadsAligner.alignQueryToReference(readSeq);
		} else {
			alignments = kmerBasedSingleStrandInexactSearchAlgorithm(readSeq);
		}
		return alignments;
	}
	/**
	 * Inexact search of kmers to an FM-index
	 * It iterates the Sequences of the genome and if there is at least MIN_ACCURACY percentage of the kmers
	 * it allow the alignment with the first an the last position of the kmers ocurrence
	 * Only tries to align the given quey in the positive strand
	 * @return List<ReadAlignment>
	 */
	private List<ReadAlignment> kmerBasedSingleStrandInexactSearchAlgorithm (String query) 
	{
		Map<Integer,CharSequence> kmersMap = KmersExtractor.extractKmersAsMap(query, kmerLength, kmerLength, true, true, true);
		List<ReadAlignment> finalAlignments =  new ArrayList<>();
		//System.out.println("Read: "+query+" length "+query.length()+" kmers: "+kmersMap.size());
		int kmersCount=kmersMap.size();
		if(kmersCount==0) return finalAlignments;
		List<UngappedSearchHit> initialKmerHits = searchKmers (kmersMap);
		List<KmerHitsCluster> clusteredKmerHits = clusterKmerHits(query, initialKmerHits);
		if(clusteredKmerHits.size()==0) return finalAlignments;
		//System.out.println("Initial kmer hits: "+initialKmerHits.size()+" Clusters: "+clusteredKmerHits.size());
		Collections.sort(clusteredKmerHits, (o1, o2) -> o2.getNumDifferentKmers()-o1.getNumDifferentKmers());
		//KmerHitsCluster cluster = clusteredKmerHits.get(0);
		//ReadAlignment readAln = createNewAlignmentFromConsistentKmers(cluster, query);
		//if(readAln!=null) finalAlignments.add(readAln);
		int kmersMaxCluster = 0;
		for (int i=0;i<clusteredKmerHits.size() && i<2*maxAlnsPerRead;i++) {
			KmerHitsCluster cluster = clusteredKmerHits.get(i);
			int numKmers = cluster.getNumDifferentKmers();
			//System.out.println("Processing cluster "+i+" spanning "+cluster.getSequenceName()+":"+cluster.getSubjectPredictedStart()+"-"+cluster.getSubjectPredictedEnd()+" Num kmers: "+cluster.getNumDifferentKmers()+" consistent: "+cluster.isAllConsistent());
			if(i==0) kmersMaxCluster = numKmers;
			else if (finalAlignments.size()>0 && (numKmers<2 || numKmers< 0.5*kmersMaxCluster)) break;
			ReadAlignment readAln = createNewAlignmentFromConsistentKmers(cluster, query);
			if(readAln!=null) finalAlignments.add(readAln);
		}
		//System.out.println("Found "+finalAlignments.size()+" alignments for query: "+query);
		return finalAlignments;
	}


	/**
	 * Searches the given kmers in the fmIndex 
	 * @param kmers to search
	 * @return List of alignments of each kmer. The read number of each alignment contains the kmer number.
	 */
	private List<UngappedSearchHit> searchKmers(Map<Integer,CharSequence> kmersMap) {
		List<UngappedSearchHit> answer = new ArrayList<>();
		for (int start:kmersMap.keySet()) {
			String kmer = kmersMap.get(start).toString();
			if(repetitiveKmers.contains(kmer)) continue;
			List<UngappedSearchHit> kmerHits=fMIndex.exactSearch(kmer);
			//System.out.println("Kmer: "+kmer+" hits: "+kmerHits.size());
			if(kmerHits.size()>50) {
				repetitiveKmers.add(kmer);
				continue;
			}
			for(UngappedSearchHit hit:kmerHits) {
				hit.setQueryIdx(start);
				answer.add(hit);
			}
		}
		return answer;
	}

	private List<KmerHitsCluster> clusterKmerHits(String query, List<UngappedSearchHit> initialKmerHits) {
		List<KmerHitsCluster> clusters = new ArrayList<>();
		Map<String,List<UngappedSearchHit>> hitsBySubjectName = new LinkedHashMap<String, List<UngappedSearchHit>>();
		for(UngappedSearchHit hit:initialKmerHits) {
			List<UngappedSearchHit> hitsSeq = hitsBySubjectName.computeIfAbsent(hit.getSequenceName(), k -> new ArrayList<>());
			hitsSeq.add(hit);
		}
		for(List<UngappedSearchHit> hitsSeq:hitsBySubjectName.values()) {
			Collections.sort(hitsSeq, (hit0,hit1)-> hit0.getStart()-hit1.getStart());
			clusters.addAll(clusterSequenceKmerAlns(query, hitsSeq));
		}
		return clusters;
	}
	
	private List<KmerHitsCluster> clusterSequenceKmerAlns(String query, List<UngappedSearchHit> sequenceHits) {
		List<KmerHitsCluster> answer = new ArrayList<>();
		//System.out.println("Alns to cluster: "+sequenceAlns.size());
		KmerHitsCluster cluster=null;
		for(UngappedSearchHit kmerHit:sequenceHits) {
			if(cluster==null || !cluster.addKmerHit(kmerHit, 0)) {
				cluster = new KmerHitsCluster(query, kmerHit);
				answer.add(cluster);
			}
		}
		return answer;
	}

	private ReadAlignment createNewAlignmentFromConsistentKmers(KmerHitsCluster cluster, String query) {
		String sequenceName = cluster.getSequenceName();
		int first = cluster.getSubjectPredictedStart()+1;
		int last = cluster.getSubjectPredictedEnd();
		int lastPerfect = first+query.length()-1;
		//System.out.println("Reference region from k-mers: "+first+"-"+last+" all consistent: "+cluster.isAllConsistent()+" lastPerfect: "+lastPerfect+" firstaln: "+cluster.isFirstKmerPresent()+" last aln: "+cluster.isLastKmerPresent());
		String cigar = ""+query.length()+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		ReadAlignment aln = buildAln(query, sequenceName, first, lastPerfect, cigar);
		if(aln!=null) {
			//System.out.println("Built alignment at "+sequenceName+":"+aln.getFirst()+"-"+aln.getLast()+" CIGAR: "+cigar);
			GenomicRegion region =findTandemRepeat(sequenceName,first,last);
			if(region!=null) {
				ReadAlignment newaln=verifyShortTandemRepeats(aln.getSequenceName(),aln.getFirst(), aln.getLast(),query,region);
				//System.out.println("Found overlapping tandem repeat at "+region.getSequenceName()+":"+region.getFirst()+"-"+region.getLast()+" new aln: "+newaln);
				if(newaln!=null) return newaln;
			}
			if(cluster.getNumDifferentKmers()>2 && cluster.isAllConsistent()) {
				int [] mismatches = countMismatches (query, aln);
				if (mismatches[1]+mismatches[2]>0) {
					aln = buildAln(query, sequenceName, first+mismatches[1], lastPerfect-mismatches[2], makeCigar(query.length(),mismatches));
				}
				aln.setAlignmentQuality((byte) Math.round(100-5*mismatches[0]));
				aln.setNumMismatches((short) mismatches[0]);
				//System.out.println("Mismatches alignment at "+aln.getSequenceName()+":"+aln.getFirst()+"-"+aln.getLast()+": "+mismatches);
				if(mismatches[0]<0.05*query.length()) return aln;
			}
		}
		if(!runFullAlignment) return null;
		//Perform smith waterman
		first = Math.max(1, first);
		last = Math.min(fMIndex.getReferenceLength(sequenceName), last);
		CharSequence refSeq = fMIndex.getSequence(sequenceName, first, last);
		if(refSeq == null) return null;
		if(refSeq.length()>1.5*query.length()) return null;
		//System.out.println("Aligning reference from "+first+" to "+last+ " to query. length: "+refSeq.length());
		completeAlns++;
		PairwiseAlignmentWithCigar pAln = new PairwiseAlignmentWithCigar(query, refSeq.toString(), false);
		//System.out.println("Pairwise alignment found from relative : "+pAln.getSubjectStartIdx()+" to "+pAln.getSubjectLastIdx()+" CIGAR:" +pAln.getCigar()+" ");
		short mismatches = (short) pAln.getMismatches();
		if(mismatches>0.1*query.length()) return null;
		//Last must be updated before first
		last = first + pAln.getSubjectLastIdx();
		first = first + pAln.getSubjectStartIdx();
		cigar = pAln.getCigar();
		//System.out.println("New genomic coordinates : "+first+"-"+last+" CIGAR:" +cigar);
		aln = buildAln(query, sequenceName, first, last, cigar);
		
		aln.setAlignmentQuality((byte) Math.round(100-5*mismatches));
		aln.setNumMismatches(mismatches);
		return aln;
	}

	private int [] countMismatches(CharSequence query, ReadAlignment aln) {
		int [] answer = {0,0,0};
		CharSequence refS = fMIndex.getSequence(aln.getSequenceName(), aln.getFirst(), aln.getLast());
		if(refS==null) return null;
		String refSeq = refS.toString();
		int lastMismatch = -1;
		for (int i=0;i<query.length() && i<refSeq.length();i++ ) {
			if(query.charAt(i)!=refSeq.charAt(i)) {
				answer[0]++;
				lastMismatch=i;
			} else if (answer[0]+3<i) answer[1]=lastMismatch+1;
		}
		if (query.length()!=refSeq.length()) {
			answer[0]+=Math.abs(query.length()-refSeq.length());
			answer[2]=Math.max(0, query.length()-refSeq.length());
		} else {
			lastMismatch=refSeq.length();
			int numM =0;
			for (int i=query.length()-1;i>=0;i-- ) {
				if (query.charAt(i)!=refSeq.charAt(i)) {
					lastMismatch = i;
					numM++;
				} else {
					int revIdx = refSeq.length()-1-i;
					if (numM+3<revIdx) {
						answer[2]=refSeq.length()-lastMismatch;
						break;
					}
				}
			}
		}
		
		return answer;
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
			if(aln.getAlignmentQuality()<=threshold) break;
			//System.out.println("read: "+aln.getReadCharacters()+" First: "+aln.getFirst()+" flags: "+aln.getFlags()+" qual: "+aln.getAlignmentQuality()+" threshold "+threshold);
			if(assignSecondary && i>0) aln.setSecondary(true);
			filteredAlignments.add(aln);
		}
		return filteredAlignments;
	}
	
	public GenomicRegion findTandemRepeat(String sequenceName, int first, int last) {
		GenomicRegionImpl region = new GenomicRegionImpl(sequenceName, first, last);
		if(knownSTRs==null) return null;
		List<GenomicRegion> l =knownSTRs.get(sequenceName);
		if(l==null) return null;
		return binaryContains(l, 0, l.size()-1, region);
	}

	private GenomicRegion binaryContains(List<GenomicRegion> l, int left, int rigth, GenomicRegion element) {
		if(rigth>=left) {
			int middle = left +(rigth-left)/2;
			GenomicRegion actual = l.get(middle);
			if(isOverlappgingSorted(actual, element)||isOverlappgingSorted(element, actual	)) {
				return actual;
			}
			if(l.get(middle).getFirst()>element.getFirst()) {
				return binaryContains(l, left, middle-1, element);
			}
			return binaryContains(l, middle+1,rigth, element);
		}
		return null;
	}
	
	/**
	 * Creates an alignment taking into account that the region overlaps with a tandem repeat
	 * @param aln
	 * @param read
	 * @param qualityScores
	 * @param region
	 * @return
	 */
	public ReadAlignment verifyShortTandemRepeats(String sequenceName, int first, int last, String read, GenomicRegion region) {
		PairwiseAlignmentWithCigar leftPart = null;
		PairwiseAlignmentWithCigar rightPart = null;
		int firstLeftPart = Math.max(first,1);
		int softClipLeft = 0;
		int softClipRight = 0;
		String cigarLeft=null;
		String cigarRight=null;
		char softClipChar = ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_SKIPFROMREAD);
		if(first<region.getFirst()-5) {
			CharSequence refSeq = fMIndex.getSequence(sequenceName, firstLeftPart, region.getFirst()-1).toString();
			if(refSeq!=null) {
				int endReadSegment= Math.min(read.length(), region.getFirst()-first+5);
				String readSegment = read.substring(0,endReadSegment);
				//System.out.println(refSeq);
				//System.out.println(readSegment);
				leftPart = new PairwiseAlignmentWithCigar(readSegment, refSeq.toString(), false);
				softClipLeft = readSegment.length()-1 - leftPart.getQueryLastIdx();
				String softClipSubstr = "";
				if(softClipLeft>0) softClipSubstr = ""+softClipLeft+""+softClipChar;
				//System.out.println("Left part query start: "+leftPart.getQueryStartIdx()+" subject start: "+leftPart.getSubjectStartIdx()+" mismatches: "+leftPart.getMismatches()+" cigar: "+leftPart.getCigar()+" soft clip: "+softClipLeft);
				if(leftPart.getMismatches()>3 || firstLeftPart+leftPart.getSubjectLastIdx()!=region.getFirst()-1 || !leftPart.getCigar().endsWith(softClipSubstr)) leftPart = null;
				else {
					//Remove soft clip part and update soft clip total with removed segment of the read
					cigarLeft = leftPart.getCigar().substring(0,leftPart.getCigar().length()-softClipSubstr.length());
					softClipLeft+=(read.length()-readSegment.length());
				}
			}	
		}
		if(last>region.getLast()+5) {
			CharSequence refSeq = fMIndex.getSequence(sequenceName, region.getLast()+1, last);
			if(refSeq!=null) {
				int startReadSegment= Math.max(0, read.length()-(last-region.getLast())-5);
				String readSegment = read.substring(startReadSegment);
				//System.out.println(refSeq);
				//System.out.println(readSegment);
				rightPart = new PairwiseAlignmentWithCigar(readSegment, refSeq.toString(), false);
				softClipRight = rightPart.getQueryStartIdx();
				String softClipSubstr = "";
				if(softClipRight>0) softClipSubstr = ""+softClipRight+""+softClipChar;
				//System.out.println("Right part query start: "+rightPart.getQueryStartIdx()+" subject start: "+rightPart.getSubjectStartIdx()+" mismatches: "+rightPart.getMismatches()+" cigar: "+rightPart.getCigar()+" soft clip: "+softClipRight);
				if(rightPart.getMismatches()>3 || rightPart.getSubjectStartIdx()!=0 || !rightPart.getCigar().startsWith(softClipSubstr)) rightPart = null;
				else {
					cigarRight = rightPart.getCigar().substring(softClipSubstr.length());
					softClipRight+=startReadSegment;
				}
			}	
		}
		if(leftPart==null && rightPart ==null) {
			return null;
		}
		if(rightPart==null) {
			//Left alignment with right soft clip
			first = firstLeftPart + leftPart.getSubjectStartIdx();
			last = region.getFirst()-1;
			String cigar = cigarLeft;
			if(softClipLeft>0) cigar+=""+softClipLeft+""+softClipChar;
			short mismatches = (short) leftPart.getMismatches();
			//System.out.println("Left alignment new genomic coordinates : "+first+"-"+last+" CIGAR:" +cigar+" quality: "+alnQual);
			ReadAlignment aln = buildAln(read, sequenceName, first, last, cigar);
			aln.setAlignmentQuality((byte)(90-5*mismatches));
			aln.setNumMismatches(mismatches);
			return aln;
		}
		if(leftPart==null) {
			//Right alignment with left soft clip
			first = region.getLast()+1;
			last = first + rightPart.getSubjectLastIdx();
			String cigar = "";
			if(softClipRight>0) cigar+=softClipRight+""+softClipChar;
			cigar+=cigarRight;
			short mismatches = (short) rightPart.getMismatches();
			//System.out.println("Right alignment new genomic coordinates : "+first+"-"+last+" CIGAR:" +cigar+" quality: "+alnQual);
			ReadAlignment aln = buildAln(read, sequenceName, first, last, cigar);
			aln.setAlignmentQuality((byte)(90-5*mismatches));
			aln.setNumMismatches(mismatches);
			return aln;
		}
		first = firstLeftPart + leftPart.getSubjectStartIdx();
		last = region.getLast()+1 + rightPart.getSubjectLastIdx();
		int alignedLeft = read.length()-softClipLeft;
		int alignedRight = read.length()-softClipRight;
		int middleLength = read.length()-alignedLeft-alignedRight;
		if(middleLength<0) return null;
		int difference = region.length()-middleLength;
		//System.out.println("Aligned left: "+alignedLeft+" aligned right: "+alignedRight+" middle length: "+middleLength+" difference: "+difference);
		String cigar = cigarLeft;
		if(difference>0) {
			// Region length > middle length. Add deletion
			cigar+=(difference)+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_DELETION);
			if(middleLength>0) cigar+= (middleLength)+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		} else if (difference<0) {
			// Region length < middle length. Add insertion
			cigar+=(-difference)+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_INSERTION);
			if(middleLength>0) cigar+= (region.length())+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		} else {
			if(middleLength>0) cigar+= (middleLength)+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_MATCH);
		}
		
		cigar+=cigarRight;
		short mismatches = (short) (leftPart.getMismatches()+rightPart.getMismatches());
		//System.out.println("Building alignment from first "+first+" last: "+last+" softClipLeft: "+softClipLeft+" softClip right "+softClipRight+" cigar "+cigar);
		ReadAlignment aln = buildAln(read, sequenceName, first, last, cigar);
		aln.setAlignmentQuality((byte)(100-5*mismatches));
		aln.setNumMismatches(mismatches);
		return aln;
	}
}

class PairwiseAlignmentWithCigar {
	private static PairwiseAlignmentAffineGap aligner = new PairwiseAlignmentAffineGap(1, 2, 1, 1, 1000);
	private String query;
	private String subject;
	private int subjectStartIdx=-1;
	private int subjectLastIdx=-1;
	private int queryStartIdx=-1;
	private int queryLastIdx=-1;
	private String cigar;
	private int mismatches=0;
	
	public PairwiseAlignmentWithCigar (String query, String subject, boolean includeEndSubject) {
		//System.out.println("Query length: "+query.length()+" subject length: "+subject.length());
		this.query = query;
		this.subject = subject;
		String [] alignment;
		synchronized (aligner) {
			alignment = aligner.getAlignment(query, subject);
		}
		 
		String alnQuery = alignment[0];
		//System.out.println(alnQuery);
		String alnSubject = alignment[1];
		//System.out.println(alnSubject);
		StringBuilder cigarBuilder = new StringBuilder();
		int subjectIdx = 0;
		int queryIdx = 0;
		byte nextOpCode = -1;
		int nextOpLength = 0;
		int lastMismatches = 0;
		boolean addDeletions = includeEndSubject;
		for(int i=0;i<alnQuery.length();i++) {
			char q = alnQuery.charAt(i);
			char s = alnSubject.charAt(i);
			boolean match = q==s;
			byte opCode = ReadAlignment.ALIGNMENT_MATCH;
			if(q!=LimitedSequence.GAP_CHARACTER) {
				if(s!=LimitedSequence.GAP_CHARACTER) {
					if(subjectStartIdx==-1) {
						if (match) {
							subjectStartIdx = subjectIdx;
							queryStartIdx = queryIdx;
						}
						else opCode = ReadAlignment.ALIGNMENT_SKIPFROMREAD;
					}
					if (queryStartIdx>=0) {
						subjectLastIdx = subjectIdx;
						queryLastIdx = queryIdx;
					}
					if(!match) {
						mismatches++;
						lastMismatches++;
					} else lastMismatches=0;
					subjectIdx++;
				} else {
					if(subjectStartIdx<0 || subjectIdx>=subject.length()) {
						opCode = ReadAlignment.ALIGNMENT_SKIPFROMREAD;
					} else {
						opCode = ReadAlignment.ALIGNMENT_INSERTION;
					}
				}
				queryIdx++;
			} else if (s!=LimitedSequence.GAP_CHARACTER) {
				opCode = ReadAlignment.ALIGNMENT_DELETION;
				subjectIdx++;
			}
			if(opCode!=nextOpCode) {
				if( nextOpCode!=ReadAlignment.ALIGNMENT_DELETION || addDeletions) {
					int lastSkip = 0;
					if(nextOpCode == ReadAlignment.ALIGNMENT_MATCH && lastMismatches>0 && (subjectIdx>=subject.length()||queryIdx>=query.length())) {
						//TODO: Improve end management
						lastSkip=Math.min(lastMismatches, nextOpLength-1);
						nextOpLength-=lastSkip;
						subjectLastIdx-=lastSkip;
						queryLastIdx-=lastSkip;
					}
					if(nextOpLength>0) cigarBuilder.append(nextOpLength+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(nextOpCode));
					if (lastSkip>0) cigarBuilder.append(lastSkip+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_SKIPFROMREAD));
					if(nextOpCode>=0 && nextOpCode!=ReadAlignment.ALIGNMENT_MATCH && nextOpCode!=ReadAlignment.ALIGNMENT_SKIPFROMREAD) mismatches+=2;
					if(nextOpLength>0 && nextOpCode!=ReadAlignment.ALIGNMENT_SKIPFROMREAD) addDeletions=true;
				}
				nextOpCode = opCode;
				nextOpLength = 0;
				lastMismatches = 0;
			}
			nextOpLength++;
		}
		if( nextOpCode!=ReadAlignment.ALIGNMENT_DELETION || includeEndSubject) {
			int lastSkip = 0;
			if(nextOpCode == ReadAlignment.ALIGNMENT_MATCH && lastMismatches>0) {
				lastSkip=Math.min(lastMismatches, nextOpLength-1);
				nextOpLength-=lastSkip;
				subjectLastIdx-=lastSkip;
				queryLastIdx-=lastSkip;
			}
			if(nextOpLength>0) cigarBuilder.append(nextOpLength+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(nextOpCode));
			if (lastSkip>0) cigarBuilder.append(lastSkip+""+ReadAlignment.ALIGNMENT_CHAR_CODES.charAt(ReadAlignment.ALIGNMENT_SKIPFROMREAD));
			if(nextOpCode>=0 && nextOpCode!=ReadAlignment.ALIGNMENT_MATCH && nextOpCode!=ReadAlignment.ALIGNMENT_SKIPFROMREAD) mismatches+=2;
		}
		cigar = cigarBuilder.toString();
	}

	public String getQuery() {
		return query;
	}

	public String getSubject() {
		return subject;
	}


	public int getSubjectStartIdx() {
		return subjectStartIdx;
	}
	
	public int getSubjectLastIdx() {
		return subjectLastIdx;
	}
	
	public int getQueryStartIdx() {
		return queryStartIdx;
	}

	public int getQueryLastIdx() {
		return queryLastIdx;
	}

	public String getCigar() {
		return cigar;
	}

	public int getMismatches() {
		return mismatches;
	}
}