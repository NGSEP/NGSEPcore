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

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.math.Distribution;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.sequences.io.FastqFileReader;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class KmersCounter {
	
	public static final int DEFAULT_KMER_SIZE = 15;
	private Logger log = Logger.getLogger(KmersCounter.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private KmersMap kmersMap = new ByteArrayKmersMapImpl((byte) DEFAULT_KMER_SIZE);
	private boolean bothStrands = false;
	private boolean fasta = false;
	private int kmerSize = DEFAULT_KMER_SIZE;
	
	
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
	 
	public boolean isBothStrands() {
		return bothStrands;
	}
	public void setBothStrands(boolean bothStrands) {
		this.bothStrands = bothStrands;
	}
	public void setBothStrands(Boolean bothStrands) {
		this.setBothStrands(bothStrands.booleanValue());
	}
	
	public boolean isFasta() {
		return fasta;
	}
	public void setFasta(boolean fasta) {
		this.fasta = fasta;
	}
	public void setFasta(Boolean fasta) {
		this.setFasta(fasta.booleanValue());
	}
	
	public int getKmerSize() {
		return kmerSize;
	}
	public void setKmerSize(int kmerSize) {
		this.kmerSize = kmerSize;
		if(kmerSize<=15) kmersMap = new ByteArrayKmersMapImpl((byte) kmerSize);
		else kmersMap = new DefaultKmersMapImpl();
	}
	public void setKmerSize(Integer kmerSize) {
		this.setKmerSize(kmerSize.intValue());
	}
	
	/**
	 * @return the hashKmers
	 */
	public KmersMap getKmersMap() {
		return kmersMap;
	}
	/**
	 * Receives the parameters from the command line interface and distributes the duties
	 * @param args
	 * @throws Exception 
	 */
	public static void main (String [ ] args) throws Exception {

		KmersCounter kmersCounter = new KmersCounter();
		//Parameters
		int k=CommandsDescriptor.getInstance().loadOptions(kmersCounter, args);
		List<String> files = new ArrayList<>();
		for(;k<args.length;k++) {
			files.add(args[k]);
		}
		
		if(files.size()>=1 && "-".equals(files.get(0))) kmersCounter.processFastqFile(System.in); 
		else kmersCounter.processFiles(files);
		kmersCounter.printResults(System.out);
		
	}
	
	/**
	 * Processes a list of input files as fasta or fastq and updates the kmers table
	 * @param files List of names of the files to process.
	 * @throws IOException If a file can not be read
	 */
	public void processFiles(List<String> files) throws IOException {
		for(String filename:files) processFile(filename);
	}
	
	/**
	 * Processes the file with the given name as fasta or fastq and updates the kmers table
	 * @param sequenceFileName Name of the file with the sequences to process.
	 * @throws IOException If the file can not be read
	 */
	public void processFile(String filename) throws IOException {
		
		//Is fasta or fastq? and read it
		if(isFasta()){
			processFastaFile(filename);
		} else {
			processFastqFile(filename);	
		}
	}
	
	/**
	 * Processes the file with the given name as fastq and updates the kmers table
	 * @param filename Name of the file with the sequences to process.
	 * @throws IOException If the file can not be read
	 */
    public void processFastqFile(String filename) throws IOException { 
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				countSequenceKmers (read);
			}
		}
	 }
    
	/**
     * Process the given input stream to get kmers count sequence by sequence
     * @param fis Input stream in fastq format
     * @throws IOException if there is an error reading the stream
     */
	public void processFastqFile(InputStream fis) throws IOException {
		try (FastqFileReader reader = new FastqFileReader(fis)) {
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				countSequenceKmers (read);
			}
		}
	}
	private void countSequenceKmers(RawRead read) {
    	String sequence = read.getCharacters().toString();
		//Kmers Counter Per Sequence
		//Forward		
		countSequenceKmers(sequence);
		//Reverse complement
		if(isBothStrands()){
			String reverseSequence = DNAMaskedSequence.getReverseComplement(sequence);
			countSequenceKmers(reverseSequence);
		}
	}
	/**
	 * Processes the file with the given name as fasta and updates the kmers table
	 * @param filename Name of the file with the sequences to process.
	 * @throws IOException If the file can not be read
	 */
    private void processFastaFile(String filename) throws IOException {
    	FastaSequencesHandler fastaSequencesHandler = new FastaSequencesHandler();
		QualifiedSequenceList sequences = fastaSequencesHandler.loadSequences(filename);
		//Kmer Count Per File
		for(QualifiedSequence seq:sequences){
			log.info("Processing sequence "+seq.getName());
			//TODO: Process in chuncks if too big
			//Forward		
			String sequence = seq.getCharacters().toString();
			countSequenceKmers(sequence);
			//Reverse complement
			if(isBothStrands()){
				String reverseSequence = DNAMaskedSequence.getReverseComplement(sequence);
				countSequenceKmers(reverseSequence);
			}
			log.info("Processed sequence "+seq.getName()+" total k-mers: "+kmersMap.size());
		}
	}
	
	/**
	 * Updates the k-mers table using the information of the given sequence
	 * @param seq CharSequence object to extract the k-mers
	 */
	public void countSequenceKmers(CharSequence seq)
	{
		
		int seqLength = seq.length();
		
		if(seqLength < kmerSize) {
			log.warning("Sequence "+seq+" smaller than k-mer size");
			return;
		}
		//TODO: Create option to process non DNA k-mers
		CharSequence [] kmers = extractKmers(seq, kmerSize, true);
		for(CharSequence kmer:kmers) {
			if(kmer!=null) kmersMap.addOcurrance(kmer);
		}	
	}
	/**
	 * Extracts the k-mers present in the given sequence
	 * @param source Sequence to process
	 * @param kmerSize Size of hte output sequences. It must be less or equal than the length of the sequence
	 * @param onlyDNA Tells if only k-mers within the DNA alphabet should be considered
	 * @return CharSequence [] Array of k-mers within the source sequence. The index in the array corresponds
	 * to the index in the sequence of the start of the k-mer
	 */
	public static CharSequence [] extractKmers(CharSequence source, int kmerSize, boolean onlyDNA) {
		return extractKmers(source, kmerSize, 0, source.length(), onlyDNA);
	}
	/**
	 * Extracts the k-mers present in the given sequence
	 * @param source Sequence to process
	 * @param kmerSize Size of the output sequences. It must be less or equal than the length of the sequence
	 * @param first position to extract k-mers. values for smaller indexes will be null
	 * @param last position to extract k-mers. The length of the array will be Math.min(last, source.length() - kmerSize) 
	 * @param onlyDNA Tells if only k-mers within the DNA alphabet should be considered
	 * @return CharSequence [] Array of k-mers within the source sequence. The index in the array corresponds
	 * to the index in the sequence of the start of the k-mer
	 */
	public static CharSequence [] extractKmers(CharSequence source, int kmerSize, int first, int last, boolean onlyDNA) {
		int n = source.length();
		if(n<kmerSize) return new CharSequence[0];
		int lastKmerStart = Math.min(last, n - kmerSize); 
		CharSequence [] kmers = new CharSequence [lastKmerStart+1];
		Arrays.fill(kmers, null);
		for(int i = Math.max(0, first); i <=lastKmerStart; i++)
		{
			String kmerStr = source.subSequence(i,kmerSize + i).toString();
			kmerStr = kmerStr.toUpperCase();
			CharSequence kmer = kmerStr;
			try {
				if(kmerSize<=31) {
					kmer = new DNAShortKmer(kmer);
				} else {
					kmer = new DNASequence(kmer);
				}
			} catch (IllegalArgumentException e) {
				if(onlyDNA) continue;
			}
			kmers[i]=kmer;
		}
		return kmers;
	}
	public void printResults (PrintStream out) {
		log.info("Calculating distribution of abundances from "+kmersMap.size()+" k-mers");
		Distribution kmerSpectrum = kmersMap.calculateAbundancesDistribution();
		out.println("Kmer_frequency\tNumber_of_distinct_kmers");
		kmerSpectrum.printDistributionInt(out);
	}	
}
