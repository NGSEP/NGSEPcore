package ngsep.sequences;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Hashtable;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.math.Distribution;
import ngsep.sequences.io.FastaSequencesHandler;

 

public class KmersCounter {
	
	private Logger log = Logger.getLogger(KmersCounter.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private Map<String, Integer> hashKmers = new Hashtable<>();
	private boolean bothStrands = false;
	private boolean fasta = false;
	private int kmerSize = 21;
	private Distribution kmerSpectrum = new Distribution(1, 200, 1);
	
	
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
	
	
	/**
	 * Process a fastq file to get kmers count sequence by sequence until the end of file.
	 * @param filename Fastq filename.
	 * @throws IOException
	 */
    public void processFastqFile(String filename) throws IOException {
		 
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(filename);
			if(filename.toLowerCase().endsWith(".gz")) {
				//log.info("Processing compressed fastq file: "+filename);
				processFastqFile (new ConcatGZIPInputStream(fis));
			} else {
				processFastqFile(fis);
			}
		} finally {
			if(fis!=null) fis.close();
		}
	 }
    /**
     * Process the given input stream to get kmers count sequence by sequence
     * @param fis Input stream in fastq format
     * @throws IOException
     */
	public void processFastqFile(InputStream fis) throws IOException {
		BufferedReader in = new BufferedReader(new InputStreamReader (fis));
		//Read sequences
		QualifiedSequence qsequence = RawRead.load(in);
		while(qsequence!=null){
			
			String sequence = qsequence.getCharacters().toString();
			
			//Kmers Counter Per Sequence
			//Forward		
			countSequenceKmers(sequence);
			//Reverse complement
			if(isBothStrands()){
				String reverseSequence = DNAMaskedSequence.getReverseComplement(sequence);
				countSequenceKmers(reverseSequence);
			}
			qsequence = RawRead.load(in);
		}
	}
    /**
	 * Processes a fasta file to get kmers count sequence by sequence until the end of file.
	 * @param filename A fasta filename.
	 * @throws IOException
	 */
    private void processFastaFile(String filename) throws IOException {
    	FastaSequencesHandler fastaSequencesHandler = new FastaSequencesHandler();
		QualifiedSequenceList sequences = fastaSequencesHandler.loadSequences(filename);
		//Kmer Count Per File
		for(QualifiedSequence seq:sequences){
			//Forward		
			String sequence = seq.getCharacters().toString();
			countSequenceKmers(sequence);
			//Reverse complement
			if(isBothStrands()){
				String reverseSequence = DNAMaskedSequence.getReverseComplement(sequence);
				countSequenceKmers(reverseSequence);
			}
		}
	}
	 
	 /**
	 * Read the input as fasta or fastq and get the kmers count, then make the kmers spectrum as a distribution
	 * @param sequenceFileName Filename with the sequences to process.
	 * @throws IOException
	 */
	public void processFile(String sequenceFileName) throws IOException {
		
		//Is fasta or fastq? and read it
		if(isFasta()){
			processFastaFile(sequenceFileName);
		} else {
			processFastqFile(sequenceFileName);	
		}
		
		
		//Make kmer spectrum
		
		for (Map.Entry<String, Integer> entry : hashKmers.entrySet()) {
		    Integer value = entry.getValue();
		    kmerSpectrum.processDatapoint(value);
		}
	}
	
	public void printResults (PrintStream out) {
		out.println("Kmer_frequency\tNumber_of_distinct_kmers");
		kmerSpectrum.printDistributionInt(out);
		out.println("More:\t"+kmerSpectrum.getOutliers().size());
	}
	
	/**
	 * Get the Kmers count of a sequence
	 * @param seq Genomic sequence
	 * @param length Mer length
	 * @return Hashtable<String, Integer> hashtable with the kmers count of the given sequence
	 */
	//Count kmer per sequence
	public void countSequenceKmers(String seq)
	{
		int seqLength = seq.length();
		
		if(seqLength > kmerSize) {
			for(int i = 0; i < seqLength - kmerSize + 1; i++)
			{
				String kmer = seq.substring(i,kmerSize + i);
				if(hashKmers.containsKey(kmer)) {
					hashKmers.put(kmer, hashKmers.get(kmer) + 1);
				} else {
					hashKmers.put(kmer, 1);
				}
			
			}
		} else {
			log.warning("Sequence "+seq+" smaller than k-mer size");
		}
		
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
		
		String sequenceFile = args[k++];
		kmersCounter.processFile(sequenceFile);
		kmersCounter.printResults(System.out);
		
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
	}
	public void setKmerSize(Integer kmerSize) {
		this.setKmerSize(kmerSize.intValue());
	}	
}
