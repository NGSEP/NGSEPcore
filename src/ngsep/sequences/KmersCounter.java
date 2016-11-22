package ngsep.sequences;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Hashtable;
import java.util.Map;

import ngsep.main.CommandsDescriptor;
import ngsep.math.Distribution;
import ngsep.sequences.io.FastaSequencesHandler;

 

public class KmersCounter {
	
	Map<String, Integer> hashKmers;
	private boolean bothStrands = false;
	private boolean fasta = false;
	private int kmerSize = 31;
	private int omitFirstBases = 0;
	QualifiedSequenceList sequenceList;
	int handicap = 0;
	Distribution kmerSpectrum;
	
	public KmersCounter(){
		
		hashKmers = new Hashtable<>();
		sequenceList = new QualifiedSequenceList();
		kmerSpectrum = new Distribution(1, 10000, 1);
		
		
	}
	/**
	 * Add the kmers of a sequence to an object that will store kmers of all sequences
	 * @param newHashKmers Kmers of a sequence.
	 */
	public void addKmersCount(Hashtable<String, Integer> newHashKmers){
		
		for (Map.Entry<String, Integer> entry : newHashKmers.entrySet()) {
		    String key = entry.getKey();
		    Integer value = entry.getValue();

		    if(hashKmers.containsKey(key)){
				hashKmers.put(key, hashKmers.get(key) + value);
			} else {
				hashKmers.put(key, value);
			}
		}
		
	}
	/**
	 * Process a fastq file to get kmers count sequence by sequence until the end of file.
	 * @param filename Fastq filename.
	 * @throws IOException
	 */
    public void processFastqSequences(String filename) throws IOException {
		 
		FileInputStream fis = null;
		BufferedReader in = null;
		try {
			fis = new FileInputStream(filename);
			//Open buffer
			in = new BufferedReader(new InputStreamReader (fis));
			//Read sequences
			while(in.ready()){
				QualifiedSequence qsequence = RawRead.load(in);
				String sequence = qsequence.getCharacters().toString();
				
				//Kmers Counter Per Sequence
				//Forward		
				Hashtable<String, Integer> newHashKmers = getSequenceKmers(sequence, getKmerSize());
				addKmersCount(newHashKmers);
				//Reverse complement
				if(isBothStrands()){
					String reverseSequence = DNAMaskedSequence.getReverseComplement(sequence);
					newHashKmers = getSequenceKmers(reverseSequence, getKmerSize());
					addKmersCount(newHashKmers);
				}
			}
		} finally {
			if(fis!=null) fis.close();
		}	
		
		
	 }
	 //fasta file reader
    /**
	 * Upload in memory the sequences of a fasta file
	 * @param filename A fasta filename.
	 * @throws IOException
	 * @return QualifiedSequenceList The list of the fasta sequence given in the file.
	 */
	 public QualifiedSequenceList loadFastaSequences(String filename) throws IOException {
		 
		FastaSequencesHandler fastaSequencesHandler = new FastaSequencesHandler();
		QualifiedSequenceList sequenceList;
		
		sequenceList = fastaSequencesHandler.loadSequences(filename);
		
		return sequenceList;	
		
	 }
	 
	 /**
	 * Read the input as fasta or fastq and get the kmers count, then make the kmers spectrum as a distribution
	 * @param sequenceFileName Filename with the sequences to process.
	 * @throws IOException
	 */
	public void processFile(String sequenceFileName) throws IOException {
		
		//Is fasta or fastq? and read it
		if(!isFasta()){
			processFastqSequences(sequenceFileName);
		} else {
			sequenceList = loadFastaSequences(sequenceFileName);
			//Kmer Count Per File
			for(int i=0; i < sequenceList.size(); i++){
				//Forward		
				String sequence = sequenceList.get(i).getCharacters().toString();
				Hashtable<String, Integer> newHashKmers = getSequenceKmers(sequence, getKmerSize());
				addKmersCount(newHashKmers);
				//Reverse complement
				if(isBothStrands()){
					String reverseSequence = DNAMaskedSequence.getReverseComplement(sequence);
					newHashKmers = getSequenceKmers(reverseSequence, getKmerSize());
					addKmersCount(newHashKmers);
				}
			}	
		}
		
		
		//Make kmer spectrum
		
		for (Map.Entry<String, Integer> entry : hashKmers.entrySet()) {
		    Integer value = entry.getValue();
		    if(value <= 10000){
		    	kmerSpectrum.processDatapoint(value);	
		    }
		    
		}
		
		double maxValue = Math.round(kmerSpectrum.getMaxValueDistribution());
		if(kmerSpectrum.getMaxValueDistribution() < 100){
			handicap = 100;
		} else {
			handicap =  (int) (maxValue + Math.round(maxValue * 1.5));
		}
		
	}
	
	/**
	 * Get the Kmers count of a sequence
	 * @param seq Genomic sequence
	 * @param length Mer length
	 * @return Hashtable<String, Integer> hashtable with the kmers count of the given sequence
	 */
	//Count kmer per sequence
	public Hashtable<String, Integer> getSequenceKmers(String seq, int length)
	{
		Hashtable<String, Integer> hashKmers = new Hashtable<>();
		int seqLength = seq.length();
		
		if(seqLength > length) {
			for(int i = 0; i < seqLength - length + 1; i++)
			{
				String kmer = seq.substring(i,length + i);
				if(hashKmers.containsKey(kmer)) {
					hashKmers.put(kmer, hashKmers.get(kmer) + 1);
				} else {
					hashKmers.put(kmer, 1);
				}
			
			}
		} else {
			System.out.println(seq);
		}
		
		return hashKmers;
	}
	 
	/**
	 * Receives the parameters from the command line interface and distributes the duties
	 * @param args
	 * @throws IOException 
	 */
	public static void main (String [ ] args) throws IOException {

		KmersCounter kmersCounter = new KmersCounter();
		
		if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
			CommandsDescriptor.getInstance().printHelp(KmersCounter.class);
			return;
		}
		//Parameters
		int k=CommandsDescriptor.getInstance().loadOptions(kmersCounter, args);
		
		String sequenceFile = args[k++];
		kmersCounter.processFile(sequenceFile);
			
		PrintStream out = System.out;
		out.println("Kmer_frequency\tNumber_of_distinct_kmers");
		kmersCounter.kmerSpectrum.printDistributionInt(out, kmersCounter.handicap);
		
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
	
	public int getOmitFirstBases() {
		return omitFirstBases;
	}
	public void setOmitFirstBases(int omitFirstBases) {
		this.omitFirstBases = omitFirstBases;
	}
	public void setOmitFirstBases(Integer omitFirstBases) {
		this.setOmitFirstBases(omitFirstBases.intValue());
	}

	 
	
}
