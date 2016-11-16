package ngsep.kmerscount;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Hashtable;
import java.util.Map;

import ngsep.main.CommandsDescriptor;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaSequencesHandler;

 

public class KmersCount {
	
	Hashtable<String, Integer> hashKmers;
	Hashtable<Integer, Integer> kmerDistribution;
	private boolean bothStrands = false;
	private boolean fasta = false;
	private int kmerSize = 31;
	private int omitFirstBases = 0;
	QualifiedSequenceList sequenceList;
	
	public KmersCount(){
		
		hashKmers = new Hashtable<>();
		kmerDistribution = new Hashtable<>();
		sequenceList = new QualifiedSequenceList();
		
		
	}
	//add kmer file or strand count
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
	//fastq file reader
    public static QualifiedSequenceList loadFastqSequences(String filename) throws IOException {
		 
		QualifiedSequenceList answer = new QualifiedSequenceList();
		FileInputStream fis = null;
		BufferedReader in = null;
		try {
			fis = new FileInputStream(filename);
			//Open buffer
			in = new BufferedReader(new InputStreamReader (fis));
			//Read sequences
			while(in.ready()){
				QualifiedSequence sequences = RawRead.load(in);
				answer.add(sequences);
			}
		} finally {
			if(fis!=null) fis.close();
		}	
		return answer;
		
	 }
	 //fasta file reader
	 public static QualifiedSequenceList loadFastaSequences(String filename) throws IOException {
		 
		FastaSequencesHandler fastaSequencesHandler = new FastaSequencesHandler();
		QualifiedSequenceList sequenceList;
		
		sequenceList = fastaSequencesHandler.loadSequences(filename);
		
		return sequenceList;	
		
	 }
	//read the fasta or fastq file using ngsep native libraries
	public void processFile(String sequenceFileName){
		
		try {
			if(!isFasta()){
				sequenceList = loadFastqSequences(sequenceFileName);
			} else {
				sequenceList = loadFastaSequences(sequenceFileName);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	//get all kmers container object
	public Hashtable<String, Integer> getTotalKmers(){
		
		return hashKmers;
		
	}
	//Count kmer in file per sequence
	public Hashtable<String, Integer> getFileKmers(String seq, int length)
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
	 
	public static void main (String [ ] args) {

		KmersCount kmerSpectrum = new KmersCount();
		
		if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
			CommandsDescriptor.getInstance().printHelp(KmersCount.class);
			return;
		}
		//Parameters
		int k=CommandsDescriptor.getInstance().loadOptions(kmerSpectrum, args);
		
		//Is fasta or fastq? and read it
		String sequenceFile = args[k];
		kmerSpectrum.processFile(sequenceFile);
	
		//Kmer Count Per File
		for(int i=0; i < kmerSpectrum.sequenceList.size(); i++){
			//Forward		
			String sequence = kmerSpectrum.sequenceList.get(i).getCharacters().toString();
			Hashtable<String, Integer> newHashKmers = kmerSpectrum.getFileKmers(sequence, kmerSpectrum.getKmerSize());
			kmerSpectrum.addKmersCount(newHashKmers);
			//Reverse complement
			if(kmerSpectrum.isBothStrands()){
				String reverseSequence = DNAMaskedSequence.getReverseComplement(sequence);
				newHashKmers = kmerSpectrum.getFileKmers(reverseSequence, kmerSpectrum.getKmerSize());
				kmerSpectrum.addKmersCount(newHashKmers);
			}
		}
		
		int maxValue = 0;
		int handicap = 0;
		
		//Make kmer spectrum
		for (Map.Entry<String, Integer> entry : kmerSpectrum.getTotalKmers().entrySet()) {
		    Integer value = entry.getValue();
	
		    if(kmerSpectrum.kmerDistribution.containsKey(value)){
		    	kmerSpectrum.kmerDistribution.put(value, kmerSpectrum.kmerDistribution.get(value) + 1);
		    	
		    	if(maxValue < kmerSpectrum.kmerDistribution.get(value) && value > 100){
		    		maxValue = kmerSpectrum.kmerDistribution.get(value);
		    		handicap = value;
		    	}
		    	
			} else {
				kmerSpectrum.kmerDistribution.put(value, 1);
			}
		}
		
		if(handicap < 100){
			handicap = 100;
		} else {
			handicap =  handicap + (int)Math.round(handicap * 1.5);
		}
		
		PrintStream out = System.out;
		out.println("Kmer_frequency\tNumber_of_distinct_kmers");
		
		for(int i= kmerSpectrum.getOmitFirstBases(); i < handicap; i++){
			if(kmerSpectrum.kmerDistribution.get(i) == null){
				out.println(i+"\t"+0);
			} else {
				out.println(i+"\t"+kmerSpectrum.kmerDistribution.get(i));
			}
			
		}
		
		
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
	public int getKmerSize() {
		return kmerSize;
	}
	public void setKmerSize(int kmerSize) {
		this.kmerSize = kmerSize;
	}
	public int getOmitFirstBases() {
		return omitFirstBases;
	}
	public void setOmitFirstBases(int omitFirstBases) {
		this.omitFirstBases = omitFirstBases;
	}

	 
	
}
