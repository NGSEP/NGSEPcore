package ngsep.kmerscount;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Hashtable;
import java.util.Map;

import ngsep.main.CommandsDescriptor;
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
				
		fis = new FileInputStream(filename);
		//Open buffer
		in = new BufferedReader(new InputStreamReader (fis));
		//Read sequences
		while(in.ready()){
			QualifiedSequence sequences = RawRead.load(in);
			answer.add(sequences);
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
		int k=0;
		
		while(k<args.length && args[k].charAt(0)=='-') {
			if ("-C".equals(args[k])) {
				kmerSpectrum.setBothStrands(true);
			} else if ("-mer".equals(args[k])) {
				k++;
				kmerSpectrum.setKmerSize(Integer.parseInt(args[k]));
			} else if ("-fasta".equals(args[k])) {
				kmerSpectrum.setFasta(true);
			} else if ("-omitFirstBases".equals(args[k])) {
				k++;
				kmerSpectrum.setOmitFirstBases(Integer.parseInt(args[k]));
			}  else {
				System.err.println("Unrecognized option: "+args[k]);
				CommandsDescriptor.getInstance().printHelp(KmersCount.class);
				return;
			}
			k++;
		}
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
				String reverseSequence = "";
				ReverseComplement reverseComplement = new ReverseComplement(sequence);
				reverseSequence = reverseComplement.makeReverseComplement();
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
		    	
		    	if(maxValue < kmerSpectrum.kmerDistribution.get(value) && value > 5){
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
			out.println(i+"\t"+kmerSpectrum.kmerDistribution.get(i));
		}
		
		
	}
	 
	boolean isBothStrands() {
		return bothStrands;
	}
	void setBothStrands(boolean bothStrands) {
		this.bothStrands = bothStrands;
	}
	boolean isFasta() {
		return fasta;
	}
	void setFasta(boolean fasta) {
		this.fasta = fasta;
	}
	int getKmerSize() {
		return kmerSize;
	}
	void setKmerSize(int kmerSize) {
		this.kmerSize = kmerSize;
	}
	int getOmitFirstBases() {
		return omitFirstBases;
	}
	void setOmitFirstBases(int omitFirstBases) {
		this.omitFirstBases = omitFirstBases;
	}

	 
	
}
