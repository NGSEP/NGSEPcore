package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.ReadAlignment;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNASequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMap;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.io.FastaSequencesHandler;

public class TransposableElementsFinder {
	
	// Logging and progress
	private Logger log = Logger.getLogger(TransposableElementsFinder.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters 
	private String inputFile = null;
	private String outputFile = null;
	private String transposonsDatabaseFile = null;
	
	// Model attributes
	
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
	public String getTransposonsDatabaseFile() {
		return transposonsDatabaseFile;
	}
	public void setTransposonsDatabaseFile(String transposonsDatabaseFile) {
		this.transposonsDatabaseFile = transposonsDatabaseFile;
	}
	
	public static void main(String[] args) throws Exception {
		TransposableElementsFinder instance = new TransposableElementsFinder();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	public void run() throws IOException, InterruptedException {
		long time = System.currentTimeMillis();
		logParameters ();
		if(inputFile==null) throw new IOException("The input genome is a required parameter");
		if(outputFile==null) throw new IOException("The output file is a required parameter");
		ReferenceGenome genome = new ReferenceGenome(inputFile);
		List<TransposableElementAnnotation> transposonAnnotations = findTransposons(genome);
		saveTransposons(transposonAnnotations,outputFile);
		double seconds = (System.currentTimeMillis()-time);
		seconds /=1000;
		log.info("Process finished in "+seconds+" seconds");
	}
	public void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		if(transposonsDatabaseFile!=null) out.println("Database of transposable elements: "+ transposonsDatabaseFile);
		log.info(os.toString());
	}
	public List<TransposableElementAnnotation> findTransposons(ReferenceGenome genome) throws InterruptedException, IOException {
		GenomicRegionSortedCollection<TransposableElementAnnotation> annotations = new GenomicRegionSortedCollection<TransposableElementAnnotation>(genome.getSequencesMetadata());
		annotations.addAll(findTransposonsDeNovo(genome));
		if(transposonsDatabaseFile!=null) annotations.addAll(findTransposonsBySimilarity(genome));
		return removeRedundantAnnotations(annotations);
	}
	
	private List<TransposableElementAnnotation> findTransposonsDeNovo(ReferenceGenome genome) throws InterruptedException {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		
		QualifiedSequenceList list = genome.getSequencesList();
		KmersExtractor tableKmer = new KmersExtractor();
		tableKmer.processQualifiedSequences(list);
		KmersMap kmap= tableKmer.getKmersMap();
		int totalTranspoSize= 0;
		
		
	    //traverses the genome, chromosome by chromosome 
		for(QualifiedSequence chromosome:list)
		{
			//gives the chromosome sequence
			CharSequence chromosomeSequence = chromosome.getCharacters();
			int n = chromosome.getLength();
			int regionStart=0;
			int regionEnd=1;
			int regionFrequency=0;			
			//extract kmer by kmer of the entire genome
			for(int i=0;i<n-15;i++) 
			{
				int end = i+15;
				
				String kmer = chromosomeSequence.subSequence(i, end).toString();
				//check if kmer is DNA
				if(DNASequence.isDNA(kmer)) 
				{
					int kmerFrequency = kmap.getCount(kmer); 
					
					//filter kmers with given frequency
					if (kmerFrequency >=10)
					{ 
						int dif= Math.abs(regionEnd-i);
						// report regions
						// does the kmer belongs to the region?
						if (dif<=50)
						{
							regionEnd=end;
							
							// the frequency of the region will be the lowest frequency of the kmers
							if(regionFrequency==0)
							{
								regionFrequency= kmerFrequency;
							}
							if (kmerFrequency<regionFrequency)
							{
								regionFrequency= kmerFrequency;
							}
							// validate if the kmer is the last int the chromosome 
							if(regionEnd==n-1)
							{
								int sizeRegion= regionEnd-regionStart;
								totalTranspoSize+= sizeRegion;
								TransposableElementAnnotation transposon = new TransposableElementAnnotation(chromosome.getName(), regionStart, regionEnd); 
								answer.add(transposon);
								//System.out.println("" + chromosomeList.getName()+ "   " + regionStart + " "  + regionEnd + "   " +  regionFrequency + "  " + sizeRegion) ;
									
								// update 
								regionStart=0;
								regionEnd= 1;
								regionFrequency=0;	
							}
						}
						// find a new region	
						else 
						{
							int sizeRegion= regionEnd-regionStart;
							totalTranspoSize+= sizeRegion;
							TransposableElementAnnotation transposon = new TransposableElementAnnotation(chromosome.getName(), regionStart, regionEnd); 
							answer.add(transposon);
							//System.out.println("" + chromosomeList.getName()+ "   " + regionStart + " "  + regionEnd + "   " +  regionFrequency + "  " + sizeRegion) ;
								
							// update
							regionStart=i;
							regionEnd= end;
							regionFrequency=kmerFrequency;	
						}
					
					}
				}
			
			}
			System.out.println(totalTranspoSize); 	
		}

		return answer;
	}

	private List<TransposableElementAnnotation> findTransposonsBySimilarity(ReferenceGenome genome) throws InterruptedException, IOException {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
	
		System.out.println("creating minimizerTable");
		MinimizersTableReadAlignmentAlgorithm minimizerTable = new MinimizersTableReadAlignmentAlgorithm();
		minimizerTable.loadGenome(genome, 15, 30, 1);
		List<ReadAlignment> alignQueryToRef;
		System.out.println("loading known transposions");
		//load the fasta
		FastaSequencesHandler load = new FastaSequencesHandler();
		//loading known transposions
		List<QualifiedSequence> knownTransposons = load.loadSequences(transposonsDatabaseFile);
		System.out.println("searching transposions db");
		for (QualifiedSequence transposon:knownTransposons) {
			RawRead read = new RawRead(transposon.getName(), transposon.getCharacters(), null);
			 List<ReadAlignment> listRead = minimizerTable.alignRead(read);
			for (ReadAlignment ReadAlig: listRead) {
				System.out.println(ReadAlig);
			}
			
		}
		
		return answer;
	}
	private List<TransposableElementAnnotation> removeRedundantAnnotations(GenomicRegionSortedCollection<TransposableElementAnnotation> annotations) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		// TODO implement
		return answer;
	}
	public void saveTransposons(List<TransposableElementAnnotation> transposonAnnotations, String outputFile) {
		// TODO Auto-generated method stub
		
	}
	

}
