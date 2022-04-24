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
import ngsep.alignments.io.ReadAlignmentFileWriter;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
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
	/**
	 * Adds the transposons found by deNovo and by similarity
	 * @param genome to which the transposons will be found
	 * @return GenomicRegionSortedCollection object with the transposons information found by both methods
	 * @throws InterruptedException  if this exception is caught,it is because for the kmer table the genome could not be loaded, the transposons found by both methods could not be added
	 * @throws IOException catches the IOException thrown if the specified part of the transposon database file is locked or the path does not exist
	 */
	public List<TransposableElementAnnotation> findTransposons(ReferenceGenome genome) throws InterruptedException, IOException {
		GenomicRegionSortedCollection<TransposableElementAnnotation> annotations = new GenomicRegionSortedCollection<TransposableElementAnnotation>(genome.getSequencesMetadata());
		//annotations.addAll(findTransposonsDeNovo(genome));
		if(transposonsDatabaseFile!=null) annotations.addAll(findTransposonsBySimilarity(genome));
		return removeRedundantAnnotations(annotations);
	}
	/**
	 * Find deNovo transposons given a genome
	 * @param genome to which the transposons will be found 
	 * @return List<TransposableElementAnnotation> list of transposons found with chromosome, starting and ending position
	 * @throws InterruptedException if this exception is caught, it is because for the kmer table the genome could not be loaded
	 */
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
								//System.out.println("" + chromosome.getName()+ "   " + regionStart + " "  + regionEnd + "   " +  regionFrequency + "  " + sizeRegion) ;
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
							//System.out.println("" + chromosome.getName()+ "   " + regionStart + " "  + regionEnd + "   " +  regionFrequency + "  " + sizeRegion) ;	
							// update
							regionStart=i;
							regionEnd= end;
							regionFrequency=kmerFrequency;	
						}
					
					}
				}
			
			}
 	
		}
		System.out.println(totalTranspoSize); 	
		return answer;
	}
	/**
	 * find transposons by similarity given a genome and a transposon database
	 * @param genome to which the transposons will be found 
	 * @return List<TransposableElementAnnotation> list of transposons found with chromosome, starting and ending position
	 * @throws InterruptedException if this exception is caught, it is because for the minimizers table the genome could not be loaded
	 * @throws IOException catches the IOException thrown if the specified part of the transposon database file is locked or the path does not exist
	 */
	private List<TransposableElementAnnotation> findTransposonsBySimilarity(ReferenceGenome genome) throws InterruptedException, IOException {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
	
		System.out.println("creating minimizerTable");
		MinimizersTableReadAlignmentAlgorithm minimizerTable = new MinimizersTableReadAlignmentAlgorithm();
		minimizerTable.loadGenome(genome, 15, 20, 1);
		System.out.println("loading known transposions");
		//load the fasta
		FastaSequencesHandler load = new FastaSequencesHandler();
		//loading known transposions
		List<QualifiedSequence> knownTransposons = load.loadSequences(transposonsDatabaseFile);
		System.out.println("searching transposions db");
		for (QualifiedSequence transposon:knownTransposons) {
			RawRead read = new RawRead(transposon.getName(), transposon.getCharacters(), null);
			List<ReadAlignment> alignments = minimizerTable.alignRead(read);
			System.out.println("Transposon sequence id: "+transposon.getName()+" alignments: "+alignments);
			for (ReadAlignment aln: alignments) {
				TransposableElementAnnotation alignedTransposon = new TransposableElementAnnotation(aln.getSequenceName(),aln.getFirst(), aln.getLast());
				alignedTransposon.setTaxonomy(aln.getReadName());
				answer.add(alignedTransposon);
			}
		}		
		return answer;
	}
	/**
	 * removes redundant transposons found by deNovo and similarity methods
	 * @param annotations the transposons found by deNovo and by similarity
	 * @return List<TransposableElementAnnotation> final list of unique transposable elements with chromosome, starting an ending position
	 */
	private List<TransposableElementAnnotation> removeRedundantAnnotations(GenomicRegionSortedCollection<TransposableElementAnnotation> annotations) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		// TODO implement
		answer.addAll(annotations);
		return answer;
	}
	/**
	 * Save found the transposons 
	 * @param transposonAnnotations list of transposons found deNovo and similarity methods
	 * @param outputFile name of the output file where the transposon annotation will be saved
	 * @throws IOException catch the IOException thrown if the specified part of the file is locked or does not exist
	 */
	public void saveTransposons(List<TransposableElementAnnotation> transposonAnnotations, String outputFile) throws IOException {
		try (PrintStream outTransposon =  new PrintStream(outputFile)) {	
			for(TransposableElementAnnotation t:transposonAnnotations) 
			{
				outTransposon.print(t.getSequenceName()+"\t");
				outTransposon.print(t.getFirst()+"\t");
				outTransposon.print(t.getLast()+"\t");
				outTransposon.print(t.getTaxonomy());
				outTransposon.println();
			}
		}
	}
}
	