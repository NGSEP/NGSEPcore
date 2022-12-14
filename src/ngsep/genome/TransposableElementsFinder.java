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
package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.alignments.MinimizersTableReadAlignmentAlgorithm;
import ngsep.alignments.UngappedSearchHitsCluster;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNASequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

/**
 * 
 * @author Daniela Lozano
 * @author Laura Gonzalez
 * @author Jorge Duitama
 *
 */
public class TransposableElementsFinder {
	
	// Logging and progress
	private Logger log = Logger.getLogger(TransposableElementsFinder.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	public static final int DEF_MIN_TE_LENGTH = 200;
	public static final int DEF_ROUNDS = 2;
	public static final int DEF_NUM_THREADS = 1;
	
	
	// Parameters 
	private String inputFile = null;
	private String outputFile = null;
	private String transposonsDatabaseFile = null;
	private int minTELength = DEF_MIN_TE_LENGTH;
	private int rounds = DEF_ROUNDS;
	private int numThreads = DEF_NUM_THREADS;
	
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
	
	public int getMinTELength() {
		return minTELength;
	}
	public void setMinTELength(int minTELength) {
		this.minTELength = minTELength;
	}
	public void setMinTELength(String value) {
		this.setMinTELength(Integer.parseInt(value));
	}
	
	public int getRounds() {
		return rounds;
	}
	public void setRounds(int rounds) {
		this.rounds = rounds;
	}
	public void setRounds(String value) {
		this.setRounds(Integer.parseInt(value));
	}
	
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	public void setNumThreads(String value) {
		this.setNumThreads(Integer.parseInt(value));
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
		if(transposonsDatabaseFile==null) throw new IOException("The transposons database file is a required parameter");
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
		out.println("Minimum TE length: "+minTELength);
		out.println("Number of search rounds: "+rounds);
		out.println("Number of threads: "+numThreads);
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
		System.out.println("Removing redundancies from final dataset with "+annotations.size()+" annotations");
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
	
		log.info("Creating reference index");
		MinimizersTableReadAlignmentAlgorithm minimizerTable = new MinimizersTableReadAlignmentAlgorithm();
		minimizerTable.loadGenome(genome, 15, 20, numThreads);
		minimizerTable.setMaxAlnsPerRead(10000);
		//minimizerTable.setMaxAlnsPerRead(100);
		//minimizerTable.setMaxAlnsPerRead(10);
		minimizerTable.setMinProportionBestCount(0);
		minimizerTable.setMinProportionReadLength(0);
		log.info("Loading known transposons");
		//load the fasta
		FastaSequencesHandler load = new FastaSequencesHandler();
		//loading known transposons
		List<QualifiedSequence> knownTransposons = load.loadSequences(transposonsDatabaseFile);
		log.info("Searching transposons db");
		List<TransposableElementAnnotation> elements = alignTransposonSequences(genome, minimizerTable, knownTransposons);
		answer.addAll(elements);
		if(answer.size() ==0 ) {
			log.info("Finished first round. No TEs found with the given database.");
			return answer;
		}
		log.info("Finished first round. Identified "+elements.size()+" regions.");
		for(int i=2;i<=rounds;i++) {
			//Second round
			log.info("Starting round "+i);
			GenomicRegionSortedCollection<TransposableElementAnnotation> sortedElements = new GenomicRegionSortedCollection<>(genome.getSequencesMetadata());
			sortedElements.addAll(elements);
			elements = removeRedundantAnnotations(sortedElements);
			log.info("Regions after removing redundancies: "+elements.size());
			List<QualifiedSequence> foundSequences = extractSequences(genome,elements);
			elements = alignTransposonSequences(genome, minimizerTable, foundSequences);
			answer.addAll(elements);
			log.info("Finished round "+i+". Identified "+elements.size()+" regions. Total: "+answer.size());
		}
		return answer;
	}
	private List<QualifiedSequence> extractSequences(ReferenceGenome genome, List<TransposableElementAnnotation> elements) {
		List<QualifiedSequence> sequences = new ArrayList<>();
		for(TransposableElementAnnotation ann:elements) {
			CharSequence sequence = genome.getReference(ann.getSequenceName(), ann.getFirst(), ann.getLast());
			QualifiedSequence qseq = new QualifiedSequence(ann.getTaxonomy(),sequence);
			sequences.add(qseq);
		}
		return sequences;
	}
	private List<TransposableElementAnnotation> alignTransposonSequences(ReferenceGenome genome, MinimizersTableReadAlignmentAlgorithm minimizerTable, List<QualifiedSequence> knownTransposons) {
		List<List<TransposableElementAnnotation>> predictions = new ArrayList<>();
		for(int i=0;i<knownTransposons.size();i++) {
			predictions.add(new ArrayList<>());
		}
		ThreadPoolExecutor pool = new ThreadPoolExecutor(numThreads, numThreads, knownTransposons.size(), TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for (int i=0;i<knownTransposons.size();i++) {
			QualifiedSequence transposon = knownTransposons.get(i);
			final int x = i;
			pool.execute(()->predictions.set(x,alignTransposonSequence(genome, minimizerTable, x, transposon)));
		}
		pool.shutdown();
		try {
			pool.awaitTermination(knownTransposons.size(), TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
    	if(!pool.isShutdown()) {
			throw new RuntimeException("The ThreadPoolExecutor was not shutdown after an await termination call");
		}
    	List<TransposableElementAnnotation> answer = new ArrayList<>();
    	for (int i=0;i<knownTransposons.size();i++) {
    		answer.addAll(predictions.get(i));
    	}
		return answer;
	}
	private List<TransposableElementAnnotation>  alignTransposonSequence(ReferenceGenome genome, MinimizersTableReadAlignmentAlgorithm minimizerTable,int seqId, QualifiedSequence transposon) {
		List<TransposableElementAnnotation> answer = new ArrayList<>();
		List<UngappedSearchHitsCluster> clusters= minimizerTable.buildHitClusters(transposon);
		//if(transposon.getName().contains("Chr2_20635304")) logClusters(genome, transposon, clusters);
		for (UngappedSearchHitsCluster cluster:clusters) {
			int sequenceIdx = cluster.getSubjectIdx();
			QualifiedSequence refSeq = genome.getSequenceByIndex(sequenceIdx);
			TransposableElementAnnotation alignedTransposon = new TransposableElementAnnotation(refSeq.getName(),cluster.getSubjectEvidenceStart(), cluster.getSubjectEvidenceEnd());
			if(alignedTransposon.length()<minTELength) continue;
			alignedTransposon.setTaxonomy(transposon.getName());
			answer.add(alignedTransposon);
		}
		if(seqId%100==0) log.info("Processed element index "+seqId+" name: "+transposon.getName()+" hits: "+answer.size());
		return answer;
	}
	private void logClusters(ReferenceGenome genome, QualifiedSequence transposon, List<UngappedSearchHitsCluster> clusters) {
		System.out.println("Transposon sequence id: "+transposon.getName()+" length: "+transposon.getLength()+" clusters: "+clusters.size());
		for (UngappedSearchHitsCluster cluster:clusters) {
			int sequenceIdx = cluster.getSubjectIdx();
			QualifiedSequence refSeq = genome.getSequenceByIndex(sequenceIdx);
			System.out.println("SeqName: "+refSeq.getName()+" expected: "+cluster.getSubjectPredictedStart() +" - "+ cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" - "+cluster.getSubjectEvidenceEnd()+" weight: "+cluster.getWeightedCount());
			
		}
		
	}
	/**
	 * Removes redundant transposons found by deNovo and similarity methods
	 * @param annotations the transposons found by deNovo and by similarity
	 * @return List<TransposableElementAnnotation> final list of unique transposable elements with chromosome, starting an ending position
	 */
	private List<TransposableElementAnnotation> removeRedundantAnnotations(GenomicRegionSortedCollection<TransposableElementAnnotation> annotations) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		TransposableElementAnnotation next = null;
		for(TransposableElementAnnotation ann:annotations) {
			//System.out.println("Next annotation. "+ann.getSequenceName()+":"+ann.getFirst()+" "+ann.getLast()+" family: "+ann.getTaxonomy());
			if(next==null) next = ann;
			else if (merge(next,ann)) {
				if(ann.length()>next.length()) next.setTaxonomy(ann.getTaxonomy());
				next.setLast(Math.max(next.getLast(),ann.getLast()));
			} else {
				answer.add(next);
				next = ann;
			}
		}
		if(next!=null)answer.add(next);
		
		return answer;
	}
	private boolean merge(TransposableElementAnnotation next, TransposableElementAnnotation ann) {
		if(!next.getSequenceName().equals(ann.getSequenceName())) return false;
		int spanLength = GenomicRegionSpanComparator.getInstance().getSpanLength(next.getFirst(), next.getLast(), ann.getFirst(), ann.getLast()); 
		return spanLength > next.length()/4 || spanLength > ann.length()/4;
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
	