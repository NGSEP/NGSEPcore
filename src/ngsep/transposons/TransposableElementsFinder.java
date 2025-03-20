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
package ngsep.transposons;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.alignments.MinimizersUngappedSearchHitsClustersFinder;
import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.alignments.PairwiseAlignment;
import ngsep.alignments.UngappedSearchHitsCluster;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.transposons.io.TransposableElementLibraryHandler;

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
	private double minWeightedCount = 10;
	private int maxAlnsPerTransposon = 10000;
	private double minProportionBestCount = 0;
	
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
		saveTransposonAnnotations(transposonAnnotations,outputFile);
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
		DeNovoTransposableElementsFinder deNovoFinder = new DeNovoTransposableElementsFinderKmerAbundance();
		List<TransposableElement> knownElements = deNovoFinder.findTransposons(genome);
		if(transposonsDatabaseFile!=null) knownElements.addAll(loadKnownTransposons());
		return findTransposonsBySimilarity(genome,knownElements);
	}
	private List<TransposableElement> loadKnownTransposons() throws IOException {
		log.info("Loading known transposons");
		TransposableElementLibraryHandler handler = new TransposableElementLibraryHandler();
		return handler.load(transposonsDatabaseFile);
	}
	/**
	 * find transposons by similarity given a genome and a transposon database
	 * @param genome to which the transposons will be found 
	 * @return List<TransposableElementAnnotation> list of transposons found with chromosome, starting and ending position
	 * @throws InterruptedException if this exception is caught, it is because for the minimizers table the genome could not be loaded
	 * @throws IOException catches the IOException thrown if the specified part of the transposon database file is locked or the path does not exist
	 */
	private List<TransposableElementAnnotation> findTransposonsBySimilarity(ReferenceGenome genome, List<TransposableElement> knownTEs) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
	
		log.info("Creating reference index");
		MinimizersUngappedSearchHitsClustersFinder minimizerTable = new MinimizersUngappedSearchHitsClustersFinder();
		minimizerTable.loadGenome(genome, 15, 20, numThreads);
		minimizerTable.setMinProportionReadLength(0);
		
		log.info("Finished first round. Identified "+answer.size()+" regions.");
		for(int i=1;i<=rounds;i++) {
			log.info("Starting round "+i);
			log.info("Searching transposons db");
			answer = alignTransposonSequences(genome, minimizerTable, knownTEs);
			if(answer.size() ==0 ) {
				log.info("Finished round. No TEs found with the given database.");
				return answer;
			}
			log.info("Finished round "+i+". Identified "+answer.size()+" regions. Total events: "+answer.size());
			if(i<rounds) knownTEs = extractTEs(genome,answer);
		}
		return answer;
	}
	private List<TransposableElement> extractTEs(ReferenceGenome genome, List<TransposableElementAnnotation> elements) {
		List<TransposableElement> sequences = new ArrayList<>();
		for(TransposableElementAnnotation ann:elements) {
			CharSequence sequence = genome.getReference(ann.getSequenceName(), ann.getFirst(), ann.getLast());
			if(sequence==null) {
				log.warning("No sequence for annotation at "+ann.getSequenceName()+": "+ann.getFirst()+" "+ann.getLast());
				continue;
			}
			TransposableElement te = new TransposableElement(ann.getSequenceName()+"_"+ann.getFirst()+"_"+ann.getLast()+"#"+ann.getTaxonomy(),sequence);
			if(ann.getInferredFamily()!=null) te.setFamily(ann.getInferredFamily());
			else te.setFamily(ann.getSourceFamily());
			sequences.add(te);
		}
		return sequences;
	}
	private List<TransposableElementAnnotation> alignTransposonSequences(ReferenceGenome genome, MinimizersUngappedSearchHitsClustersFinder minimizerTable, List<TransposableElement> knownTransposons) {
		List<List<TransposableElementAnnotation>> predictions = new ArrayList<>();
		for(int i=0;i<knownTransposons.size();i++) {
			predictions.add(new ArrayList<>());
		}
		try (ThreadPoolExecutor pool = new ThreadPoolExecutor(numThreads, numThreads, knownTransposons.size(), TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>())) {
			for (int i=0;i<knownTransposons.size();i++) {
				TransposableElement transposon = knownTransposons.get(i);
				final int x = i;
				pool.execute(()->predictions.set(x,alignTransposonSequence(genome, minimizerTable, x, transposon)));
			}
			pool.shutdown();
			if(!pool.isShutdown()) {
				throw new RuntimeException("The ThreadPoolExecutor was not shutdown after an await termination call");
			}
		}
		GenomicRegionSortedCollection<TransposableElementAnnotation> answer = new GenomicRegionSortedCollection<TransposableElementAnnotation>(genome.getSequencesMetadata());
    	for (int i=0;i<knownTransposons.size();i++) {
    		answer.addAll(predictions.get(i));
    	}
		return removeRedundantAnnotations(genome,answer);
	}
	private List<TransposableElementAnnotation>  alignTransposonSequence(ReferenceGenome genome, MinimizersUngappedSearchHitsClustersFinder minimizerTable,int seqId, TransposableElement transposon) {
		boolean debug = false;
		//boolean debug = transposon.getName().contains("Pl05_19834351_19842561");
		List<TransposableElementAnnotation> rawHits = new ArrayList<>();
		String readSeq = transposon.getSequence().toString();
		List<UngappedSearchHitsCluster> forwardClusters = minimizerTable.buildHitClusters(readSeq,true);
		double maxCount = summarize(forwardClusters);
		if(debug) logClusters(genome, transposon, forwardClusters);
		for (UngappedSearchHitsCluster cluster:forwardClusters) rawHits.add(buildTEAnnotation(genome, cluster, false,transposon));
		String reverseComplement =  DNAMaskedSequence.getReverseComplement(readSeq).toString();
		List<UngappedSearchHitsCluster> reverseClusters = minimizerTable.buildHitClusters(reverseComplement,true);
		maxCount = Math.max(maxCount, summarize(reverseClusters));
		if(debug) logClusters(genome, transposon, reverseClusters);
		for (UngappedSearchHitsCluster cluster:reverseClusters) rawHits.add(buildTEAnnotation(genome, cluster, true,transposon));
		Collections.sort(rawHits, (o1,o2)-> ((int)o2.getCount())-((int)o1.getCount()));
		List<TransposableElementAnnotation> answer = new ArrayList<>();
		for (int i=0;i<rawHits.size() && i<maxAlnsPerTransposon;i++) {
			TransposableElementAnnotation ann = rawHits.get(i);
			if(ann.length()<minTELength) continue;
			double wc = ann.getCount();
			if(debug) System.out.println("Next candidate TE "+ann.getSequenceName()+": "+ann.getFirst()+" "+ann.getLast()+" count: "+wc+" min: "+minWeightedCount);
			if(wc<minWeightedCount || wc<minProportionBestCount*maxCount) break;
			
			answer.add(ann);
		}
		
		if(seqId%100==0) log.info("Processed element index "+seqId+" name: "+transposon.getId()+" hits: "+answer.size());
		return answer;
	}
	private double summarize(List<UngappedSearchHitsCluster> clusters) {
		double maxCount = 0;
		for (UngappedSearchHitsCluster cluster:clusters) {
			cluster.summarize();
			maxCount = Math.max(maxCount,cluster.getCountKmerHitsCluster());
			//System.out.println("Summarizing clusters. Next cluster "+cluster.getSubjectIdx()+": "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" hits: "+cluster.getNumDifferentKmers()+" count: "+cluster.getWeightedCount()+" maxCount: "+maxCount);
		}
		return maxCount;
	}
	private TransposableElementAnnotation buildTEAnnotation(ReferenceGenome genome, UngappedSearchHitsCluster cluster, boolean negativeStrand, TransposableElement transposon) {
		int sequenceIdx = cluster.getSubjectIdx();
		QualifiedSequence refSeq = genome.getSequenceByIndex(sequenceIdx);
		int first = Math.max(1, cluster.getSubjectEvidenceStart());
		int last = Math.min(refSeq.getLength()-1, cluster.getSubjectEvidenceEnd());
		TransposableElementAnnotation alignedTransposon = new TransposableElementAnnotation(refSeq.getName(),first, last);
		
		alignedTransposon.setNegativeStrand(negativeStrand);
		alignedTransposon.setCount(cluster.getCountKmerHitsCluster());
		alignedTransposon.setSource(transposon);
		return alignedTransposon;
	}
	private void logClusters(ReferenceGenome genome, TransposableElement transposon, List<UngappedSearchHitsCluster> clusters) {
		System.out.println("Transposon sequence id: "+transposon.getId()+" length: "+transposon.getSequence().length()+" clusters: "+clusters.size());
		for (UngappedSearchHitsCluster cluster:clusters) {
			int sequenceIdx = cluster.getSubjectIdx();
			QualifiedSequence refSeq = genome.getSequenceByIndex(sequenceIdx);
			System.out.println("SeqName: "+refSeq.getName()+" expected: "+cluster.getSubjectPredictedStart() +" - "+ cluster.getSubjectPredictedEnd()+" evidence: "+cluster.getSubjectEvidenceStart()+" - "+cluster.getSubjectEvidenceEnd()+" count: "+cluster.getCountKmerHitsCluster());
			
		}
		
	}
	/**
	 * Removes redundant transposons found by deNovo and similarity methods
	 * @param annotations the transposons found by deNovo and by similarity
	 * @return List<TransposableElementAnnotation> final list of unique transposable elements with chromosome, starting an ending position
	 */
	private List<TransposableElementAnnotation> removeRedundantAnnotations(ReferenceGenome genome, GenomicRegionSortedCollection<TransposableElementAnnotation> annotations) {
		List<TransposableElementAnnotation> answer = new ArrayList<TransposableElementAnnotation>();
		TransposableElementAnnotation next = null;
		for(TransposableElementAnnotation ann:annotations) {
			//System.out.println("Next annotation. "+ann.getSequenceName()+":"+ann.getFirst()+" "+ann.getLast()+" family: "+ann.getTaxonomy());
			if(next==null) next = ann;
			else if (merge(next,ann)) {
				if(next.getInferredFamily()==null && ann.length()>next.length()) {
					next.setSource(ann.getSource());
				}
				next.setLast(Math.max(next.getLast(),ann.getLast()));
				//assignFamily(next, genome);
			} else {
				answer.add(next);
				next = ann;
			}
		}
		if(next!=null)answer.add(next);
		
		return answer;
	}
	private void assignFamily(TransposableElementAnnotation ann, ReferenceGenome genome) {
		if(ann.length()<5000) return;
		HMMTransposonDomainsFinder domainsFinder = new HMMTransposonDomainsFinder();
		domainsFinder.loadHMMsFromClasspath();
		CharSequence seq = genome.getReference(ann);
		TransposableElementFamily family = domainsFinder.assignFamily((DNAMaskedSequence) seq);
		if(family==null) return;
		ann.setInferredFamily(family);
	}
	private boolean merge(TransposableElementAnnotation next, TransposableElementAnnotation ann) {
		if(!next.getSequenceName().equals(ann.getSequenceName())) return false;
		int spanLength = GenomicRegionSpanComparator.getInstance().getSpanLength(next.getFirst(), next.getLast(), ann.getFirst(), ann.getLast()); 
		double p1 = 1.0*spanLength/next.length();
		double p2 = 1.0*spanLength/ann.length();
		if(p1>0.95 || p2>0.95) return true;
		return p1 > 0.7 && p2>0.7 && Math.max(ann.getLast(), next.getLast())-Math.min(ann.getFirst(), next.getFirst())<20000;
	}
	/**
	 * Save transposons 
	 * @param transposonAnnotations list of transposons found deNovo and similarity methods
	 * @param outputFile name of the output file where the transposon annotation will be saved
	 * @throws IOException catch the IOException thrown if the specified part of the file is locked or does not exist
	 */
	public void saveTransposonAnnotations(List<TransposableElementAnnotation> transposonAnnotations, String outputFile) throws IOException {
		try (PrintStream outTransposon =  new PrintStream(outputFile)) {	
			for(TransposableElementAnnotation t:transposonAnnotations) 
			{
				outTransposon.print(t.getSequenceName());
				outTransposon.print("\t"+t.getFirst());
				outTransposon.print("\t"+t.getLast()+"\t");
				outTransposon.print("\t"+(t.isNegativeStrand()?"-":"+"));
				outTransposon.print("\t"+t.getTaxonomy());
				outTransposon.print("\t"+t.getQueryName());
				outTransposon.println();
			}
		}
	}
}
	