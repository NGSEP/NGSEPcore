package ngsep.assembly;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.alignments.MinimizersUngappedSearchHitsClustersFinder;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.UngappedSearchHitsCluster;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class AssemblyReferenceSorter {
	// Constants for default values
	public static final int DEF_PLOIDY = 1;
	public static final int DEF_KMER_LENGTH = 25;
	public static final int DEF_WINDOW_LENGTH = 40;
	public static final int DEF_NUM_THREADS = 1;
	
	// Logging and progress
	private Logger log = Logger.getLogger(AssemblyReferenceSorter.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	// Parameters
	private String inputFile;
	private String outputFile;
	private String referenceFile;
	private int ploidy = DEF_PLOIDY;
	private int kmerLength = DEF_KMER_LENGTH;
	private int windowLength = DEF_WINDOW_LENGTH;
	private int numThreads = DEF_NUM_THREADS;
	
	//Model attributes
	private ReferenceGenome referenceGenome;
	
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

	public String getReferenceFile() {
		return referenceFile;
	}
	public void setReferenceFile(String referenceFile) {
		this.referenceFile = referenceFile;
	}

	public int getPloidy() {
		return ploidy;
	}
	public void setPloidy(int ploidy) {
		this.ploidy = ploidy;
	}
	public void setPloidy(String value) {
		this.setPloidy((int) OptionValuesDecoder.decode(value, Integer.class));
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
		this.setNumThreads((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	public ReferenceGenome getReferenceGenome() {
		return referenceGenome;
	}
	public void setReferenceGenome(ReferenceGenome referenceGenome) {
		this.referenceGenome = referenceGenome;
	}

	public static void main(String[] args) throws Exception {
		AssemblyReferenceSorter instance = new AssemblyReferenceSorter();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}

	private void run() throws IOException {
		logParameters();
		referenceGenome = new ReferenceGenome(referenceFile);
		try (PrintStream out = new PrintStream(outputFile)) {
			sortGenome(inputFile, out);
		}
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		out.println("Reference file:"+ referenceFile);
		out.println("Kmer length: "+kmerLength);
		out.println("Window length for minimizers: "+windowLength);
		out.println("Sample ploidy: "+ploidy);
		out.println("Number of threads "+numThreads);
		log.info(os.toString());
	}
	public void sortGenome(String assemblyFile, PrintStream out) throws IOException {
		ReferenceGenome assembly = new ReferenceGenome(assemblyFile);
		ReferenceGenome answer = sortGenome(assembly);
		answer.saveGenome(out, 100);
	}

	public ReferenceGenome sortGenome(ReferenceGenome assembly) {
		MinimizersUngappedSearchHitsClustersFinder minimizerTable = new MinimizersUngappedSearchHitsClustersFinder();
		minimizerTable.loadGenome(referenceGenome, kmerLength, windowLength, numThreads);
		minimizerTable.setLog(log);
		QualifiedSequenceList refMetadata = referenceGenome.getSequencesMetadata();
		QualifiedSequenceList sequences = assembly.getSequencesList();
		Map<String,ReadAlignment> contigAlns = new Hashtable<String, ReadAlignment>();
		ThreadPoolExecutor pool = new ThreadPoolExecutor(numThreads, numThreads, 60*sequences.size(), TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for (QualifiedSequence seq:sequences) {
			pool.execute(()->mapContig(seq, refMetadata, minimizerTable, contigAlns));
		}
		pool.shutdown();
		try {
			pool.awaitTermination(60*sequences.size(), TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
    	if(!pool.isShutdown()) {
			throw new RuntimeException("The ThreadPoolExecutor was not shutdown after an await termination call");
		}
    	log.info("Mapped "+contigAlns.size()+" contigs");
    	GenomicRegionComparator cmp = new GenomicRegionComparator(referenceGenome.getSequencesMetadata());
    	List<ReadAlignment> alns = new ArrayList<ReadAlignment>(contigAlns.values());
    	Collections.sort(alns,cmp);
    	QualifiedSequenceList outSeqs = new QualifiedSequenceList();
    	Set<String> namesMapped = new HashSet<String>();
    	for(ReadAlignment aln:alns) {
    		String contigName = aln.getReadName();
    		DNAMaskedSequence seq = (DNAMaskedSequence) assembly.getSequenceByName(contigName).getCharacters();
    		if(aln.isNegativeStrand()) seq = seq.getReverseComplement();
    		outSeqs.add(new QualifiedSequence(contigName, seq));
    		namesMapped.add(contigName);
    	}
    	log.info("Sorted "+outSeqs.size()+" contigs");
    	for (QualifiedSequence seq:sequences) {
    		if(!namesMapped.contains(seq.getName())) {
    			outSeqs.add(seq);
    		}
    	}
    	log.info("Added remaining contigs. Total: "+outSeqs.size());
		return new ReferenceGenome(outSeqs);
	}

	private void mapContig(QualifiedSequence qseq, QualifiedSequenceList refMetadata, MinimizersUngappedSearchHitsClustersFinder minimizerTable, Map<String, ReadAlignment> contigAlns) {
		String seqName = qseq.getName();
		log.info("Mapping contig "+seqName);
		DNAMaskedSequence contigSeq = (DNAMaskedSequence) qseq.getCharacters();
		List<UngappedSearchHitsCluster> forwardClusters = minimizerTable.buildHitClusters(contigSeq.toString(),false,true);
		ReadAlignment alnF = buildAlignment(qseq, forwardClusters, refMetadata, false);
		log.info("Sequence: "+seqName+" Forward clusters: "+forwardClusters.size()+" alignment: "+alnF);
		String reverseComplement =  DNAMaskedSequence.getReverseComplement(contigSeq).toString();
		List<UngappedSearchHitsCluster> reverseClusters = minimizerTable.buildHitClusters(reverseComplement,false,true);
		ReadAlignment alnR = buildAlignment(qseq, reverseClusters, refMetadata, true);
		log.info("Sequence: "+seqName+" Reverse clusters: "+reverseClusters.size()+" alignment: "+alnR);
		ReadAlignment aln = alnF;
		if (alnF==null && alnR==null) return;
		else if (alnF==null || (alnR!=null && alnR.getAlignmentQuality()>alnF.getAlignmentQuality())) aln = alnR;
		contigAlns.put(qseq.getName(), aln);
	}

	private ReadAlignment buildAlignment(QualifiedSequence qseq, List<UngappedSearchHitsCluster> clusters, QualifiedSequenceList refMetadata, boolean reverse) {
		if(clusters.size()==0) return null;
		Collections.sort(clusters, (o1,o2)-> ((int)o2.getWeightedCount())-((int)o1.getWeightedCount()));
		UngappedSearchHitsCluster bestCluster = clusters.get(0);
		ReadAlignment aln = new ReadAlignment(bestCluster.getSubjectIdx(), bestCluster.getSubjectPredictedStart(), bestCluster.getSubjectPredictedEnd(), qseq.getLength(), 0);
		aln.setSequenceName(refMetadata.get(aln.getSequenceIndex()).getName());
		aln.setReadName(qseq.getName());
		if(reverse) aln.setNegativeStrand(true);
		double q = 1000*bestCluster.getWeightedCount()/qseq.getLength();
		q = Math.max(120, q);
		aln.setAlignmentQuality((byte) q);
		return aln;
	}
}
