package ngsep.assembly;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
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
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
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
	public static final int RENAME_CONTIGS_POLICY_REFNAMES = 1;
	public static final int RENAME_CONTIGS_POLICY_CONSECUTIVE = 2;
	public static final int RENAME_CONTIGS_POLICY_KEEPNAMES = 0;
	public static final int DEF_NUM_THREADS = 1;
	
	// Logging and progress
	private Logger log = Logger.getLogger(AssemblyReferenceSorter.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	// Parameters
	private String inputFile;
	private String outputFile;
	private ReferenceGenome genome = null;
	private int ploidy = DEF_PLOIDY;
	private int kmerLength = DEF_KMER_LENGTH;
	private int windowLength = DEF_WINDOW_LENGTH;
	private int renameContigsPolicy = RENAME_CONTIGS_POLICY_REFNAMES;
	private int numThreads = DEF_NUM_THREADS;
	
	
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

	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	public void setGenome(String genomeFile) throws IOException {
		if(genomeFile==null || genomeFile.length()==0) this.genome = null;
		else setGenome(OptionValuesDecoder.loadGenome(genomeFile,log));
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
	
	public int getRenameContigsPolicy() {
		return renameContigsPolicy;
	}
	public void setRenameContigsPolicy(int renameContigsPolicy) {
		this.renameContigsPolicy = renameContigsPolicy;
	}
	public void setRenameContigsPolicy(String value) {
		setRenameContigsPolicy((int)OptionValuesDecoder.decode(value, Integer.class));
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

	public static void main(String[] args) throws Exception {
		AssemblyReferenceSorter instance = new AssemblyReferenceSorter();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}

	public void run() throws IOException {
		logParameters();
		if(inputFile==null) throw new IOException("The input genome assembly is required");
		if(genome==null) throw new IOException("The input reference assembly is required");
		if(outputFile==null) throw new IOException("An output file is required");
		try (PrintStream out = new PrintStream(outputFile)) {
			sortGenome(inputFile, out);
		}
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		out.println("Reference genome:"+ genome.getFilename());
		out.println("Kmer length: "+kmerLength);
		out.println("Window length for minimizers: "+windowLength);
		//out.println("Sample ploidy: "+ploidy);
		out.println("Policy to rename contigs: "+renameContigsPolicy);
		out.println("Number of threads: "+numThreads);
		log.info(os.toString());
	}
	public void sortGenome(String assemblyFile, PrintStream out) throws IOException {
		ReferenceGenome assembly = new ReferenceGenome(assemblyFile);
		ReferenceGenome answer = sortGenome(assembly);
		answer.saveGenome(out, 100);
	}

	public ReferenceGenome sortGenome(ReferenceGenome assembly) {
		MinimizersUngappedSearchHitsClustersFinder minimizerTable = new MinimizersUngappedSearchHitsClustersFinder();
		minimizerTable.loadGenome(genome, kmerLength, windowLength, numThreads);
		minimizerTable.setLog(log);
		QualifiedSequenceList refMetadata = genome.getSequencesMetadata();
		QualifiedSequenceList sequences = assembly.getSequencesList();
		Map<String,ReadAlignment> contigAlns = new Hashtable<String, ReadAlignment>();
		//TODO: Improve parallel processing
		ThreadPoolExecutor pool = new ThreadPoolExecutor(1, 1, 60*sequences.size(), TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
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
		GenomicRegionComparator cmp = new GenomicRegionComparator(refMetadata);
		List<ReadAlignment> alns = new ArrayList<ReadAlignment>(contigAlns.values());
		Collections.sort(alns,cmp);
		QualifiedSequenceList primarySeqs = new QualifiedSequenceList();
		QualifiedSequenceList redundantSeqs = new QualifiedSequenceList();
		GenomicRegionSortedCollection<ReadAlignment> referenceAlnsGenome = new GenomicRegionSortedCollection<ReadAlignment>(refMetadata);
		Set<String> namesMapped = new HashSet<String>();
		int numPrimarySeqsChr = 0;
		int lastChr = -1;
		for(ReadAlignment aln:alns) {
			String contigName = aln.getReadName();
			DNAMaskedSequence seq = (DNAMaskedSequence) assembly.getSequenceByName(contigName).getCharacters();
			if(aln.isNegativeStrand()) seq = seq.getReverseComplement();
			List<ReadAlignment> genomeOverlappingAlns = referenceAlnsGenome.findSpanningRegions(aln).asList();
			boolean found = false;
			for(ReadAlignment aln2:genomeOverlappingAlns) {
				int overlap = GenomicRegionSpanComparator.getInstance().getSpanLength(aln.getFirst(), aln.getLast(), aln2.getFirst(), aln2.getLast());
				if(overlap > 0.1*seq.length()) {
					redundantSeqs.add(assembly.getSequenceByName(contigName));
					found = true;
					break;
				}
			}
			if(!found) {
				if(lastChr!=aln.getSequenceIndex()) {
					lastChr = aln.getSequenceIndex();
					numPrimarySeqsChr = 0;
				}
				numPrimarySeqsChr++;
				primarySeqs.add(new QualifiedSequence(calculateName(contigName, aln, numPrimarySeqsChr, primarySeqs.size()+1, false), seq));
				referenceAlnsGenome.add(aln);
			}
			namesMapped.add(contigName);
		}
		log.info("Sorted "+primarySeqs.size()+" primary contigs and "+redundantSeqs.size()+" redundant contigs");
		QualifiedSequenceList outSeqs = new QualifiedSequenceList();
		outSeqs.addAll(primarySeqs);
		outSeqs.addAll(redundantSeqs);
		int n=0;
		for (QualifiedSequence seq:sequences) {
			if(!namesMapped.contains(seq.getName())) {
				outSeqs.add(seq);
				n++;
			}
		}
		log.info("Added "+n+" remaining contigs. Total: "+outSeqs.size());
		return new ReferenceGenome(outSeqs);
	}

	private String calculateName(String contigName, ReadAlignment aln, int consecutiveChr, int consecutive, boolean redundant) {
		if(renameContigsPolicy==RENAME_CONTIGS_POLICY_KEEPNAMES) return contigName;
		else if (renameContigsPolicy==RENAME_CONTIGS_POLICY_CONSECUTIVE) {
			return "SortedContig"+consecutive;
		}
		String name = aln.getSequenceName();
		name+="_c"+consecutiveChr;
		return name;
	}
	private void mapContig(QualifiedSequence qseq, QualifiedSequenceList refMetadata, MinimizersUngappedSearchHitsClustersFinder minimizerTable, Map<String, ReadAlignment> contigAlns) {
		ContigClustersComparator cmp = new ContigClustersComparator();
		String seqName = qseq.getName();
		log.info("Mapping contig "+seqName);
		DNAMaskedSequence contigSeq = (DNAMaskedSequence) qseq.getCharacters();
		List<UngappedSearchHitsCluster> forwardClusters = minimizerTable.findHitClusters(contigSeq.toString());
		log.info("Sequence: "+seqName+" Forward clusters: "+forwardClusters.size());
		//printClusters(forwardClusters);
		UngappedSearchHitsCluster bestClusterForward = null;
		if(forwardClusters.size()>0) {
			Collections.sort(forwardClusters, cmp);
			bestClusterForward = forwardClusters.get(0);
		}
		
		String reverseComplement =  DNAMaskedSequence.getReverseComplement(contigSeq).toString();
		List<UngappedSearchHitsCluster> reverseClusters = minimizerTable.findHitClusters(reverseComplement);
		log.info("Sequence: "+seqName+" Reverse clusters: "+reverseClusters.size());
		//printClusters(reverseClusters);
		UngappedSearchHitsCluster bestClusterReverse = null;
		if(reverseClusters.size()>0) {
			Collections.sort(reverseClusters, cmp);
			bestClusterReverse = reverseClusters.get(0);
		}
		UngappedSearchHitsCluster bestCluster = bestClusterForward;
		if (bestClusterForward == null && bestClusterReverse == null) return;
		else if (bestClusterForward ==null) bestCluster = bestClusterReverse;
		else if (bestClusterReverse!= null && cmp.compare(bestClusterForward, bestClusterReverse) >0 ) bestCluster = bestClusterReverse;
		ReadAlignment aln = buildAlignment(qseq, bestCluster, refMetadata, bestCluster==bestClusterReverse);
		System.out.println("Alignment of contig "+qseq.getName()+" with length "+qseq.getLength()+" "+aln);
		if(aln!=null) contigAlns.put(qseq.getName(), aln);
	}
	private ReadAlignment buildAlignment(QualifiedSequence qseq, UngappedSearchHitsCluster cluster, QualifiedSequenceList refMetadata, boolean reverse) {
		
		//System.out.println("Clusters: "+clusters.size()+" best weighted count: "+bestCluster.getWeightedCount());
		ReadAlignment aln = new ReadAlignment(cluster.getSubjectIdx(), cluster.getSubjectPredictedStart(), cluster.getSubjectPredictedEnd(), qseq.getLength(), 0);
		aln.setSequenceName(refMetadata.get(aln.getSequenceIndex()).getName());
		aln.setReadName(qseq.getName());
		if(reverse) aln.setNegativeStrand(true);
		double q = 1000*cluster.getWeightedCount()/qseq.getLength();
		q = Math.max(120, q);
		aln.setAlignmentQuality((byte) q);
		return aln;
	}
	public void printClusters(List<UngappedSearchHitsCluster> clusters) {
		for(UngappedSearchHitsCluster cluster: clusters) {
			System.out.print("Aln: "+cluster.getSubjectIdx()+" "+cluster.getSubjectPredictedStart()+" "+cluster.getSubjectPredictedEnd());
			System.out.print(" ev: "+cluster.getSubjectEvidenceStart()+" "+cluster.getSubjectEvidenceEnd()+" "+cluster.getRawKmerHits());
			System.out.println(" "+cluster.getNumDifferentKmers()+" "+cluster.getWeightedCount());
		}
		
	}
}
class ContigClustersComparator implements Comparator<UngappedSearchHitsCluster> {

	@Override
	public int compare(UngappedSearchHitsCluster o1, UngappedSearchHitsCluster o2) {
		int evLength1 = (o1.getSubjectEvidenceEnd()-o1.getSubjectEvidenceStart())/100000;
		int evLength2 = (o2.getSubjectEvidenceEnd()-o2.getSubjectEvidenceStart())/100000;
		if(evLength1!=evLength2) return evLength2-evLength1;
		if(o1.getNumDifferentKmers()!=o2.getNumDifferentKmers()) return o2.getNumDifferentKmers()-o1.getNumDifferentKmers();
		return (int)o2.getWeightedCount()-(int)o1.getWeightedCount();
	}
	
}
