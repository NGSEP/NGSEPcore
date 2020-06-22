package ngsep.assembly;

import java.util.List;
import java.util.Map;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.math.Distribution;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.QualifiedSequence;

public class GraphBuilderMinimizers implements GraphBuilder {

	private Logger log = Logger.getLogger(GraphBuilderFMIndex.class.getName());
	
	public static final int DEF_WINDOW_LENGTH = 5;
	public static final int DEF_MIN_KMER_PCT = 10;
	public static final int DEF_NUM_THREADS = 1;
	
	private int kmerLength=KmersExtractor.DEF_KMER_LENGTH;
	private int windowLength=DEF_WINDOW_LENGTH;
	private int minKmerPercentage=DEF_MIN_KMER_PCT;
	private int numThreads = DEF_NUM_THREADS;
	
	private static final int TIMEOUT_SECONDS = 30;
	
	private static int idxDebug = -1;
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}

	
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
	}
	public int getWindowLength() {
		return windowLength;
	}
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}
	
	public int getMinKmerPercentage() {
		return minKmerPercentage;
	}
	public void setMinKmerPercentage(int minKmerPercentage) {
		this.minKmerPercentage = minKmerPercentage;
	}
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}

	@Override
	public AssemblyGraph buildAssemblyGraph(List<QualifiedSequence> sequences) {
		
		MinimizersTable table = new MinimizersTable(kmerLength, windowLength);
		//TODO: Make parameter
		table.setMaxAbundanceMinimizer(100);
		for(int seqId = 0; seqId < sequences.size(); seqId++) {
			CharSequence seq = sequences.get(seqId).getCharacters();
			table.addSequence(seqId, seq);
		}
		log.info("Built minimizers.");
		Distribution minimizerHitsDist = table.calculateDistributionHits();
		int meanDepth = (int) Math.round(minimizerHitsDist.getAverage());
		
		minimizerHitsDist.printDistributionInt(System.out);
		
		table.clearOverrepresentedMinimizers();
		log.info("Minimizers after removing overrepresented: "+table.getTotalMinimizers());
		
		AssemblyGraph graph = new AssemblyGraph(sequences);
		log.info("Created graph vertices. Edges: "+graph.getEdges().size());
		
		KmerHitsAssemblyEdgesFinder edgesFinder = new KmerHitsAssemblyEdgesFinder(graph);
		edgesFinder.setMinKmerPercentage(minKmerPercentage);
		edgesFinder.setMeanDepth(Math.max(5, meanDepth));
		ThreadPoolExecutor poolSearch = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		
		for (int seqId = 0; seqId < sequences.size(); seqId++) {
			CharSequence seq = sequences.get(seqId).getCharacters();
			if(numThreads==1) {
				processSequence(edgesFinder, table, seqId, seq);
				if ((seqId+1)%100==0) log.info("Processed "+(seqId+1) +" sequences. Number of edges: "+graph.getEdges().size()+ " Embedded: "+graph.getEmbeddedCount());
				continue;
			}
			Runnable task = new ProcessSequenceTask(this, edgesFinder, table, seqId, seq);
			poolSearch.execute(task);
		}
		int finishTime = 2*sequences.size();
		waitToFinish(finishTime, poolSearch);
		log.info("Built graph. Edges: "+graph.getEdges().size()+" Embedded: "+graph.getEmbeddedCount());
		return graph;
	}
	private void waitToFinish(int time, ThreadPoolExecutor poolSearch) {
		poolSearch.shutdown();
		try {
			poolSearch.awaitTermination(time, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
    	if(!poolSearch.isShutdown()) {
			throw new RuntimeException("The ThreadPoolExecutor was not shutdown after an await Termination call");
		}
	}
	void processSequence(KmerHitsAssemblyEdgesFinder finder, MinimizersTable table, int seqId, CharSequence seq) {
		Map<Integer,List<UngappedSearchHit>> hitsBySubjectIdx = table.match(seq);
		List<UngappedSearchHit> selfHits = hitsBySubjectIdx.get(seqId);
		int selfHitsCount = (selfHits!=null)?selfHits.size():1;
		finder.updateGraphWithKmerHitsMap(seqId, seq, false, selfHitsCount, hitsBySubjectIdx);
		CharSequence complement = DNAMaskedSequence.getReverseComplement(seq);
		finder.updateGraphWithKmerHitsMap(seqId, complement, true, selfHitsCount, table.match(complement));
		AssemblyGraph graph = finder.getGraph();
		if(seqId == idxDebug) log.info("Edges start: "+graph.getEdges(graph.getVertex(seqId, true)).size()+" edges end: "+graph.getEdges(graph.getVertex(seqId, false)).size()+" Embedded: "+graph.getEmbeddedBySequenceId(seqId));
	}
	
}
class ProcessSequenceTask implements Runnable {
	private GraphBuilderMinimizers parent;
	private KmerHitsAssemblyEdgesFinder finder;
	private MinimizersTable table;
	private int sequenceId;
	private CharSequence sequence;
	
	public ProcessSequenceTask(GraphBuilderMinimizers parent, KmerHitsAssemblyEdgesFinder finder, MinimizersTable table, int sequenceId, CharSequence sequence) {
		super();
		this.parent = parent;
		this.finder = finder;
		this.table = table;
		this.sequenceId = sequenceId;
		this.sequence = sequence;
	}
	@Override
	public void run() {
		// TODO Auto-generated method stub
		parent.processSequence(finder, table, sequenceId, sequence);
		AssemblyGraph graph = finder.getGraph();
		if ((sequenceId+1)%100==0) parent.getLog().info("Processed "+(sequenceId+1) +" sequences. Number of edges: "+graph.getNumEdges()+ " Embedded: "+graph.getEmbeddedCount());
	}
	
	
}