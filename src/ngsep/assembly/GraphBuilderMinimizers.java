package ngsep.assembly;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.alignments.LongReadsAligner;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.MinimizersTable;

public class GraphBuilderMinimizers implements GraphBuilder {

	private Logger log = Logger.getLogger(GraphBuilderFMIndex.class.getName());
	
	public static final int DEF_WINDOW_LENGTH = 5;
	public static final int DEF_MIN_KMER_PCT = 20;
	public static final int DEF_NUM_THREADS = 1;
	
	private int kmerLength=KmersExtractor.DEF_KMER_LENGTH;
	private int windowLength=DEF_WINDOW_LENGTH;
	private int minKmerPercentage=DEF_MIN_KMER_PCT;
	private int numThreads = DEF_NUM_THREADS;
	
	private static final int TIMEOUT_SECONDS = 30;
	
	private LongReadsAligner aligner = new LongReadsAligner();
	
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
	public AssemblyGraph buildAssemblyGraph(List<CharSequence> sequences) {
		AssemblyGraph graph = new AssemblyGraph(sequences);
		log.info("Created graph vertices. Edges: "+graph.getEdges().size());
		
		KmerHitsAssemblyEdgesFinder edgesFinder = new KmerHitsAssemblyEdgesFinder(graph, minKmerPercentage);
		MinimizersTable table = new MinimizersTable(kmerLength, windowLength);
		int time = 2*sequences.size();
		ThreadPoolExecutor poolCreate =  new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		
		for(int seqId = 0; seqId < sequences.size(); seqId++) {
			CharSequence seq = sequences.get(seqId);
			poolCreate.execute (new Runnable() {
				@Override
				public void run() {
					table.addSequence(seq.toString());
				}
			});
		}
		waitToFinish(time, poolCreate);
		log.info("Built minimizers");
		ThreadPoolExecutor poolSearch = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		
		for (int seqId = 0; seqId < sequences.size(); seqId++) {
			CharSequence seq = sequences.get(seqId);
			if(numThreads==1) {
				processSequence(edgesFinder, table, seqId, seq);
				if ((seqId+1)%100==0) log.info("Processed "+(seqId+1) +" sequences. Number of edges: "+graph.getEdges().size()+ " Embedded: "+graph.getEmbeddedCount());
				continue;
			}
			Runnable task = new ProcessSequenceTask(this, edgesFinder, table, seqId, seq);
			poolSearch.execute(task);
		}
		
		waitToFinish(time, poolSearch);
		log.info("Built graph. Edges: "+graph.getEdges().size()+" Embedded: "+graph.getEmbeddedCount()+" Prunning embedded sequences");
		graph.pruneEmbeddedSequences();
		log.info("Prunned graph. Edges: "+graph.getEdges().size());
		//Create reverse map
		//traverse map to find matches
		
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
		int queryCount = table.getMinimizerCountBySequenceId(seqId);
		Map<Integer,Integer> minimizerCounts = table.countMinimizerHits(seqId, seq);
		Set<Integer> subjectIdxs = new HashSet<Integer>();
		int totalHits = 0;
		int maxSubjectCount = 0;
		for(int subjectIdx:minimizerCounts.keySet()) {
			int subjectCount = minimizerCounts.get(subjectIdx);
			if (seqId == idxDebug) System.out.println("GraphBuilderMinimizers. Query: "+seqId+" total: "+queryCount+" Subject sequence: "+subjectIdx+" hits: "+subjectCount);
			totalHits+=subjectCount;
			int minCount = minKmerPercentage*queryCount/100;
			if(subjectIdx< seqId && subjectCount>=minCount) {
				subjectIdxs.add(subjectIdx);
				maxSubjectCount = Math.max(maxSubjectCount, subjectCount);
			}
		}
		//Aproximate average minimizer hits
		double average = totalHits / queryCount;
		if(average<1) average = 1;
			
		Set<Integer> subjectFilteredIdxs = new HashSet<Integer>();
		for(int subjectIdx:subjectIdxs) {
			int subjectCount = minimizerCounts.get(subjectIdx);
			if(subjectIdx< seqId && subjectCount > 0.5*maxSubjectCount) subjectFilteredIdxs.add(subjectIdx);
		}
		if (seqId == idxDebug || subjectFilteredIdxs.size()>10) System.out.println("GraphBuilderMinimizers. Query: "+seqId+" Subject sequences: "+subjectFilteredIdxs.size()+" hits: "+totalHits+" average: "+average);
		
		updateGraph(finder, seqId, seq, false, subjectFilteredIdxs, average);
		CharSequence complement = DNAMaskedSequence.getReverseComplement(seq);
		updateGraph(finder, seqId, complement, true, subjectIdxs, average);
		AssemblyGraph graph = finder.getGraph();
		synchronized (graph) {
			graph.filterEdgesAndEmbedded (seqId);
		}
		if(seqId == idxDebug) System.out.println("Edges start: "+graph.getEdges(graph.getVertex(seqId, true)).size()+" edges end: "+graph.getEdges(graph.getVertex(seqId, false)).size()+" Embedded: "+graph.getEmbeddedBySequenceId(seqId));
	}
	private void updateGraph(KmerHitsAssemblyEdgesFinder finder, int querySequenceId, CharSequence query, boolean queryRC, Set<Integer> subjectIdxs, double averageHits) {
		Map<CharSequence,Integer> queryKmers = aligner.extractUniqueKmers(query, 0, query.length());
		AssemblyGraph graph = finder.getGraph();
		for (int subjectId:subjectIdxs) {
			CharSequence subject = graph.getSequence(subjectId);
			Map<CharSequence,Integer> subjectKmers = aligner.extractUniqueKmers(subject, 0, subject.length());
			List<UngappedSearchHit> hits = aligner.alignUniqueKmers(subjectId, subjectKmers, queryKmers);
			if (querySequenceId == idxDebug) System.out.println("GraphBuilderMinimizers. Query: "+querySequenceId+" "+queryRC+" Subject idx: "+subjectId+" kmers: "+subjectKmers.size()+" hits: "+hits.size());
			if(hits.size()==0) continue;
			List<KmerHitsCluster> subjectClusters = LongReadsAligner.clusterSequenceKmerAlns(querySequenceId, query, hits, 0);
			if (querySequenceId == idxDebug) System.out.println("Query: "+querySequenceId+" "+queryRC+" Subject idx: "+subjectId+" clusters: "+subjectClusters.size());
			finder.updateGraphWithKmerClusters(querySequenceId, queryRC, queryKmers.size(), averageHits, subjectClusters);
		}
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