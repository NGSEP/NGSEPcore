package ngsep.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.alignments.LongReadsAligner;
import ngsep.math.Distribution;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.MinimizersTableEntry;

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
		
		for(int seqId = 0; seqId < sequences.size(); seqId++) {
			CharSequence seq = sequences.get(seqId);
			table.addSequence(seqId, seq);
		}
		log.info("Built minimizers.");
		Distribution minimizerHitsDist = table.calculateDistributionHits();
		minimizerHitsDist.printDistributionInt(System.out);
		
		table.clearSingletonMinimizers();
		log.info("Minimizers after removing singletons: "+table.getTotalMinimizers());
		
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
		int finishTime = 2*sequences.size();
		waitToFinish(finishTime, poolSearch);
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
		updateGraph(finder, seqId, seq, false, table);
		CharSequence complement = DNAMaskedSequence.getReverseComplement(seq);
		updateGraph(finder, seqId, complement, true, table);
		AssemblyGraph graph = finder.getGraph();
		synchronized (graph) {
			graph.filterEdgesAndEmbedded (seqId);
		}
		if(seqId == idxDebug) log.info("Edges start: "+graph.getEdges(graph.getVertex(seqId, true)).size()+" edges end: "+graph.getEdges(graph.getVertex(seqId, false)).size()+" Embedded: "+graph.getEmbeddedBySequenceId(seqId));
	}
	private void updateGraph(KmerHitsAssemblyEdgesFinder finder, int querySequenceId, CharSequence query, boolean queryRC, MinimizersTable table) {
		Map<Integer, CharSequence> kmers = KmersExtractor.extractKmersAsMap(query, kmerLength, 1, 0, query.length(), false, true, true);
		Map<Integer, List<MinimizersTableEntry>> minimizersQuery = table.computeSequenceMinimizers(querySequenceId, query.length(), kmers);
		int queryCount = minimizersQuery.size();
		
		int minCount = minKmerPercentage*queryCount/100;
		if (querySequenceId == idxDebug) log.info("GraphBuilderMinimizers. Counting hits for query: "+querySequenceId+" queryCount: "+queryCount+" min count: "+minCount);
		Map<Integer,List<MinimizersTableEntry>> minimizerHitsBySubject = table.calculateMinimizerHits(querySequenceId, minimizersQuery);
		Set<Integer> subjectIdxs = new HashSet<Integer>();
		int totalHits = 0;
		int maxSubjectCount = 0;
		for(int subjectIdx:minimizerHitsBySubject.keySet()) {
			int subjectCount = minimizerHitsBySubject.get(subjectIdx).size();
			totalHits+=subjectCount;
			if(subjectIdx< querySequenceId && subjectCount>=minCount) {
				if (querySequenceId == idxDebug) System.out.println("GraphBuilderMinimizers. Query: "+querySequenceId+" total: "+queryCount+" Subject sequence: "+subjectIdx+" hits: "+subjectCount);
				subjectIdxs.add(subjectIdx);
				maxSubjectCount = Math.max(maxSubjectCount, subjectCount);
			}
		}
		//Aproximate average minimizer hits
		double averageHits = totalHits / queryCount;
		if(averageHits<1) averageHits = 1;
		if (querySequenceId == idxDebug) log.info("GraphBuilderMinimizers. Query: "+querySequenceId+" Subject sequences: "+subjectIdxs.size()+" hits: "+totalHits+" average: "+averageHits);
		for(int subjectIdx:subjectIdxs) {
			List<MinimizersTableEntry> subjectMinHits = minimizerHitsBySubject.get(subjectIdx);
			//Filter sequences with less than half of the maximum hits
			if(subjectMinHits.size()<maxSubjectCount/2) continue;
			Collections.sort(subjectMinHits,(h1,h2) -> h1.getStart()-h2.getStart());
			List<UngappedSearchHit> hits = new ArrayList<UngappedSearchHit>();
			for(MinimizersTableEntry entry: subjectMinHits) {
				int minimizer = entry.getMinimizer();
				List<MinimizersTableEntry> queryHits = minimizersQuery.get(minimizer);
				if (queryHits == null || queryHits.size() != 1) continue;
				
				MinimizersTableEntry queryEntry = queryHits.get(0);
				CharSequence kmer = kmers.get(queryEntry.getStart());
				if(kmer == null) {
					//Neighbor kmers normally share minimizers
					continue;
				}
				UngappedSearchHit kmerHit = new UngappedSearchHit(kmer, subjectIdx, entry.getStart());
				kmerHit.setQueryIdx(queryEntry.getStart());
				kmerHit.setTotalHitsQuery(table.getTotalHits(minimizer));
				hits.add(kmerHit);
					
			}
			if(hits.size()==0) continue;
			List<KmerHitsCluster> subjectClusters = LongReadsAligner.clusterSequenceKmerAlns(querySequenceId, query, hits, 0);
			if (querySequenceId == idxDebug) System.out.println("GraphBuilderMinimizers. Query: "+querySequenceId+" "+queryRC+" Subject idx: "+subjectIdx+" hits: "+hits.size()+" clusters: "+subjectClusters.size());
			finder.updateGraphWithKmerClusters(querySequenceId, queryRC, queryCount, averageHits, subjectClusters);
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