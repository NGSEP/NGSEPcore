package ngsep.assembly;

import java.util.List;
import java.util.Map;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.math.Distribution;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMap;
import ngsep.sequences.KmersMapAnalyzer;
import ngsep.sequences.MinimizersTable;
import ngsep.sequences.QualifiedSequence;

public class GraphBuilderMinimizers implements GraphBuilder {

	private Logger log = Logger.getLogger(GraphBuilderFMIndex.class.getName());
	
	public static final int DEF_WINDOW_LENGTH = 10;
	public static final int DEF_NUM_THREADS = 1;
	
	private int kmerLength=KmersExtractor.DEF_KMER_LENGTH;
	private int windowLength=DEF_WINDOW_LENGTH;
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
	
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}

	@Override
	public AssemblyGraph buildAssemblyGraph(List<QualifiedSequence> sequences) {
		return buildAssemblyGraph(sequences,null);
	}
		
	public AssemblyGraph buildAssemblyGraph(final List<QualifiedSequence> sequences, final double [] compressionFactors) {
		Runtime runtime = Runtime.getRuntime();
		log.info("Calculating kmers distribution");
		KmersExtractor extractor = new KmersExtractor();
		extractor.setLog(log);
		//The conditional avoids creating twice the large array in ShortArrayKmersMapImpl
		if(extractor.getKmerLength()!=kmerLength) extractor.setKmerLength(kmerLength);
		extractor.initializeMap();
		long totalLength =  0;
		int finishTime = sequences.size();
		ThreadPoolExecutor poolKmers = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for(int seqId = 0; seqId < sequences.size(); seqId++) {
			QualifiedSequence seq = sequences.get(seqId);
			totalLength+=seq.getLength();
			if(numThreads==1) {
				countSequenceKmers(extractor, seqId, seq);
			} else {
				final int i = seqId;
				poolKmers.execute(()->countSequenceKmers(extractor, i, seq));
			}
		}
		waitToFinish(finishTime, poolKmers);
		KmersMap map = extractor.getKmersMap();
		KmersMapAnalyzer kmersAnalyzer = new KmersMapAnalyzer(map, false);
		int modeDepth = kmersAnalyzer.getMode();
		long expectedAssemblyLength = kmersAnalyzer.getExpectedAssemblyLength();
		
		if(compressionFactors!=null) {
			double averageCompression = NumberArrays.getAverage(compressionFactors);
			expectedAssemblyLength/= averageCompression;
		}
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		log.info("Total reads length: "+totalLength+" Mode: "+modeDepth+" Expected assembly length: "+expectedAssemblyLength+" Memory: "+usedMemory);
		
		MinimizersTable table = new MinimizersTable(kmersAnalyzer, kmerLength, windowLength);
		table.setLog(log);
		//table.setMaxAbundanceMinimizer(Math.max(100, 5*modeDepth));
		//int firstIdNoGraph = sequences.size();
		ThreadPoolExecutor poolMinimizers = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for(int seqId = 0; seqId < sequences.size(); seqId++) {
			QualifiedSequence qseq = sequences.get(seqId);
			CharSequence seq = qseq.getCharacters();
			if(numThreads==1) {
				addSequenceToTable(table, seqId, seq);
			} else {
				final int i = seqId;
				poolMinimizers.execute(()->addSequenceToTable(table, i, seq));
			}
		}
		waitToFinish(finishTime, poolMinimizers);
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		log.info("Built minimizers. Memory: "+usedMemory);
		Distribution minimizerHitsDist = table.calculateDistributionHits();
		minimizerHitsDist.printDistributionInt(System.out);
		
		AssemblyGraph graph = new AssemblyGraph(sequences);
		log.info("Created graph vertices. Edges: "+graph.getEdges().size());
		
		KmerHitsAssemblyEdgesFinder edgesFinder = new KmerHitsAssemblyEdgesFinder(graph);
		ThreadPoolExecutor poolSearch = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for (int seqId = 0; seqId < sequences.size(); seqId++) {
			CharSequence seq = sequences.get(seqId).getCharacters();
			double compressionFactor = compressionFactors!=null?compressionFactors[seqId]:1;
			if(numThreads==1) {
				processSequence(edgesFinder, table, seqId, seq, compressionFactor);
			} else {
				final int i = seqId;
				poolSearch.execute(()->processSequence(edgesFinder, table, i, seq, compressionFactor));
			}
		}
		waitToFinish(finishTime, poolSearch);
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		log.info("Built graph. Edges: "+graph.getEdges().size()+" Embedded: "+graph.getEmbeddedCount()+" Memory: "+usedMemory);
		return graph;
	}
	public void countSequenceKmers(KmersExtractor extractor, int seqId, QualifiedSequence seq) {
		extractor.countSequenceKmers(seq);
		if ((seqId+1)%1000==0) log.info("Kmers extracted for "+(seqId+1)+" sequences.");
	}
	public void addSequenceToTable(MinimizersTable table, int seqId, CharSequence seq) {
		table.addSequence(seqId, seq);
		if ((seqId+1)%1000==0) log.info("Processed "+(seqId+1)+" sequences. Total minimizers: "+table.size()+" total entries: "+table.getTotalEntries());
	}
	
	private void processSequence(KmerHitsAssemblyEdgesFinder finder, MinimizersTable table, int seqId, CharSequence seq, double compressionFactor) {
		Map<Integer,List<UngappedSearchHit>> hitsBySubjectIdx = table.match(seqId, seq);
		
		List<UngappedSearchHit> selfHits = hitsBySubjectIdx.get(seqId);
		int selfHitsCount = (selfHits!=null)?selfHits.size():1;
		finder.updateGraphWithKmerHitsMap(seqId, seq, false, compressionFactor, selfHitsCount, kmerLength, hitsBySubjectIdx);
		CharSequence complement = DNAMaskedSequence.getReverseComplement(seq);
		hitsBySubjectIdx = table.match(seqId, complement);
		finder.updateGraphWithKmerHitsMap(seqId, complement, true, compressionFactor, selfHitsCount, kmerLength, hitsBySubjectIdx);
		AssemblyGraph graph = finder.getGraph();
		/*synchronized (graph) {
			graph.filterEmbedded(seqId, 0.5, 10);
		}*/
		if(seqId == idxDebug) log.info("Edges start: "+graph.getEdges(graph.getVertex(seqId, true)).size()+" edges end: "+graph.getEdges(graph.getVertex(seqId, false)).size()+" Embedded: "+graph.getEmbeddedBySequenceId(seqId));
		if ((seqId+1)%1000==0) log.info("Processed "+(seqId+1) +" sequences. Number of edges: "+graph.getNumEdges()+ " Embedded: "+graph.getEmbeddedCount());
	}
	private void waitToFinish(int time, ThreadPoolExecutor pool) {
		pool.shutdown();
		try {
			pool.awaitTermination(time, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
    	if(!pool.isShutdown()) {
			throw new RuntimeException("The ThreadPoolExecutor was not shutdown after an await termination call");
		}
	}
}