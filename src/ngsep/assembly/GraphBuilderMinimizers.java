package ngsep.assembly;

import java.util.ArrayList;
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

	private Logger log = Logger.getLogger(GraphBuilderMinimizers.class.getName());
	
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
		long startTime = System.currentTimeMillis();
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
		usedMemory/=1000000000;
		long time1 = System.currentTimeMillis();
		long diff1 = (time1-startTime)/1000;
		log.info("Total reads length: "+totalLength+" Mode: "+modeDepth+" Expected assembly length: "+expectedAssemblyLength+" Memory (Gbp): "+usedMemory+" Time(s): "+diff1);
		
		MinimizersTable table = new MinimizersTable(kmersAnalyzer, kmerLength, windowLength);
		table.setLog(log);
		//table.setMaxAbundanceMinimizer(Math.max(100, 5*modeDepth));
		
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
		usedMemory/=1000000000;
		long time2 = System.currentTimeMillis();
		diff1 = (time2-time1)/1000;
		long diff2 = (time2-startTime)/1000;
		log.info("Built minimizers. Memory(Gbp): "+usedMemory+" Time minimizers (s): "+diff1+" total time (s): "+diff2);
		Distribution minimizerHitsDist = table.calculateDistributionHits();
		minimizerHitsDist.printDistributionInt(System.out);
		
		AssemblyGraph graph = new AssemblyGraph(sequences);
		log.info("Created graph vertices. Edges: "+graph.getEdges().size());
		KmerHitsAssemblyEdgesFinder edgesFinder = new KmerHitsAssemblyEdgesFinder(graph);
		edgesFinder.setExpectedAssemblyLength(expectedAssemblyLength);
		List<List<AssemblySequencesRelationship>> relationshipsPerSequence = new ArrayList<List<AssemblySequencesRelationship>>(sequences.size());
		for(int i=0;i<sequences.size();i++) relationshipsPerSequence.add(null);
		ThreadPoolExecutor poolSearch = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for (int seqId = 0; seqId < sequences.size(); seqId++) {
			CharSequence seq = sequences.get(seqId).getCharacters();
			double compressionFactor = compressionFactors!=null?compressionFactors[seqId]:1;
			final int i = seqId;
			poolSearch.execute(()->processSequence(edgesFinder, table, i, seq, compressionFactor, relationshipsPerSequence));
			//if ((seqId+1)%1000==0) log.info("Scheduled sequence "+(seqId+1));
		}
		addRelationshipsToGraph(graph, relationshipsPerSequence, runtime);
		waitToFinish(finishTime, poolSearch);
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time3 = System.currentTimeMillis();
		diff1 = (time3-time2)/1000;
		diff2 = (time3-startTime)/1000;
		log.info("Built graph. Edges: "+graph.getEdges().size()+" Embedded: "+graph.getEmbeddedCount()+" Memory: "+usedMemory+" Time graph construction (s): "+diff1+" total time (s): "+diff2);
		//log.info(" Raw hits for "+edgesFinder.getCountRawHits()+" sequences. Completed hits for "+edgesFinder.getCountCompletedHits()+" sequences");
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
	
	private void processSequence(KmerHitsAssemblyEdgesFinder finder, MinimizersTable table, int seqId, CharSequence seq, double compressionFactor, List<List<AssemblySequencesRelationship>> relationshipsPerSequence ) {
		List<AssemblySequencesRelationship> answer;
		try {
			Map<Integer, Long> codesForward = KmersExtractor.extractDNAKmerCodes(seq.toString(), kmerLength, 0, seq.length());
			Map<Integer,List<UngappedSearchHit>> hitsForward = table.match(seqId, seq.length(), codesForward);
			String complement = DNAMaskedSequence.getReverseComplement(seq).toString();
			Map<Integer, Long> codesReverse = KmersExtractor.extractDNAKmerCodes(complement, kmerLength, 0, complement.length());
			Map<Integer,List<UngappedSearchHit>> hitsReverse = table.match(seqId, complement.length(), codesReverse);
			answer = finder.inferRelationshipsFromKmerHits(seqId, seq.toString(), complement, hitsForward, hitsReverse, compressionFactor);
			relationshipsPerSequence.set(seqId, answer);
			if ((seqId)%1000==0) {
				int edges = 0;
				int embedded = 0;
				for(AssemblySequencesRelationship next:answer) {
					if(next instanceof AssemblyEmbedded) embedded++;
					if(next instanceof AssemblyEdge) edges++;
				}
				log.info("Identified relationships for sequence "+(seqId) +" Candidate edges: "+edges+"  candidate embedded hosts "+embedded);
			}
		} catch (RuntimeException e) {
			relationshipsPerSequence.set(seqId, new ArrayList<AssemblySequencesRelationship>());
			throw e;
		}
	}
	private void addRelationshipsToGraph(AssemblyGraph graph, List<List<AssemblySequencesRelationship>> relationshipsPerSequence, Runtime runtime) {
		int n = relationshipsPerSequence.size();
		int i=0;
		log.info("Adding relationships to graph");
		while(i<n) {
			try {
				Thread.sleep(1);
			} catch (InterruptedException e) {
				throw new RuntimeException("Add relationships thread interrupted", e);
			}
			for(;i<n;i++) {
				List<AssemblySequencesRelationship> nextList = relationshipsPerSequence.get(i);
				if(nextList == null) break;
				//if ((i+1)%1000==0) log.info("Adding relationships for sequence "+(i+1) +" Relationships sequence: "+nextList.size());
				for(AssemblySequencesRelationship next:nextList) {
					if(next instanceof AssemblyEmbedded) graph.addEmbedded((AssemblyEmbedded)next);
					if(next instanceof AssemblyEdge) graph.addEdge((AssemblyEdge)next);
				}
				if(i == idxDebug) log.info("Edges start: "+graph.getEdges(graph.getVertex(i, true)).size()+" edges end: "+graph.getEdges(graph.getVertex(i, false)).size()+" Embedded: "+graph.getEmbeddedBySequenceId(i));
				if ((i+1)%10000==0) {
					long usedMemory = runtime.totalMemory()-runtime.freeMemory();
					usedMemory/=1000000000;
					log.info("Processed "+(i+1) +" sequences. Number of edges: "+graph.getNumEdges()+ " Embedded: "+graph.getEmbeddedCount()+" Memory: "+usedMemory);
				}
				//if ((seqId+1)%100==0) log.info("Processed "+(seqId+1) +" sequences. Number of edges: "+graph.getNumEdges()+ " Embedded: "+graph.getEmbeddedCount());
			}	
		}
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