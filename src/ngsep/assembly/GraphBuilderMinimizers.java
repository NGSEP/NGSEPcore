package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.math.NumberArrays;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.UngappedSearchHit;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMap;
import ngsep.sequences.KmersMapAnalyzer;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.ShortKmerCodesTable;

public class GraphBuilderMinimizers implements GraphBuilder {

	private Logger log = Logger.getLogger(GraphBuilderMinimizers.class.getName());
	
	public static final int DEF_WINDOW_LENGTH = 10;
	public static final int DEF_NUM_THREADS = 1;
	
	private int kmerLength=KmersExtractor.DEF_KMER_LENGTH;
	private int windowLength=DEF_WINDOW_LENGTH;
	private int ploidy = AssemblyGraph.DEF_PLOIDY_ASSEMBLY;
	private int numThreads = DEF_NUM_THREADS;
	private KmersMap kmersMap;
	
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
	
	public int getPloidy() {
		return ploidy;
	}
	public void setPloidy(int ploidy) {
		this.ploidy = ploidy;
	}
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	
	public KmersMap getKmersMap() {
		return kmersMap;
	}
	public void setKmersMap(KmersMap kmersMap) {
		this.kmersMap = kmersMap;
	}
	@Override
	public AssemblyGraph buildAssemblyGraph(List<QualifiedSequence> sequences) {
		return buildAssemblyGraph(sequences,null);
	}
		
	public AssemblyGraph buildAssemblyGraph(final List<QualifiedSequence> sequences, final double [] compressionFactors) {
		Runtime runtime = Runtime.getRuntime();
		
		KmersMapAnalyzer kmersAnalyzer = new KmersMapAnalyzer(kmersMap, false);
		int modeDepth = kmersAnalyzer.getMode();
		long expectedAssemblyLength = kmersAnalyzer.getExpectedAssemblyLength();
		
		if(compressionFactors!=null) {
			double averageCompression = NumberArrays.getAverage(compressionFactors);
			expectedAssemblyLength/= averageCompression;
		}
		log.info("Mode: "+modeDepth+" Expected assembly length: "+expectedAssemblyLength);
		
		AssemblyGraph graph = new AssemblyGraph(sequences);
		log.info("Created graph vertices. Edges: "+graph.getEdges().size());
		graph.setExpectedAssemblyLength(expectedAssemblyLength);
		graph.setPloidy(ploidy);
		
		long time1 = System.currentTimeMillis();
		ShortKmerCodesTable table = new ShortKmerCodesTable(kmersAnalyzer, kmerLength, windowLength);
		//MinimizersTable table = new MinimizersTable(kmersAnalyzer, kmerLength, windowLength);
		table.setLog(log);
		//table.setMaxAbundanceMinimizer(Math.max(100, 5*modeDepth));
		ThreadPoolExecutor poolMinimizers1 = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		int seqIdMinimizers = 0;
		long limit = 10*ploidy*expectedAssemblyLength;
		long totalLengthMinimizers = 0;
		while( seqIdMinimizers < sequences.size() ) {
			QualifiedSequence qseq = sequences.get(seqIdMinimizers);
			totalLengthMinimizers+=qseq.getLength();
			if(totalLengthMinimizers>limit) break;
			CharSequence seq = qseq.getCharacters();
			if(numThreads==1) {
				addSequenceToTable(table, seqIdMinimizers, seq);
			} else {
				final int i = seqIdMinimizers;
				poolMinimizers1.execute(()->addSequenceToTable(table, i, seq));
			}
			seqIdMinimizers++;
		}
		waitToFinish(sequences.size(), poolMinimizers1);
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time2 = System.currentTimeMillis();
		long diff = (time2-time1)/1000;
		log.info("Built minimizers for the first 10x of sequences. Time minimizers (s): "+diff+" Memory (Gbp): "+usedMemory+" first sequence search: "+seqIdMinimizers);
		//Distribution minimizerHitsDist = table.calculateDistributionHits();
		//minimizerHitsDist.printDistributionInt(System.out);
		KmerHitsAssemblyEdgesFinder edgesFinder = new KmerHitsAssemblyEdgesFinder(graph);
		
		List<List<AssemblySequencesRelationship>> relationshipsPerSequence = new ArrayList<List<AssemblySequencesRelationship>>(sequences.size());
		for(int i=0;i<sequences.size();i++) relationshipsPerSequence.add(null);
		
		//Find first edges between the longest reads
		ThreadPoolExecutor poolSearch = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		edgesFinder.setExtensiveSearch(true);
		for (int seqId = 0; seqId < seqIdMinimizers; seqId++) {
			CharSequence seq = sequences.get(seqId).getCharacters();
			double compressionFactor = compressionFactors!=null?compressionFactors[seqId]:1;
			final int i = seqId;
			poolSearch.execute(()->processSequence(edgesFinder, table, i, seq, compressionFactor, false, relationshipsPerSequence));
		}
		waitToFinish(sequences.size(), poolSearch);
		
		//Find embedded relationships in not indexed reads
		poolSearch = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		edgesFinder.setExtensiveSearch(false);
		for (int seqId = seqIdMinimizers; seqId < sequences.size(); seqId++) {
			CharSequence seq = sequences.get(seqId).getCharacters();
			double compressionFactor = compressionFactors!=null?compressionFactors[seqId]:1;
			final int i = seqId;
			poolSearch.execute(()->processSequence(edgesFinder, table, i, seq, compressionFactor, true, relationshipsPerSequence));
		}
		waitToFinish(sequences.size(), poolSearch);
		/*int countCurrentEmbedded = 0;
		for(int i=0;i<relationshipsPerSequence.size();i++) {
			if(relationshipsPerSequence.get(i)!=null) countCurrentEmbedded++;
		}*/
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time3 = System.currentTimeMillis();
		diff = (time3-time2)/1000;
		log.info("Identified relationships with the first 10x of sequences. Time mapping (s): "+diff+". Memory: "+usedMemory);
		
		//Index non embedded reads
		ThreadPoolExecutor poolMinimizers2 = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for(int seqId = seqIdMinimizers ;seqId < sequences.size();seqId++ ) {
			//if (seqIdMinimizers%1000==0) log.info("Seqid: "+seqIdMinimizers+" Current list: "+relationshipsPerSequence.get(seqIdMinimizers));
			if(relationshipsPerSequence.get(seqId)!=null) continue;
			QualifiedSequence qseq = sequences.get(seqId);
			CharSequence seq = qseq.getCharacters();
			if(numThreads==1) {
				addSequenceToTable(table, seqId, seq);
			} else {
				final int i = seqId;
				poolMinimizers2.execute(()->addSequenceToTable(table, i, seq));
			}
		}
		waitToFinish(sequences.size(), poolMinimizers2);
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time4 = System.currentTimeMillis();
		diff = (time4-time3)/1000;
		log.info("Built minimizers for the remaining non embedded sequences. Time minimizers (s): "+diff+" Memory (Gbp): "+usedMemory);
		edgesFinder.setCompleteAlignment(false);
		edgesFinder.setExtensiveSearch(true);
		
		ThreadPoolExecutor poolSearch2 = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for (int seqId = seqIdMinimizers; seqId < sequences.size(); seqId++) {
			CharSequence seq = sequences.get(seqId).getCharacters();
			double compressionFactor = compressionFactors!=null?compressionFactors[seqId]:1;
			final int i = seqId;
			poolSearch2.execute(()->processSequence(edgesFinder, table, i, seq, compressionFactor, false, relationshipsPerSequence));
			//if ((seqId+1)%1000==0) log.info("Scheduled sequence "+(seqId+1));
		}
		addRelationshipsToGraph(graph, relationshipsPerSequence, runtime);
		waitToFinish(sequences.size(), poolSearch2);
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time5 = System.currentTimeMillis();
		diff = (time5-time4)/1000;
		log.info("Built graph. Edges: "+graph.getEdges().size()+" Embedded: "+graph.getEmbeddedCount()+" Memory: "+usedMemory+" Time graph construction (s): "+diff);
		Set<Integer> orphanEmbeddedIds = graph.calculateEmbeddedToChimeric();
		if(orphanEmbeddedIds.size()>0) {
			ShortKmerCodesTable tableOrphans = new ShortKmerCodesTable(kmersAnalyzer, kmerLength, windowLength);
			for(int seqId:orphanEmbeddedIds) {
				for(AssemblySequencesRelationship rel:relationshipsPerSequence.get(seqId)) graph.removeRelationship(rel);
				relationshipsPerSequence.set(seqId, null);
				CharSequence seq = sequences.get(seqId).getCharacters();
				addSequenceToTable(tableOrphans, seqId, seq);
				double compressionFactor = compressionFactors!=null?compressionFactors[seqId]:1;
				final int i = seqId;
				processSequence(edgesFinder, table, i, seq, compressionFactor, false, relationshipsPerSequence);
				for(AssemblySequencesRelationship rel:relationshipsPerSequence.get(seqId)) graph.addRelationship(rel);
			}
			List<List<AssemblySequencesRelationship>> orphanRelationshipsPerSequence = new ArrayList<List<AssemblySequencesRelationship>>(sequences.size());
			for(int i=0;i<sequences.size();i++) orphanRelationshipsPerSequence.add(null);
			for(int seqId = seqIdMinimizers; seqId < sequences.size(); seqId++) {
				if(graph.isEmbedded(seqId)) continue;
				CharSequence seq = sequences.get(seqId).getCharacters();
				double compressionFactor = compressionFactors!=null?compressionFactors[seqId]:1;
				final int i = seqId;
				processSequence(edgesFinder, tableOrphans, i, seq, compressionFactor, false, orphanRelationshipsPerSequence);
				for(AssemblySequencesRelationship rel:orphanRelationshipsPerSequence.get(seqId)) graph.addRelationship(rel);
			}
			usedMemory = runtime.totalMemory()-runtime.freeMemory();
			usedMemory/=1000000000;
			long time6 = System.currentTimeMillis();
			diff = (time6-time5)/1000;
			log.info("Recalculated information for "+orphanEmbeddedIds.size()+" sequences embedded into chimeric sequences. Edges: "+graph.getEdges().size()+" Embedded: "+graph.getEmbeddedCount()+" Memory: "+usedMemory+" Time reprocessing (s): "+diff);
		}
		
		//log.info(" Raw hits for "+edgesFinder.getCountRawHits()+" sequences. Completed hits for "+edgesFinder.getCountCompletedHits()+" sequences");
		return graph;
	}
	
	public void countSequenceKmers(KmersExtractor extractor, int seqId, QualifiedSequence seq) {
		extractor.countSequenceKmers(seq);
		if (seqId%1000==0) log.info("Kmers extracted for "+(seqId)+" sequences.");
	}
	//public void addSequenceToTable(MinimizersTable table, int seqId, CharSequence seq) {
	public void addSequenceToTable(ShortKmerCodesTable table, int seqId, CharSequence seq) {
		table.addSequence(seqId, seq);
		if(seqId == idxDebug) System.out.println("Added sequence: "+seqId);
		if (seqId%1000==0) log.info("Processed "+(seqId)+" sequences. Total minimizers: "+table.size()+" total entries: "+table.getTotalEntries());
	}
	//private void processSequence(KmerHitsAssemblyEdgesFinder finder, MinimizersTable table, int seqId, CharSequence seq, double compressionFactor, boolean onlyEmbedded, List<List<AssemblySequencesRelationship>> relationshipsPerSequence ) {
	private void processSequence(KmerHitsAssemblyEdgesFinder finder, ShortKmerCodesTable table, int seqId, CharSequence seq, double compressionFactor, boolean onlyEmbedded, List<List<AssemblySequencesRelationship>> relationshipsPerSequence ) {
		try {
			List<AssemblySequencesRelationship> rels = relationshipsPerSequence.get(seqId);  
			if(rels==null) {
				Map<Integer,List<UngappedSearchHit>> hitsForward = table.match(seqId, seq);
				String complement = DNAMaskedSequence.getReverseComplement(seq).toString();
				Map<Integer,List<UngappedSearchHit>> hitsReverse = table.match(seqId, complement);
				rels = finder.inferRelationshipsFromKmerHits(seqId, seq.toString(), complement, hitsForward, hitsReverse, compressionFactor);
				//if(seqId == 19584) System.out.println("Total relationships identified for sequence "+seqId+" "+rels.size()+" onlyEmbedded: "+onlyEmbedded);
				//rels = new ArrayList<AssemblySequencesRelationship>();
				if(!onlyEmbedded) relationshipsPerSequence.set(seqId, rels);
				else {
					rels = selectGoodEmbedded(rels);
					if(rels.size()>=ploidy) relationshipsPerSequence.set(seqId, rels);
				}
			}
			if ((seqId)%1000==0) {
				int edges = 0;
				int embedded = 0;
				for(AssemblySequencesRelationship next:rels) {
					if(next instanceof AssemblyEmbedded) embedded++;
					if(next instanceof AssemblyEdge) edges++;
				}
				log.info("Identified relationships for sequence "+(seqId) +" Candidate edges: "+edges+"  candidate embedded hosts "+embedded);
				//if (onlyEmbedded) log.info("List: "+relationshipsPerSequence.get(seqId));
			}
		} catch (RuntimeException e) {
			if(!onlyEmbedded) relationshipsPerSequence.set(seqId, new ArrayList<AssemblySequencesRelationship>());
			throw e;
		}
	}
	private List<AssemblySequencesRelationship> selectGoodEmbedded(List<AssemblySequencesRelationship> rels) {
		List<AssemblySequencesRelationship> answer = new ArrayList<AssemblySequencesRelationship>();
		for(AssemblySequencesRelationship rel:rels) {
			if(rel instanceof AssemblyEmbedded) {
				AssemblyEmbedded embedded = (AssemblyEmbedded)rel;
				//if(embedded.getAlignment()!=null ) return true;
				if(embedded.getEvidenceProportion()>0.99 && embedded.getIndelsPerKbp()<10 && embedded.getWeightedCoverageSharedKmers()>0.5*embedded.getRead().getLength()) answer.add(embedded);
			}
		}
		return answer;
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
				for(AssemblySequencesRelationship next:nextList) graph.addRelationship(next);
				//nextList.clear();
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