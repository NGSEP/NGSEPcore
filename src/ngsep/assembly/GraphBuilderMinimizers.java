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
package ngsep.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmerSearchResultsCompressedTable;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMapAnalyzer;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.ShortKmerCodesSampler;
import ngsep.sequences.ShortKmerCodesTable;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class GraphBuilderMinimizers implements GraphBuilder {

	private Logger log = Logger.getLogger(GraphBuilderMinimizers.class.getName());
	
	public static final int DEF_WINDOW_LENGTH = 10;
	public static final int DEF_NUM_THREADS = 1;
	
	
	private int ploidy = AssemblyGraph.DEF_PLOIDY_ASSEMBLY;
	private int numThreads = DEF_NUM_THREADS;
	private int depthCompleteIndexing = 10;
	private KmersMapAnalyzer kmersAnalyzer;
	private ShortKmerCodesSampler sampler; 
	
	
	private static final int TIMEOUT_SECONDS = 30;
	
	private static int idxDebug = -1;
	
	public GraphBuilderMinimizers(KmersMapAnalyzer kmersAnalyzer, ShortKmerCodesSampler sampler) {
		this.kmersAnalyzer = kmersAnalyzer;
		this.sampler = sampler;
	}
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
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
	
	public KmersMapAnalyzer getKmersAnalyzer() {
		return kmersAnalyzer;
	}
	public void setKmersAnalyzer(KmersMapAnalyzer kmersAnalyzer) {
		this.kmersAnalyzer = kmersAnalyzer;
	}
	@Override	
	public AssemblyGraph buildAssemblyGraph(final List<QualifiedSequence> sequences) {
		Runtime runtime = Runtime.getRuntime();
		int modeDepth = Math.max(1,kmersAnalyzer.getMode());
		long expectedAssemblyLength = kmersAnalyzer.getExpectedAssemblyLength();
		
		log.info("Mode: "+modeDepth+" Expected assembly length: "+expectedAssemblyLength);
		
		AssemblyGraph graph = new AssemblyGraph(sequences);
		log.info("Created graph vertices. Edges: "+graph.getEdges().size());
		graph.setExpectedAssemblyLength(expectedAssemblyLength);
		graph.setPloidy(ploidy);
		
		long time1 = System.currentTimeMillis();
		
		
		//Create the structures with appropriate initial capacity
		int capacity = kmersAnalyzer.getKmersMap().size()/10;
		ShortKmerCodesTable table = new ShortKmerCodesTable(sampler,capacity,false);
		table.setLog(log);
		table.setMode(modeDepth);
		table.setKmerDistModeLocalSD(kmersAnalyzer.getModeLocalSD());
		
		
		
		//table.setMaxAbundanceMinimizer(Math.max(100, 5*modeDepth));
		ThreadPoolExecutor poolMinimizers1 = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		int seqIdMinimizers = 0;
		long limit = depthCompleteIndexing*ploidy*expectedAssemblyLength;
		long totalLengthMinimizers = 0;
		int n = sequences.size();
		while( seqIdMinimizers < n) {
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
		waitToFinish(n, poolMinimizers1);
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time2 = System.currentTimeMillis();
		long diff = (time2-time1)/1000;
		log.info("Built minimizers for the first "+depthCompleteIndexing+"x of sequences. Time minimizers (s): "+diff+" Memory (Gbp): "+usedMemory+" first sequence not indexed: "+seqIdMinimizers);
		//Distribution minimizerHitsDist = table.calculateDistributionHits();
		//minimizerHitsDist.printDistributionInt(System.out);
		KmerHitsAssemblyEdgesFinder edgesFinder = new KmerHitsAssemblyEdgesFinder(graph);
		
		List<List<AssemblySequencesRelationship>> relationshipsPerSequence = new ArrayList<List<AssemblySequencesRelationship>>(sequences.size());
		for(int i=0;i<n;i++) relationshipsPerSequence.add(null);
		
		//Find first edges between the longest reads
		ThreadPoolExecutor poolSearch = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		edgesFinder.setExtensiveSearch(true);
		int firstSeqNonExtensiveSearch = seqIdMinimizers;
		//int firstSeqNonExtensiveSearch = 2*seqIdMinimizers;
		for (int seqId = 0; seqId < firstSeqNonExtensiveSearch; seqId++) {
			CharSequence seq = sequences.get(seqId).getCharacters();
			final int i = seqId;
			poolSearch.execute(()->processSequence(edgesFinder, table, i, seq, false, relationshipsPerSequence));
		}
		waitToFinish(n, poolSearch);
		
		//Find embedded relationships in not indexed reads
		poolSearch = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		edgesFinder.setExtensiveSearch(false);
		for (int seqId = firstSeqNonExtensiveSearch; seqId < n; seqId++) {
			CharSequence seq = sequences.get(seqId).getCharacters();
			final int i = seqId;
			poolSearch.execute(()->processSequence(edgesFinder, table, i, seq, true, relationshipsPerSequence));
		}
		waitToFinish(n, poolSearch);
		boolean [] added = new boolean[n];
		Arrays.fill(added, false);
		addRelationshipsToGraph(graph, relationshipsPerSequence, 0, added, false, runtime);
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time3 = System.currentTimeMillis();
		diff = (time3-time2)/1000;
		log.info("Identified relationships with the first "+depthCompleteIndexing+"x of sequences. Repetitive: "+repetitiveSequenceIds.size()+" Time mapping (s): "+diff+". Memory: "+usedMemory);
		
		//Remove embedded relationships for embedded to chimeric sequences
		Set<Integer> orphanEmbeddedIds = graph.calculateEmbeddedToChimeric();
		log.info("Identified "+orphanEmbeddedIds.size()+" sequences embedded to chimeric sequences");
		for(int seqId:orphanEmbeddedIds) {
			if(seqId<seqIdMinimizers) continue;
			List<AssemblyEmbedded> embeddedRels = graph.getEmbeddedBySequenceId(seqId);
			if(seqId == idxDebug) System.out.println("Processing orphan sequence: "+seqId+" embedded: "+graph.isEmbedded(seqId)+" current rels: "+relationshipsPerSequence.get(seqId)+" embedded rels graph: "+embeddedRels);
			for(AssemblySequencesRelationship rel:relationshipsPerSequence.get(seqId)) graph.removeRelationship(rel);
			relationshipsPerSequence.set(seqId, null);
			added[seqId] = false;
		}
		
		//Index non embedded reads
		ThreadPoolExecutor poolMinimizers2 = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for(int seqId = seqIdMinimizers ;seqId < sequences.size();seqId++ ) {
			List<AssemblyEmbedded> embeddedRels = graph.getEmbeddedBySequenceId(seqId);
			if (seqId == idxDebug) log.info("Seqid: "+seqId+" Current list: "+relationshipsPerSequence.get(seqId)+" embedded rels: "+embeddedRels);
			if(embeddedRels.size()>0) continue;
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
		
		//Search again non embedded reads
		ThreadPoolExecutor poolSearch2 = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		for (int seqId = seqIdMinimizers; seqId < sequences.size(); seqId++) {
			List<AssemblyEmbedded> embeddedRels = graph.getEmbeddedBySequenceId(seqId);
			if (seqId == idxDebug) log.info("Seqid: "+seqId+" Current list: "+relationshipsPerSequence.get(seqId)+" embedded rels: "+embeddedRels);
			if(embeddedRels.size()>0) continue;
			List<AssemblySequencesRelationship> currentRels = relationshipsPerSequence.get(seqId);
			if(currentRels!=null) {
				for(AssemblySequencesRelationship rel:currentRels) graph.removeRelationship(rel);
				relationshipsPerSequence.set(seqId, null);
			}
			CharSequence seq = sequences.get(seqId).getCharacters();
			final int i = seqId;
			poolSearch2.execute(()->processSequence(edgesFinder, table, i, seq, false, relationshipsPerSequence));
			//if ((seqId+1)%1000==0) log.info("Scheduled sequence "+(seqId+1));
		}
		addRelationshipsToGraph(graph, relationshipsPerSequence, seqIdMinimizers, added, true, runtime);
		waitToFinish(sequences.size(), poolSearch2);
		
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time5 = System.currentTimeMillis();
		diff = (time5-time4)/1000;
		log.info("Built graph. Edges: "+graph.getEdges().size()+" Embedded: "+graph.getEmbeddedCount()+"  Repetitive: "+repetitiveSequenceIds.size()+" Memory: "+usedMemory+" Time graph construction (s): "+diff);
		//updateRelationshipsRepeats(graph,repetitiveSequenceIds);
		log.info("Finished graph construction");
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
	private int totalRels = 0;
	private Set<Integer> repetitiveSequenceIds = new HashSet<>();
	private void processSequence(KmerHitsAssemblyEdgesFinder finder, ShortKmerCodesTable table, int seqId, CharSequence seq, boolean onlyEmbedded, List<List<AssemblySequencesRelationship>> relationshipsPerSequence ) {
		try {
			List<AssemblySequencesRelationship> rels = relationshipsPerSequence.get(seqId);
			boolean debug = seqId == idxDebug; 
			if(debug) System.out.println("Identifying relationships for sequence "+seqId+" current: "+rels);
			if(rels==null) {
				long time0 = System.currentTimeMillis();
				if(debug) log.info("GraphBuilderMinimizers. Sequence "+seqId+" total rels: "+totalRels+" Memory: "+KmerHitsAssemblyEdgesFinder.calculateMemoryGbp());
				KmerSearchResultsCompressedTable hitsForward = table.matchCompressed(seqId, seq, seqId);
				long time1 = System.currentTimeMillis();
				long diff1 = (time1-time0)/1000;
				debug = debug || diff1>10;
				//if(debug) log.info("GraphBuilderMinimizers. Sequence "+seqId+" total rels: "+totalRels+" forward: "+hitsForward.getTotalHits()+" time1: "+diff1+" Memory: "+KmerHitsAssemblyEdgesFinder.calculateMemoryGbp());
				String complement = DNAMaskedSequence.getReverseComplement(seq).toString();
				KmerSearchResultsCompressedTable hitsReverse = table.matchCompressed(seqId, complement, seqId);
				long time2 = System.currentTimeMillis();
				long diff2 = (time2-time1)/1000;
				debug = debug || diff2>10;
				//if(debug) log.info("GraphBuilderMinimizers. Sequence "+seqId+" total rels: "+totalRels+" forward: "+hitsForward.getTotalHits()+" reverse: "+hitsReverse.getTotalHits()+" times12: "+diff1+" "+diff2+" Memory: "+KmerHitsAssemblyEdgesFinder.calculateMemoryGbp());
				rels = finder.inferRelationshipsFromKmerHits(seqId, seq.toString(), complement, hitsForward, hitsReverse);
				long time3 = System.currentTimeMillis();
				long diff3 = (time3-time2)/1000;
				debug = debug || diff3>10;
				int repStatus = calculateRepetitiveStatus(seqId, hitsForward); 
				debug = debug || repStatus>3;
				if(debug) log.info("GraphBuilderMinimizers. Rels for seq "+seqId+" "+rels.size()+" len: "+seq.length()+" totalRels: "+totalRels+" onlyEmbedded: "+onlyEmbedded+" repStatus: " +repStatus+" InternalRepKmers "+hitsForward.getInternalMultihitUsedKmerCodes().size()+" used codes "+hitsForward.getUsedCodesCount()+" selfHits: "+hitsForward.getKmerHitCount(seqId)+" countsF: "+hitsForward.getInputCodesCount()+" "+hitsForward.getDistinctCodesCount()+" "+ hitsForward.getTotalHits()+" " +hitsForward.getMultisequenceCodesCount()+" "+hitsForward.getMultihitCodesCount()+" "+hitsForward.getNotFoundCodesCount()+" "+hitsForward.getDistinctUsedCodesCount()+" countsR: "+hitsReverse.getInputCodesCount()+" "+hitsReverse.getDistinctCodesCount()+" "+hitsReverse.getTotalHits()+" "+ hitsReverse.getMultisequenceCodesCount()+" "+hitsReverse.getNotFoundCodesCount()+" subjects: "+hitsForward.getNumSubjects()+" "+hitsReverse.getNumSubjects()+" timesAll: "+diff1+" "+diff2+" "+diff3+" Memory: "+KmerHitsAssemblyEdgesFinder.calculateMemoryGbp());
				if(debug) {
					int [] countsF = hitsForward.getNormalizedCountsDist();
					int [] countsR = hitsReverse.getNormalizedCountsDist();
					System.out.println("Sequence "+seqId+" Profile Forward: "+printArray(countsF));
					System.out.println("Sequence "+seqId+" Profile Reverse: "+printArray(countsR));
				}
				//rels = new ArrayList<AssemblySequencesRelationship>();
				if(!onlyEmbedded) relationshipsPerSequence.set(seqId, rels);
				else {
					rels = selectGoodEmbedded(rels);
					if(rels.size()>=1) relationshipsPerSequence.set(seqId, rels);
				}
				synchronized (this) {
					totalRels+=rels.size();
					if(repStatus>0) repetitiveSequenceIds.add(seqId);
				}
				
			}
			if ((seqId)%1000==0) {
				int edges = 0;
				int embedded = 0;
				for(AssemblySequencesRelationship next:rels) {
					if(next instanceof AssemblyEmbedded) embedded++;
					if(next instanceof AssemblyEdge) edges++;
				}
				log.info("Identified relationships for sequence "+(seqId) +" Candidate edges: "+edges+"  candidate embedded hosts "+embedded+" Memory: "+KmerHitsAssemblyEdgesFinder.calculateMemoryGbp());
				//if (onlyEmbedded) log.info("List: "+relationshipsPerSequence.get(seqId));
			}
		} catch (RuntimeException e) {
			if(!onlyEmbedded) relationshipsPerSequence.set(seqId, new ArrayList<AssemblySequencesRelationship>());
			throw e;
		}
	}
	private int calculateRepetitiveStatus(int seqId, KmerSearchResultsCompressedTable hits) {
		int answer = 0;
		Set<Long> internalRepetitiveKmerCodes = hits.getInternalMultihitUsedKmerCodes();
		//for(Map.Entry<Long,Set<Integer>> entry:internalRepetitiveKmers.entrySet()) {
			//System.err.println("Seq id: "+seqId+" nextRepKmer: "+entry.getKey()+" times: "+entry.getValue().size()+" positions: "+entry.getValue());
		//}
		if(internalRepetitiveKmerCodes.size()>0.05*hits.getDistinctUsedCodesCount()) answer = 1;
		int [] profileRep = hits.getNormalizedCountsDist();
		int sum = 0;
		for(int i=2;i<profileRep.length;i++) sum+=profileRep[i];
		if(profileRep[1]+sum>0.6*profileRep[0] && sum>0.05*profileRep[0]) answer += 2;
		  
		return answer;
	}
	private String printArray(int[] normalizedCountsDist) {
		StringBuilder answer = new StringBuilder();
		for(int i=0;i<normalizedCountsDist.length;i++) answer.append(" "+normalizedCountsDist[i]);
		
		return answer.toString();
	}
	private List<AssemblySequencesRelationship> selectGoodEmbedded(List<AssemblySequencesRelationship> rels) {
		List<AssemblySequencesRelationship> answer = new ArrayList<AssemblySequencesRelationship>();
		for(AssemblySequencesRelationship rel:rels) {
			if(rel instanceof AssemblyEmbedded) {
				AssemblyEmbedded embedded = (AssemblyEmbedded)rel;
				//if(embedded.getAlignment()!=null ) return true;
				if(KmerHitsAssemblyEdgesFinder.isGoodEmbedded(embedded)) answer.add(embedded);
			}
		}
		return answer;
	}
	private void addRelationshipsToGraph(AssemblyGraph graph, List<List<AssemblySequencesRelationship>> relationshipsPerSequence, int firstIndex, boolean [] added, boolean waitForNull, Runtime runtime) {
		int n = relationshipsPerSequence.size();
		int i=firstIndex;
		log.info("Adding relationships to graph");
		while(i<n) {
			if(waitForNull) {
				try {
					Thread.sleep(1);
				} catch (InterruptedException e) {
					throw new RuntimeException("Add relationships thread interrupted", e);
				}
			}
			
			for(;i<n;i++) {
				if(added[i]) continue;
				List<AssemblySequencesRelationship> nextList = relationshipsPerSequence.get(i);
				if(nextList == null) {
					if(waitForNull) break;
					else continue;
				}
				//if ((i+1)%1000==0) log.info("Adding relationships for sequence "+(i+1) +" Relationships sequence: "+nextList.size());
				for(AssemblySequencesRelationship next:nextList) graph.addRelationship(next);
				added[i] = true;
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
	//private void updateRelationshipsRepeats(AssemblyGraph graph, Set<Integer> repetitiveSequenceIds) {
		//List<Set<Integer>> connectedComponents = graph.findConnectedComponentsSequences(repetitiveSequenceIds);
		//log.info("RemoveEdgesRepeats. Num connected components: "+connectedComponents.size());
		/*KmerHitsAssemblyEdgesFinder edgesFinder = new KmerHitsAssemblyEdgesFinder(graph);
		for(Set<Integer> cluster: connectedComponents) {
			log.info("RemoveEdgesRepeats. Next connected component. Size: "+cluster.size()+" Elements: "+cluster);
			if(cluster.size()<100) continue;
			List<Integer> clusterList = new ArrayList<>(cluster);
			Collections.sort(clusterList);
			KmersMap localMap = calculateKmersMap(graph, cluster);
			//Distribution kmersDist = localMap.calculateAbundancesDistribution();
			KmersMapAnalyzer kmersAnalyzer = new KmersMapAnalyzer(localMap, false);
			System.out.println("Mode depth: "+kmersAnalyzer.getMode());
			ShortKmerCodesTable table = new ShortKmerCodesTable(kmersAnalyzer, kmerLength, windowLength);
			for(int i:clusterList) {
				QualifiedSequence seq = graph.getSequence(i);
				System.out.println("Calculating local minimizers for local repetitive sequence: "+i);
				table.addSequence(i, seq.getCharacters());
			}
			int n = graph.getNumSequences();
			List<List<AssemblySequencesRelationship>> relationshipsPerSequence = new ArrayList<List<AssemblySequencesRelationship>>(n);
			for(int i=0;i<n;i++) relationshipsPerSequence.add(null);
			for(int i:clusterList) {
				System.out.println("Updating relationships for repetitive sequence: "+i);
				QualifiedSequence seq = graph.getSequence(i);
				processSequence(edgesFinder, table, i, seq.getCharacters(), false, relationshipsPerSequence);
			}
			
			AssemblyGraph subgraph = graph.buildSubgraph(cluster);
			List<AssemblyEdge> edges = subgraph.getEdges();
			List<AssemblyEmbedded> embeddedRels = subgraph.getAllEmbedded();
			log.info("RemoveEdgesRepeats. Next connected component subgraph. Edges: "+edges.size()+" embedded: "+embeddedRels.size());
			for(AssemblyEdge edge:edges) {
				if(!edge.isSameSequenceEdge()) graph.removeEdge(edge);
			}
			for(AssemblyEmbedded embedded: embeddedRels) graph.removeEmbedded(embedded);
			for(int i:clusterList) {
				for(AssemblySequencesRelationship rel:relationshipsPerSequence.get(i)) {
					graph.addRelationship(rel);
				}
			}*/
			/*Distribution ikbpDist = new Distribution(0, 10, 0.1);
			for(AssemblyEdge edge:edges) {
				//System.out.println("Next edge "+edge);
				ikbpDist.processDatapoint(edge.getIndelsPerKbp());
			}
			for(AssemblyEmbedded embedded:embeddedRels) {
				//System.out.println("Next embedded "+embedded);
				ikbpDist.processDatapoint(embedded.getIndelsPerKbp());
				//if(embedded.getIndelsPerKbp()>3) graph.removeEmbedded(embedded); 
			}
			System.out.println("Dist IKBP");
			ikbpDist.printDistribution(System.out);*/
		//}	
	//}
	/*private KmersMap calculateKmersMap(AssemblyGraph graph, Set<Integer> cluster) {
		KmersExtractor extractor = new KmersExtractor();
		extractor.setKmerLength(kmerLength);
		for(int i:cluster) {
			QualifiedSequence seq = graph.getSequence(i);
			extractor.countSequenceKmers(seq);
		}
		return extractor.getKmersMap();
	}*/
}