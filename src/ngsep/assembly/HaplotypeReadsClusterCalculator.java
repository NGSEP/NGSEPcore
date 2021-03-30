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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileWriter;
import ngsep.discovery.AlignmentsPileupGenerator;
import ngsep.discovery.PileupListener;
import ngsep.discovery.PileupRecord;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.haplotyping.SingleIndividualHaplotyper;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.SNV;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class HaplotypeReadsClusterCalculator {

	private Logger log = Logger.getLogger(HaplotypeReadsClusterCalculator.class.getName());
	public static final int DEF_NUM_THREADS = 1;
	private static final int TIMEOUT_SECONDS = 30;
	
	private int numThreads = DEF_NUM_THREADS;
	
	private int debugIdx = -1;
	
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	public int getNumThreads() {
		return numThreads;
	}

	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}

	public List<Set<Integer>> clusterReads(AssemblyGraph graph, int ploidy) {
		//System.out.println("HaplotypeReadsClustering. Edges printed vertex "+graph.getEdges(graph.getVertex(181, true)).size());
		ThreadPoolExecutor poolClustering = new ThreadPoolExecutor(numThreads, numThreads, TIMEOUT_SECONDS, TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		List<ClusterReadsTask> tasksList = new ArrayList<ClusterReadsTask>();
		List<List<AssemblyEdge>> paths = graph.getPaths();
		for(int i = 0; i < paths.size(); i++)
		{
			List<AssemblyEdge> path = paths.get(i);
			ClusterReadsTask task = new ClusterReadsTask(this, graph, path, i+1, ploidy);
			poolClustering.execute(task);
			tasksList.add(task);	
		}
		int finishTime = 10*graph.getNumSequences();
		poolClustering.shutdown();
		try {
			poolClustering.awaitTermination(finishTime, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
    	if(!poolClustering.isShutdown()) {
			throw new RuntimeException("The ThreadPoolExecutor was not shutdown after an await Termination call");
		}
    	return mergePathClusters(graph,tasksList,ploidy);
	}
	
	private List<Set<Integer>> mergePathClusters(AssemblyGraph graph, List<ClusterReadsTask> tasksList, int ploidy) {
		List<List<AssemblyEdge>> paths = graph.getPaths();
		log.info("Merging reads from "+paths.size()+" diploid paths");
		List<Set<Integer>> inputClusters = new ArrayList<Set<Integer>>();
    	//Reads of each cluster
		Map<Integer,Integer> readsClusters = new HashMap<Integer, Integer>();
		List<AssemblyVertex> verticesPaths = new ArrayList<AssemblyVertex>();
    	//Clusters that can not go in the same supercluster
    	Map<Integer,List<Integer>> clusterRestrictions = new HashMap<Integer, List<Integer>>();
    	int inputClusterId = 0;
    	for(int i=0;i<paths.size();i++) {
    		
    		int firstIdCluster = inputClusterId;
    		List<AssemblyEdge> path = paths.get(i);
    		ClusterReadsTask task = tasksList.get(i);
    		int pathId = task.getPathIdx();
    		
    		List<PhasedPathBlock> haplotypeBlocks = task.getPhasedBlocks();
    		log.info("Calculated "+haplotypeBlocks.size()+" blocks from haplotyping for path: "+pathId);
    		for(PhasedPathBlock block:haplotypeBlocks) {
    			int firstIdClusterBlock = inputClusterId;
    			List<Set<Integer>> inputClustersBlock = block.getPhasedReadIds();
    			System.out.println("Path: "+pathId+ " next block with "+inputClustersBlock.size() +" phased clusters");
    			for(Set<Integer> inputCluster:inputClustersBlock) {
    				for(int readId:inputCluster) { 
        				readsClusters.put(readId, inputClusterId);
        				//if (readId == 3692) System.out.println("Read "+readId+" "+graph.getSequence(readId).getName()+" input path: "+i+" input cluster: "+inputClusterId);
        				//if(inputCluster.size()<20 || readId%5==0) System.out.println("Cluster: "+inputClusterId+" read: "+readId+" "+graph.getSequence(readId).getName());
        				if(pathId==debugIdx) System.out.println("Cluster: "+inputClusterId+" read: "+readId+" "+graph.getSequence(readId).getName());
        			}
    				inputClusters.add(inputCluster);
    				clusterRestrictions.put(inputClusterId, new ArrayList<Integer>());
        			inputClusterId++;
    			}
    			int lastIdClusterBlock = inputClusterId-1;
    			
    			for(int j=firstIdClusterBlock;j<=lastIdClusterBlock;j++) {
    				for(int k=firstIdClusterBlock;k<=lastIdClusterBlock;k++) {
    					if(j!=k) {
    						clusterRestrictions.get(j).add(k);
    						clusterRestrictions.get(k).add(j);
    					}
    				}
    			}
    		}
    		
    		int lastIdCluster = inputClusterId-1;
    		System.out.println("Path: "+pathId+" first id cluster: "+firstIdCluster+" last id cluster: "+lastIdCluster);
    		
    		List<AssemblyVertex> verticesPath = extractVerticesPath(path, readsClusters);
    		System.out.println("Extracted "+verticesPath.size()+" vertices for path: "+pathId);
    		verticesPaths.addAll(verticesPath);
    	}
    	//Build connections graph
    	Map<String,ReadsClusterEdge> clusterEdgesMap = new HashMap<String, ReadsClusterEdge>();
    	for (AssemblyVertex vertex:verticesPaths) {
    		int readId = vertex.getSequenceIndex();
    		Integer clusterId = readsClusters.get(readId);
    		if(clusterId ==null) continue;
    		if(clusterId == 20) System.out.println("MergeHaplotypeClusters. Vertex: "+vertex+" cluster id: "+clusterId);
    		List<AssemblyEdge> edges = graph.getEdges(vertex);
    		for(AssemblyEdge edge:edges) {
    			if(edge.isSameSequenceEdge()) continue;
    			AssemblyVertex v2 = edge.getConnectingVertex(vertex);
    			Integer id2 = readsClusters.get(v2.getSequenceIndex());
    			if(clusterId == 20) System.out.println("MergeHaplotypeClusters. Connecting vertex: "+v2+" cluster id: "+id2+" restrictions: "+clusterRestrictions.get(clusterId));
    			if(id2!=null && clusterId!=id2 && !clusterRestrictions.get(clusterId).contains(id2) && ! clusterRestrictions.get(id2).contains(clusterId)) {
    				int idMin = Math.min(clusterId, id2);
    				int idMax = Math.max(clusterId, id2);
    				String key = ReadsClusterEdge.getKey(idMin,idMax);
    				ReadsClusterEdge clusterEdge = clusterEdgesMap.computeIfAbsent(key, (v)->new ReadsClusterEdge(idMin, idMax));
    				clusterEdge.addAssemblyEdge(edge);
    				if((idMin==13 || idMin == 14) && idMax == 20) System.out.println("Using edge "+edge+" for clusters joining. Current cluster edge:  "+clusterEdge);
    			}
    		}
    	}
    	
    	//Sort edges by total score
    	int n = clusterEdgesMap.size();
    	log.info("Created "+n+" edges between clusters");
    	List<ReadsClusterEdge> clusterEdgesList = new ArrayList<ReadsClusterEdge>(n);
    	clusterEdgesList.addAll(clusterEdgesMap.values());
    	Collections.sort(clusterEdgesList,(c1,c2)-> c2.getScore()-c1.getScore());
    	
    	//Perform clustering taking into account restrictions
    	Map<Integer,Integer> inputClustersAssignment = new HashMap<Integer,Integer>();
    	for(int i=0;i<n;i++) {
    		//Find next unused cluster
    		boolean change = false;
    		for(int j=i;j<n && !change;j++) {
    			ReadsClusterEdge edge = clusterEdgesList.get(j);
    			int c1 = edge.getClusterId1();
        		int c2 = edge.getClusterId2();
        		Integer assignment1 = inputClustersAssignment.get(c1);
        		Integer assignment2 = inputClustersAssignment.get(c2);
        		if(assignment1==null && assignment2==null) {
        			System.out.println("Joining clusters "+c1+" "+c2+ " score: "+edge.getScore()+" assignment: 0 rankings: "+i+" "+j);
        			assignCluster(inputClustersAssignment, c1, 0, clusterRestrictions.get(c1), ploidy);
        			assignCluster(inputClustersAssignment, c2, 0, clusterRestrictions.get(c2), ploidy);
        			change = true;
        		}
    		}
    		if(!change) break;
    		//Add clusters connected with already assigned clusters
    		while(change) {
    			change = false;
    			for(int j=i;j<n && !change;j++) {
    				ReadsClusterEdge edge = clusterEdgesList.get(j);
        			int c1 = edge.getClusterId1();
            		int c2 = edge.getClusterId2();
            		Integer assignment1 = inputClustersAssignment.get(c1);
            		Integer assignment2 = inputClustersAssignment.get(c2);
            		if(assignment1!=null && assignment2!=null) continue;
            		else if(assignment1==null && assignment2==null) continue;
            		else if (assignment1==null) {
        				System.out.println("Joining cluster "+c1+" with assigned cluster "+c2+ " score: "+edge.getScore()+" assignment: "+assignment2+" rankings: "+i+" "+j);
            			assignCluster(inputClustersAssignment, c1, assignment2, clusterRestrictions.get(c1), ploidy);
            			change = true;
            		} else {
        				System.out.println("Joining cluster "+c2+" with assigned cluster "+c1+ " score: "+edge.getScore()+" assignment: "+assignment1+" rankings: "+i+" "+j);
            			assignCluster(inputClustersAssignment, c2, assignment1, clusterRestrictions.get(c2), ploidy);
            			change = true;
            		}
    			}
    		}
    	}
    	log.info("Finished haplotype reads clustering. Assigned "+inputClustersAssignment.size()+" input clusters");
    	//Build answer from clusters
    	List<Set<Integer>> answer = new ArrayList<Set<Integer>>(ploidy);
    	for(int i=0;i<ploidy;i++) {
    		Set<Integer> hapSuperCluster = new HashSet<Integer>();
        	answer.add(hapSuperCluster);
    	}
    	for(Map.Entry<Integer, Integer> entry:inputClustersAssignment.entrySet()) {
    		Set<Integer> inputCluster = inputClusters.get(entry.getKey());
    		Set<Integer> outputCluster = answer.get(entry.getValue());
    		//System.out.println("MergeHaplotypeClusters. Input cluster "+entry.getKey()+" assigned to output cluster "+entry.getValue());
    		outputCluster.addAll(inputCluster); 
    	}
    	return answer;
	}
	
	private List<AssemblyVertex> extractVerticesPath(List<AssemblyEdge> path, Map<Integer,Integer> readsClusters) {
		List<AssemblyVertex> answer = new ArrayList<AssemblyVertex>(path.size()+1);
		if(path.size()==0) return answer;
		if(path.size()==1) {
			answer.add(path.get(0).getVertex1());
			answer.add(path.get(0).getVertex2());
			return answer;
		}
		for(AssemblyEdge edge:path) {
			if(edge.isSameSequenceEdge()) {
				answer.add(edge.getVertex1());
				answer.add(edge.getVertex2());
			}
		}
		//Number of vertices to select per cluster
		/*int n = 5;
		Map<Integer,Integer> countsPerCluster = new HashMap<Integer, Integer>();
		int m = path.size();
		AssemblyEdge edge0 = path.get(0);
		AssemblyEdge edge1 = path.get(1);
		AssemblyVertex vInternal = edge0.getSharedVertex(edge1);
		AssemblyVertex v0External = edge0.getConnectingVertex(vInternal);
		
		
		answer.add(v0External);
		int clusterId = readsClusters.getOrDefault(v0External.getSequenceIndex(),-1);
		if(clusterId>=0) countsPerCluster.compute(clusterId, (k,v)->(v==null)?1:v+1);
		int i=1;
		while(i<path.size()-1) {
			AssemblyEdge edge = path.get(i);
			AssemblyVertex vExt = edge.getConnectingVertex(vInternal);
			clusterId = readsClusters.getOrDefault(vExt.getSequenceIndex(),-1);
			if(clusterId>=0 && countsPerCluster.getOrDefault(clusterId, 0) < n) {
				answer.add(vExt);
				countsPerCluster.compute(clusterId, (k,v)->(v==null)?1:v+1);
			}
			
			i++;
			edge = path.get(i);
			vInternal = edge.getConnectingVertex(vExt);
			i++;
		}
		edge0 = path.get(m-1);
		edge1 = path.get(m-2);
		vInternal = edge0.getSharedVertex(edge1);
		v0External = edge0.getConnectingVertex(vInternal);
		answer.add(v0External);
		clusterId = readsClusters.getOrDefault(v0External.getSequenceIndex(),-1);
		if(clusterId>=0) countsPerCluster.compute(clusterId, (k,v)->(v==null)?1:v+1);
		i=m-2;
		while(i>1) {
			AssemblyEdge edge = path.get(i);
			AssemblyVertex vExt = edge.getConnectingVertex(vInternal);
			clusterId = readsClusters.getOrDefault(vExt.getSequenceIndex(),-1);
			if(clusterId>=0 && countsPerCluster.getOrDefault(clusterId, 0) < n) {
				answer.add(vExt);
				countsPerCluster.compute(clusterId, (k,v)->(v==null)?1:v+1);
			}
			i--;
			edge = path.get(i);
			vInternal = edge.getConnectingVertex(vExt);
			i--;
		}*/
		return answer;
	}

	private void assignCluster (Map<Integer,Integer> inputClustersAssignment, int clusterId, int assignment, List<Integer> restrictions, int ploidy) {
		//This should work fine with diploids 
		inputClustersAssignment.put(clusterId,assignment);
		
		if(restrictions==null) return;
		int i=0;
		for(int j=0;i<restrictions.size() && j<ploidy;j++) {
			if(j!=assignment) {
				int c3 = restrictions.get(i);
				System.out.println("Assigning by restrictions cluster "+c3+" to "+j);
				inputClustersAssignment.put(c3,j);
				i++;
			}
		}
	}
	
	
	List<PhasedPathBlock> clusterReadsPath(AssemblyGraph graph, List<AssemblyEdge> path, int pathIdx, int ploidy) {
		AssemblyPathReadsAligner aligner = new AssemblyPathReadsAligner();
		aligner.setLog(log);
		aligner.setAlignEmbedded(true);
		aligner.alignPathReads(graph, path, pathIdx);
		StringBuilder rawConsensus = aligner.getConsensus();
		List<ReadAlignment> alignments = aligner.getAlignedReads();
		Set<Integer> unalignedReadIds = aligner.getUnalignedReadIds();
		List<List<ReadAlignment>> clusters = null;
		String sequenceName = "diploidPath_"+pathIdx;
		int countHetSNVs = 0;
		if(alignments.size()>0) {
			
			for(ReadAlignment aln:alignments) aln.setSequenceName(sequenceName);
			Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
			List<CalledGenomicVariant> hetSNVs = findHeterozygousSNVs(rawConsensus, alignments, sequenceName);
			countHetSNVs = hetSNVs.size();
			if(pathIdx == debugIdx) savePathFiles("debug_"+pathIdx,sequenceName, rawConsensus.toString(),alignments,hetSNVs);
			if(hetSNVs.size()>0) {
				SingleIndividualHaplotyper sih = new SingleIndividualHaplotyper();
				sih.setAlgorithmName(SingleIndividualHaplotyper.ALGORITHM_NAME_REFHAP);
				try {
					clusters = sih.phaseSequenceVariants(sequenceName, hetSNVs, alignments);
				} catch (IOException e) {
					throw new RuntimeException (e);
				}
			}
		}
		List<PhasedPathBlock> answer = new ArrayList<PhasedPathBlock>();
		if(clusters == null) {
			//System.out.println("No clusters for path: "+pathIdx+". hetSNVs: "+countHetSNVs+" alignments: "+alignments.size());
			Set<Integer> sequenceIds = new HashSet<Integer>();
			for(ReadAlignment aln:alignments) sequenceIds.add(aln.getReadNumber());
			//Add not aligned reads within the path
			sequenceIds.addAll(unalignedReadIds);
			PhasedPathBlock block = new PhasedPathBlock();
			block.addPhasedReadIds(sequenceIds);
			answer.add(block);
			return answer;
		}
		log.info("Path: "+sequenceName+". hetSNVs: "+countHetSNVs+" alignments: "+alignments.size()+" clusters from haplotyping: "+clusters.size());
		Set<Integer> readIdsInPhasedBlocks = new HashSet<Integer>();
		List<GenomicRegion> phasedPathBoundaries = new ArrayList<GenomicRegion>(clusters.size()/2);
		Map<Integer,PhasedPathBlock> phasedBlocksByFirst = new HashMap<Integer, PhasedPathBlock>();
		for(int i=0;i<clusters.size();i+=2) {
			List<ReadAlignment> alnsHap0 = clusters.get(i);
			List<ReadAlignment> alnsHap1 = clusters.get(i+1);
			//Not real phasing
			if(alnsHap0.size()==0 || alnsHap1.size()==0) continue;
			ReadAlignment firstAln0 = alnsHap0.get(0);
			ReadAlignment firstAln1 = alnsHap1.get(0);
			GenomicRegionImpl region = new GenomicRegionImpl(firstAln0.getSequenceName(), Math.min(firstAln0.getFirst(), firstAln1.getFirst()), Math.max(firstAln0.getLast(), firstAln1.getLast()));
			Set<Integer> sequenceIdsHap0 = new HashSet<Integer>();
			for(ReadAlignment aln:alnsHap0) {
				sequenceIdsHap0.add(aln.getReadNumber());
				region.setFirst(Math.min(region.getFirst(), aln.getFirst()));
				region.setLast(Math.max(region.getLast(), aln.getLast()));
			}
			Set<Integer> sequenceIdsHap1 = new HashSet<Integer>();
			for(ReadAlignment aln:alnsHap1) {
				sequenceIdsHap1.add(aln.getReadNumber());
				region.setFirst(Math.min(region.getFirst(), aln.getFirst()));
				region.setLast(Math.max(region.getLast(), aln.getLast()));
			}
			PhasedPathBlock block = new PhasedPathBlock();
			block.addPhasedReadIds(sequenceIdsHap0);
			block.addPhasedReadIds(sequenceIdsHap1);
			answer.add(block);
			readIdsInPhasedBlocks.addAll(sequenceIdsHap0);
			readIdsInPhasedBlocks.addAll(sequenceIdsHap1);
			phasedBlocksByFirst.put(region.getFirst(), block);
			phasedPathBoundaries.add(region);
		}
		Collections.sort(phasedPathBoundaries, GenomicRegionPositionComparator.getInstance());
		int i=0;
		Set<Integer> nextUnphasedBlock = new HashSet<Integer>();
		for(ReadAlignment aln:alignments) {
			if(readIdsInPhasedBlocks.contains(aln.getReadNumber())) continue;
			boolean addRead = true;
			while(i<phasedPathBoundaries.size()) {
				GenomicRegion region = phasedPathBoundaries.get(i);
				if(region.getLast()>aln.getLast()) {
					if(aln.getFirst()>=region.getFirst()) {
						addRead = false;
						//Alignment contained in block but not assigned
						/*PhasedPathBlock block = phasedBlocksByFirst.get(region.getFirst());
						if(block!=null) {
							block.addHomozygousReadId(aln.getReadNumber());
						}*/
					}
					break;
				}
				if(nextUnphasedBlock.size()>0) {
					PhasedPathBlock block = new PhasedPathBlock();
					block.addPhasedReadIds(nextUnphasedBlock);
					answer.add(block);
					if(pathIdx == debugIdx) System.out.println("Path: "+pathIdx+" Adding unphased block with "+nextUnphasedBlock.size()+" reads");
					nextUnphasedBlock = new HashSet<Integer>();
				}
				i++;
			}
			if(addRead) {
				if(pathIdx == debugIdx) System.out.print("Next aln: "+aln);
				if(pathIdx == debugIdx && i< phasedPathBoundaries.size()) System.out.println("Next region: "+phasedPathBoundaries.get(i).getFirst()+"-"+phasedPathBoundaries.get(i).getLast());
				else if(pathIdx == debugIdx) System.out.println("Region end");
				nextUnphasedBlock.add(aln.getReadNumber());
			}
		}
		if(nextUnphasedBlock.size()>0) {
			PhasedPathBlock block = new PhasedPathBlock();
			block.addPhasedReadIds(nextUnphasedBlock);
			answer.add(block);
			if(pathIdx == debugIdx) System.out.println("Path: "+pathIdx+" Adding unphased block with "+nextUnphasedBlock.size()+" reads");
		}
		return answer;
	}

	private void savePathFiles(String outPrefix, String sequenceName, String consensus, List<ReadAlignment> alignments, List<CalledGenomicVariant> hetSNVs) {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList sequences = new QualifiedSequenceList();
		sequences.add(new QualifiedSequence(sequenceName,consensus));
		try (PrintStream out=new PrintStream(outPrefix+".fa")) {
			handler.saveSequences(sequences, out, 100);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		}
		try (PrintStream out=new PrintStream(outPrefix+".bam");
			 ReadAlignmentFileWriter readAlnsWriter = new ReadAlignmentFileWriter(sequences, out);) {
			for(ReadAlignment aln:alignments) readAlnsWriter.write(aln);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		}
		VCFFileWriter vcfWriter = new VCFFileWriter();
		VCFFileHeader header =VCFFileHeader.makeDefaultEmptyHeader();
		header.addDefaultSample("sample");
		try (PrintStream out=new PrintStream(outPrefix+".vcf");
			 ) {
			vcfWriter.printHeader(header, out);
			for(CalledGenomicVariant call:hetSNVs) {
				vcfWriter.printVCFRecord(new VCFRecord(call, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, call, header), out);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		}
	}
	
	private List<Set<Integer>> mergeClustersNaive(AssemblyGraph graph, int pathIdx, List<AssemblyEdge> path, List<ReadAlignment> alignments, List<List<ReadAlignment>> clusters) {
		Set<Integer> cluster0 = new HashSet<Integer>();
		Set<Integer> cluster1 = new HashSet<Integer>();
		List<Set<Integer>> answer = new ArrayList<Set<Integer>>();
		answer.add(cluster0);
		answer.add(cluster1);
		for(int i=0;i<clusters.size();i++) {
			List<ReadAlignment> clusterAlns = clusters.get(i);
			for(ReadAlignment aln:clusterAlns) {
				if(pathIdx == debugIdx) System.out.println("Path: "+aln.getSequenceName()+" input cluster: "+i+" Read "+aln.getReadNumber()+" "+aln.getReadName()+" added to cluster: "+(i%2));
				answer.get(i%2).add(aln.getReadNumber());
			}
		}
		return answer;
	}

	private List<Set<Integer>> mergeClustersWithPath(AssemblyGraph graph, int pathId, List<AssemblyEdge> path, List<ReadAlignment> allAlignments, List<List<ReadAlignment>> clusters) {
		//Index read assignments
		Map<Integer,Integer> clustersByReadId = new HashMap<Integer, Integer>();
		for(int i=0;i<clusters.size();i++) {
			List<ReadAlignment> cluster = clusters.get(i);
			for(ReadAlignment aln:cluster) {
				clustersByReadId.put(aln.getReadNumber(), i);
				
			}
		}
		
		Set<Integer> cluster0 = new HashSet<Integer>();
		Set<Integer> cluster1 = new HashSet<Integer>();
		Map<Integer,Integer> pathReadsAssignments = new HashMap<Integer, Integer>();
		Map<Integer,Integer> inputClusterAssignments = new HashMap<Integer, Integer>();
		int lastAssignment = -1;
		for(AssemblyEdge edge:path) {
			if(!edge.isSameSequenceEdge()) continue;
			int readId = edge.getVertex1().getSequenceIndex();
			int clusterId = clustersByReadId.getOrDefault(readId,-1);
			Integer assignment;
			if(inputClusterAssignments.size()==0) {
				//First edge
				if(pathId==debugIdx) System.out.println("Adding read "+graph.getSequence(readId).getName()+" to cluster 0");
				cluster0.add(readId);
				if(clusterId>=0) inputClusterAssignments.put(clusterId, 0);
				if(clusterId%2==0) inputClusterAssignments.put(clusterId+1, 1);
				else if(clusterId>=0) inputClusterAssignments.put(clusterId-1, 1);
				assignment=0;
			} else {
				assignment = inputClusterAssignments.get(clusterId);
				if(assignment==null) {
					assignment = lastAssignment;
					int opposite = 1-lastAssignment;
					if(clusterId>=0) inputClusterAssignments.put(clusterId, assignment);
					if(clusterId%2==0) inputClusterAssignments.put(clusterId+1, opposite);
					else if(clusterId>=0) inputClusterAssignments.put(clusterId-1, opposite);
				}
				if (assignment==0) cluster0.add(readId);
				else cluster1.add(readId);
				if(pathId==debugIdx) System.out.println("Adding path read "+graph.getSequence(readId).getName()+" to cluster "+assignment+" input cluster: "+clusterId);
			}
			pathReadsAssignments.put(readId, assignment);
			lastAssignment = assignment;
		}
		//Look for unassigned embedded reads
		for(AssemblyEdge edge:path) {
			if(!edge.isSameSequenceEdge()) continue;
			int hostAssignment = pathReadsAssignments.get(edge.getVertex1().getSequenceIndex());
			List<AssemblyEmbedded> embeddedList = graph.getAllEmbedded(edge.getVertex1().getSequenceIndex());
			//List<AssemblyEmbedded> embeddedList = graph.getEmbeddedByHostId(vertexPreviousEdge.getSequenceIndex());
			for(AssemblyEmbedded embedded:embeddedList) {
				int readId = embedded.getSequenceId();
				int clusterId = clustersByReadId.getOrDefault(readId,-1);
				Integer assignment = inputClusterAssignments.get(clusterId);
				if(assignment==null) {
					assignment = hostAssignment;
					int opposite = 1-lastAssignment;
					if(clusterId>=0) inputClusterAssignments.put(clusterId, assignment);
					if(clusterId%2==0) inputClusterAssignments.put(clusterId+1, opposite);
					else if(clusterId>=0) inputClusterAssignments.put(clusterId-1, opposite);
				}
				if(pathId==debugIdx) System.out.println("Adding embedded read "+graph.getSequence(readId).getName()+" to cluster "+assignment+" input cluster: "+clusterId);
				if (assignment==0) cluster0.add(readId);
				else cluster1.add(readId);
			}
		}
		
		List<Set<Integer>> answer = new ArrayList<Set<Integer>>();
		answer.add(cluster0);
		answer.add(cluster1);
		return answer;
	}

	private List<CalledGenomicVariant> findHeterozygousSNVs(StringBuilder consensus, List<ReadAlignment> alignments, String sequenceName) {
		//List<GenomicRegion> activeSegments = ConsensusBuilderBidirectionalWithPolishing.calculateActiveSegments(sequenceName, alignments);
		List<GenomicRegion> activeSegments = new ArrayList<>();
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		generator.setLog(log);
		QualifiedSequenceList metadata = new QualifiedSequenceList();
		metadata.add(new QualifiedSequence(sequenceName,consensus.length()));
		generator.setSequencesMetadata(metadata);
		generator.setMaxAlnsPerStartPos(0);
		SimpleHeterozygousSNVsDetectorPileupListener hetSNVsListener = new SimpleHeterozygousSNVsDetectorPileupListener(consensus, activeSegments);
		generator.addListener(hetSNVsListener);
		
		int count = 0;
		for(ReadAlignment aln:alignments) {
			generator.processAlignment(aln);
			count++;
			if(count%1000==0) log.info("Sequence: "+sequenceName+". identified heterozygous SNVs from "+count+" alignments"); 
		}
		generator.notifyEndOfAlignments();
		log.info("Called variants in sequence: "+sequenceName+". Total heterozygous SNVs: "+hetSNVsListener.getHeterozygousSNVs().size()+" alignments: "+count);
		return hetSNVsListener.getHeterozygousSNVs();
	}

}
class ClusterReadsTask implements Runnable {
	private HaplotypeReadsClusterCalculator parent;
	private AssemblyGraph graph;
	private List<AssemblyEdge> path;
	private int pathIdx;
	private int ploidy;
	private List<PhasedPathBlock> phasedBlocks;
	public ClusterReadsTask(HaplotypeReadsClusterCalculator parent, AssemblyGraph graph, List<AssemblyEdge> path, int pathIdx, int ploidy) {
		super();
		this.parent = parent;
		this.graph = graph;
		this.path = path;
		this.pathIdx = pathIdx;
		this.ploidy = ploidy;
	}
	@Override
	public void run() {
		phasedBlocks = parent.clusterReadsPath(graph,path, pathIdx, ploidy);
	}
	
	public List<PhasedPathBlock> getPhasedBlocks() {
		return phasedBlocks;
	}
	public int getPathIdx() {
		return pathIdx;
	}
	
}
class SimpleHeterozygousSNVsDetectorPileupListener implements PileupListener {

	private StringBuilder consensus;
	private List<GenomicRegion> indelRegions;
	private List<CalledGenomicVariant> heterozygousSNVs = new ArrayList<CalledGenomicVariant>();
	private int nextIndelPos = 0;

	public SimpleHeterozygousSNVsDetectorPileupListener(StringBuilder consensus, List<GenomicRegion> indelRegions) {
		super();
		this.consensus = consensus;
		this.indelRegions = indelRegions;
	}
	
	public List<CalledGenomicVariant> getHeterozygousSNVs() {
		return heterozygousSNVs;
	}
	@Override
	public void onPileup(PileupRecord pileup) {
		int idxDebug = -1;
		int pos = pileup.getPosition();
		//Check if pileup is located within an indel region
		while(nextIndelPos<indelRegions.size()) {
			GenomicRegion region = indelRegions.get(nextIndelPos);
			if(region.getFirst()<=pos && pos<=region.getLast()) return;
			else if (pos<region.getFirst()) break;
			nextIndelPos++;
		}
		List<ReadAlignment> alns = pileup.getAlignments();
		if(pos==idxDebug) System.out.println("SimpleHetSNVs. Alignments: "+alns.size());
		//Index alignments per nucleotide call
		int n = DNASequence.BASES_STRING.length();
		Map<Character,List<ReadAlignment>> alnsPerNucleotide = new HashMap<Character, List<ReadAlignment>>(n);
		for(int i=0;i<n;i++) {
			alnsPerNucleotide.put(DNASequence.BASES_STRING.charAt(i), new ArrayList<ReadAlignment>(alns.size()));
		}
		for(ReadAlignment aln:alns) {
			CharSequence call = aln.getAlleleCall(pos);
			if(call == null) continue;
			char c = call.charAt(0);
			List<ReadAlignment> alnsAllele = alnsPerNucleotide.get(c);
			if(alnsAllele==null) continue;
			alnsAllele.add(aln);
		}
		//Extract counts from map of allele calls
		int [] acgtCounts = new int [n];
		for(int i=0;i<n;i++) {
			char c = DNASequence.BASES_STRING.charAt(i);
			acgtCounts[i] = alnsPerNucleotide.get(c).size();
			if(pos==idxDebug) System.out.println("SimpleHetSNVs. Alignments nucleotide "+c+": "+acgtCounts[i]);
		}
		int maxIdx = NumberArrays.getIndexMaximum(acgtCounts);
		int secondMaxIdx = NumberArrays.getIndexMaximum(acgtCounts, maxIdx);
		int maxCount = acgtCounts[maxIdx];
		int secondCount = acgtCounts[secondMaxIdx];
		char maxBp = DNASequence.BASES_STRING.charAt(maxIdx);
		char secondBp = DNASequence.BASES_STRING.charAt(secondMaxIdx);
		char refBase = consensus.charAt(pileup.getPosition()-1);
		char altBase = (maxBp==refBase)?secondBp:maxBp;
		if(pos==idxDebug) System.out.println("SimpleHetSNVs. Max count: "+maxCount+" secondCount: "+secondCount+" total alns: "+alns.size()+" ref: "+refBase+" maxBp: "+maxBp+" secondBp: "+secondBp);
		if(maxCount+secondCount>=alns.size()-1 && secondCount>=5 && (refBase==maxBp || refBase == secondBp)) {
			heterozygousSNVs.add(new CalledSNV(new SNV(pileup.getSequenceName(), pileup.getPosition(), refBase, altBase), CalledGenomicVariant.GENOTYPE_HETERO));
		}
	}

	@Override
	public void onSequenceStart(QualifiedSequence sequence) {
	}

	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {
	}
	
}
class ReadsClusterEdge {
	private int clusterId1;
	private int clusterId2;
	private int numEdges=0;
	private long totalScore=0;
	private int maxScore = 0;
	public ReadsClusterEdge(int clusterId1, int clusterId2) {
		super();
		this.clusterId1 = clusterId1;
		this.clusterId2 = clusterId2;
	}
	public void addAssemblyEdge(AssemblyEdge edge) {
		numEdges++;
		int score = (int) (edge.getScore()/(edge.getIndelsPerKbp()+1));
		totalScore+=score;
		maxScore = Math.max(maxScore, score);
	}
	public int getClusterId1() {
		return clusterId1;
	}
	public int getClusterId2() {
		return clusterId2;
	}
	public int getNumEdges() {
		return numEdges;
	}
	public int getScore() {
		//return (int) (totalScore/numEdges);
		return maxScore;
	}
	public static String getKey(int id1, int id2) {
		return ""+id1+" "+id2;
	}
	public String toString() {
		return ""+clusterId1+" "+clusterId2+" edges: "+numEdges+" score: "+getScore();
	}
	
}
class PhasedPathBlock {
	private List<Set<Integer>> phasedReadIds = new ArrayList<Set<Integer>>();

	public void addPhasedReadIds(Set<Integer> readIds ) {
		phasedReadIds.add(readIds);
	}
	public void addHomozygousReadId(int readId) {
		for(Set<Integer> hap:phasedReadIds) hap.add(readId);
		
	}
	public List<Set<Integer>> getPhasedReadIds() {
		return phasedReadIds;
	}
}
