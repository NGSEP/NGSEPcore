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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileWriter;
import ngsep.discovery.AlignmentsPileupGenerator;
import ngsep.discovery.PileupListener;
import ngsep.discovery.PileupRecord;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.haplotyping.HaplotypeBlock;
import ngsep.haplotyping.SingleIndividualHaplotyper;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariantImpl;
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
		List<AssemblyPath> paths = graph.getPaths();
		Map<Integer,List<PathReadsCluster>> pathBlocks = new TreeMap<Integer, List<PathReadsCluster>>();
		for(int i = 0; i < paths.size(); i++)
		{
			AssemblyPath path = paths.get(i);
			int pathId = i+1;
			path.setPathId(pathId);
			List<PathReadsCluster> blocks = clusterReadsPath(graph,path, ploidy);
			pathBlocks.put(pathId,blocks);
		}
    	return mergePathClusters(graph,pathBlocks,ploidy);
	}
	
	public Map<Integer,ReadPathPhasingData> calculatePathReadsPhasingData(AssemblyGraph graph, int ploidy) {
		//System.out.println("HaplotypeReadsClustering. Edges printed vertex "+graph.getEdges(graph.getVertex(181, true)).size());
		List<AssemblyPath> paths = graph.getPaths();
		Map<Integer,List<PathReadsCluster>> pathBlocks = new TreeMap<Integer, List<PathReadsCluster>>();
		Map<Integer,ReadPathPhasingData> answer = new HashMap<>(graph.getNumSequences());
		for(int i = 0; i < paths.size(); i++) {
			AssemblyPath path = paths.get(i);
			int pathId = i+1;
			path.setPathId(pathId);
			List<PathReadsCluster> blocks = clusterReadsPath(graph,path, ploidy);
			pathBlocks.put(pathId,blocks);
		}
		double averageHaploidRd = calculateAverageHaploidRD(pathBlocks);
		log.info("Global haploid read depth: "+averageHaploidRd);
		for(Map.Entry<Integer,List<PathReadsCluster>> entry:pathBlocks.entrySet()) {
			int pathId = entry.getKey();
			List<PathReadsCluster> blocksPath = entry.getValue();
			int blockNumber = 0;
			for(PathReadsCluster block:blocksPath) {
				int minReads = block.getMinReads();
				double rd = block.getReadDepth();
				double proportion = block.calculateProportion();
				List<Set<Integer>> phasedReadIdsBlock = block.getPhasedReadIds();
				boolean homozygousBlock = phasedReadIdsBlock.size()==1 && rd>1.4*averageHaploidRd;
    			boolean falsePhasedBlock = phasedReadIdsBlock.size()>1 && rd<1.4*averageHaploidRd && (rd<averageHaploidRd || (Math.abs(0.5-proportion)>0.2 && minReads<20));
    			System.out.println("Path: "+pathId+ " next block with "+phasedReadIdsBlock.size() +" clusters. Total read depth: "+rd+" proportion: "+proportion+" homozygousBlock: "+homozygousBlock+" falsePhasedBlock: "+falsePhasedBlock);
    			for(int i=0;i<phasedReadIdsBlock.size();i++) {
    				Set<Integer> readIds = phasedReadIdsBlock.get(i);
    				for(int readId:readIds) {
    					ReadPathPhasingData data = new ReadPathPhasingData(readId, pathId, blockNumber);
    					data.setInHomozygousRegion(homozygousBlock);
    					data.setReadDepth(rd);
    					if(phasedReadIdsBlock.size()>1 && !falsePhasedBlock) data.setPhaseWithinBlock(i);
    					answer.put(readId,data);
    				}
    				
    			}
    			
    			blockNumber++;
			}
		}
		
    	return answer;
	}
	
	List<PathReadsCluster> clusterReadsPath(AssemblyGraph graph, AssemblyPath path, int ploidy) {
		int pathIdx = path.getPathId();
		AssemblyPathReadsAligner aligner = new AssemblyPathReadsAligner();
		aligner.setLog(log);
		aligner.setHaploid(false);
		aligner.setBuildUnalignedReadRecords(true);
		aligner.calculateConsensus(path);
		List<ReadAlignment> alignmentRecords = aligner.alignPathReads(graph, path, numThreads);
		List<HaplotypeBlock> haplotypeBlocks = null;
		String sequenceName = "diploidPath_"+pathIdx;
		int countHetVars = 0;
		List<ReadAlignment> alignments = new ArrayList<>(alignmentRecords.size());
		Map<Integer,ReadAlignment> alnsByReadId = new HashMap<>();
		Set<Integer> unalignedReadIds = new HashSet<>();
		for(ReadAlignment aln:alignmentRecords) {
			if(aln.isReadUnmapped()) unalignedReadIds.add(aln.getReadNumber());
			else {
				aln.setSequenceName(sequenceName);
				alignments.add(aln);
				alnsByReadId.put(aln.getReadNumber(),aln);
			}
			//if(aln.getReadNumber()==61) System.out.println("Read "+aln.getReadName()+" aligned to path "+pathIdx);
		}
		
		if(alignments.size()>0) {
			
			List<CalledGenomicVariant> hetVars = findHeterozygousVariants(path, alignments, sequenceName);
			countHetVars = hetVars.size();
			if(pathIdx == debugIdx) savePathFiles("debug_"+pathIdx,sequenceName, path.getConsensus(),alignments,hetVars);
			if(countHetVars>0) {
				SingleIndividualHaplotyper sih = new SingleIndividualHaplotyper();
				sih.setAlgorithmName(SingleIndividualHaplotyper.ALGORITHM_NAME_REFHAP);
				//sih.setAlgorithmName(SingleIndividualHaplotyper.ALGORITHM_NAME_DGS);
				try {
					haplotypeBlocks = sih.phaseSequenceVariants(sequenceName, hetVars, alignments);
				} catch (IOException e) {
					throw new RuntimeException (e);
				}
			}
		}
		
		List<PathReadsCluster> answer = new ArrayList<PathReadsCluster>();
		if(haplotypeBlocks == null) {
			//System.out.println("No clusters for path: "+pathIdx+". hetSNVs: "+countHetSNVs+" alignments: "+alignments.size());
			Set<Integer> sequenceIds = new HashSet<Integer>();
			double totalBasePairs = 0;
			double length = path.getConsensus().length();
			for(ReadAlignment aln:alignmentRecords) {
				sequenceIds.add(aln.getReadNumber());
				totalBasePairs+=aln.getReadLength();
			}
			PathReadsCluster block = new PathReadsCluster(totalBasePairs/length,0);
			block.addReadIds(sequenceIds);
			answer.add(block);
			return answer;
		}
		log.info("Path: "+sequenceName+". het vars: "+countHetVars+" alignments: "+alignments.size()+" haplotype blocks: "+haplotypeBlocks.size());
		
		//Collect phased blocks first
		Set<Integer> readIdsInPhasedBlocks = new HashSet<Integer>();
		List<GenomicRegion> phasedPathBoundaries = new ArrayList<GenomicRegion>(haplotypeBlocks.size());
		Map<Integer,PathReadsCluster> phasedBlocksByFirst = new HashMap<Integer, PathReadsCluster>();
		for(HaplotypeBlock block:haplotypeBlocks) {
			List<List<Integer>> readIdsClusters = block.getClusteredFragmentIds();
			//Not real phasing
			if(readIdsClusters.size()<2) continue;
			List<Integer> readIdsHap0 = readIdsClusters.get(0);
			List<Integer> readIdsHap1 = readIdsClusters.get(1);
			//Not real phasing
			if(readIdsHap0.size()==0 || readIdsHap1.size()==0) continue;
			ReadAlignment firstAln0 = alnsByReadId.get(readIdsHap0.get(0));
			ReadAlignment firstAln1 = alnsByReadId.get(readIdsHap1.get(0));
			double totalBasePairs = 0;
			GenomicRegionImpl region = new GenomicRegionImpl(firstAln0.getSequenceName(), Math.min(firstAln0.getFirst(), firstAln1.getFirst()), Math.max(firstAln0.getLast(), firstAln1.getLast()));
			Set<Integer> sequenceIdsHap0 = new HashSet<Integer>();
			for(int readId:readIdsHap0) {
				ReadAlignment aln = alnsByReadId.get(readId);
				sequenceIdsHap0.add(aln.getReadNumber());
				region.setFirst(Math.min(region.getFirst(), aln.getFirst()));
				region.setLast(Math.max(region.getLast(), aln.getLast()));
				totalBasePairs+=aln.getReadLength();
				//if(aln.getReadNumber()==61) System.out.println("Adding read "+aln.getReadName()+" to first haplotype");
			}
			Set<Integer> sequenceIdsHap1 = new HashSet<Integer>();
			for(int readId:readIdsHap1) {
				ReadAlignment aln = alnsByReadId.get(readId);
				sequenceIdsHap1.add(aln.getReadNumber());
				region.setFirst(Math.min(region.getFirst(), aln.getFirst()));
				region.setLast(Math.max(region.getLast(), aln.getLast()));
				totalBasePairs+=aln.getReadLength();
				//if(aln.getReadNumber()==61) System.out.println("Adding read "+aln.getReadName()+" to second haplotype");
			}
			//TODO Guess better 
			double totalLength  = Math.max(10000, region.length()-10000);
			PathReadsCluster cluster = new PathReadsCluster(totalBasePairs/totalLength,block.getCallsLenght());
			cluster.addReadIds(sequenceIdsHap0);
			cluster.addReadIds(sequenceIdsHap1);
			answer.add(cluster);
			if(pathIdx == debugIdx) System.out.println("Path: "+pathIdx+" Adding phased block with "+sequenceIdsHap0.size()+" and "+sequenceIdsHap1.size()+" reads and estimated length: "+totalLength);
			readIdsInPhasedBlocks.addAll(sequenceIdsHap0);
			readIdsInPhasedBlocks.addAll(sequenceIdsHap1);
			phasedBlocksByFirst.put(region.getFirst(), cluster);
			phasedPathBoundaries.add(region);
		}
		Collections.sort(phasedPathBoundaries, GenomicRegionPositionComparator.getInstance());
		//Build unphased (haploid) blocks with reads not embedded in phased blocks
		int i=0;
		double totalBasePairs = 0;
		GenomicRegionImpl unphasedRegion = null;
		Set<Integer> nextUnphasedBlock = new HashSet<Integer>();
		for(ReadAlignment aln:alignments) {
			if(readIdsInPhasedBlocks.contains(aln.getReadNumber())) continue;
			boolean addRead = true;
			//Find next phased block
			while(i<phasedPathBoundaries.size()) {
				GenomicRegion phasedRegion = phasedPathBoundaries.get(i);
				if(phasedRegion.getLast()>aln.getLast()) {
					if(aln.getFirst()>=phasedRegion.getFirst()) {
						addRead = false;
						if(pathIdx == debugIdx) System.out.print("Unphased aln within phased region "+phasedRegion.getFirst()+"-"+phasedRegion.getLast()+": "+aln);
						//Alignment contained in block but not assigned
						/*PhasedPathBlock block = phasedBlocksByFirst.get(region.getFirst());
						if(block!=null) {
							block.addHomozygousReadId(aln.getReadNumber());
						}*/
					}
					break;
				}
				if(nextUnphasedBlock.size()>0) {
					double totalLength  = Math.max(10000, unphasedRegion.length()-10000);
					PathReadsCluster block = new PathReadsCluster(totalBasePairs/totalLength,0);
					block.addReadIds(nextUnphasedBlock);
					answer.add(block);
					if(pathIdx == debugIdx) System.out.println("Path: "+pathIdx+" Adding unphased block with "+nextUnphasedBlock.size()+" reads and estimated length: "+totalLength);
					nextUnphasedBlock = new HashSet<Integer>();
					totalBasePairs = 0;
					unphasedRegion = null;
				}
				i++;
			}
			if(addRead) {
				if(pathIdx == debugIdx) System.out.print("Next aln: "+aln);
				if(pathIdx == debugIdx && i< phasedPathBoundaries.size()) System.out.println("Next region: "+phasedPathBoundaries.get(i).getFirst()+"-"+phasedPathBoundaries.get(i).getLast());
				else if(pathIdx == debugIdx) System.out.println("Region end");
				nextUnphasedBlock.add(aln.getReadNumber());
				totalBasePairs+=aln.getReadLength();
				if(unphasedRegion == null) unphasedRegion = new GenomicRegionImpl(aln.getSequenceName(), aln.getFirst(), aln.getLast());
				else {
					unphasedRegion.setFirst(Math.min(unphasedRegion.getFirst(), aln.getFirst()));
					unphasedRegion.setLast(Math.max(unphasedRegion.getLast(), aln.getLast()));
				}
			}
		}
		if(nextUnphasedBlock.size()>0) {
			PathReadsCluster block = new PathReadsCluster(totalBasePairs/unphasedRegion.length(),0);
			block.addReadIds(nextUnphasedBlock);
			answer.add(block);
			if(pathIdx == debugIdx) System.out.println("Path: "+pathIdx+" Adding unphased block with "+nextUnphasedBlock.size()+" reads");
		}
		return answer;
	}
	
	private List<CalledGenomicVariant> findHeterozygousVariants(AssemblyPath path, List<ReadAlignment> alignments, String sequenceName) {
		String consensus = path.getConsensus();
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		generator.setLog(log);
		QualifiedSequenceList metadata = new QualifiedSequenceList();
		metadata.add(new QualifiedSequence(sequenceName,consensus.length()));
		generator.setSequencesMetadata(metadata);
		generator.setMaxAlnsPerStartPos(0);
		SimpleHeterozygousVariantsDetectorPileupListener hetVarsListener = new SimpleHeterozygousVariantsDetectorPileupListener(consensus);
		hetVarsListener.setCallIndels(true);
		generator.addListener(hetVarsListener);
		
		double readDepth = 0;
		int count = 0;
		for(ReadAlignment aln:alignments) {
			generator.processAlignment(aln);
			readDepth+=aln.getReadLength();
			count++;
			if(count%1000==0) log.info("Sequence: "+sequenceName+". identified heterozygous SNVs from "+count+" alignments"); 
		}
		generator.notifyEndOfAlignments();
		readDepth/=consensus.length();
		List<CalledGenomicVariant> hetVars = hetVarsListener.getHeterozygousVariants();
		List<CalledGenomicVariant> filteredVars = new ArrayList<CalledGenomicVariant>();
		int countSNVs = 0;
		int n = hetVars.size();
		for(int i=0;i<n;i++) {
			CalledGenomicVariant hetVar = hetVars.get(i);
			if(path.getPathId() == debugIdx) System.out.println("Next raw heterozygous variant at "+hetVar.getSequenceName()+":"+hetVar.getFirst()+" alleles: "+hetVar.getAlleles()[0]+" "+hetVar.getAlleles()[1]);
			if(hetVar.getTotalReadDepth()<0.7*readDepth) {
				if(path.getPathId() == debugIdx) System.out.println("Removing heterozygous variant at "+hetVar.getSequenceName()+":"+hetVar.getFirst()+" with depth: "+hetVar.getTotalReadDepth()+" average: "+readDepth);
				continue;
			}
			int posBefore = (i==0)?-50:hetVars.get(i-1).getLast();
			int posAfter = (i==n-1)?consensus.length()+50:hetVars.get(i+1).getFirst();
			if(hetVar.getFirst()-50>posBefore && hetVar.getLast()+50<posAfter) {
				if(hetVar instanceof CalledSNV) {
					countSNVs++;
					filteredVars.add(hetVar);
				}
			} else if(path.getPathId() == debugIdx) System.out.println("Removing heterozygous variant at "+hetVar.getSequenceName()+":"+hetVar.getFirst()+" close to other variant. Limits next: "+posBefore+" "+posAfter);
		}
		//if(countSNVs==0) return new ArrayList<CalledGenomicVariant>();
		log.info("Called variants in sequence: "+sequenceName+". Total heterozygous variants: "+hetVars.size()+" alignments: "+count+" filtered variants: "+filteredVars.size()+" SNVs: "+countSNVs);
		return filteredVars;
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
				vcfWriter.printVCFRecord(new VCFRecord(call, VCFRecord.DEF_FORMAT_ARRAY_NGSEP_SNV, call, header), out);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		}
	}
	
	private List<Set<Integer>> mergePathClusters(AssemblyGraph graph, Map<Integer,List<PathReadsCluster>> pathBlocks, int ploidy) {
		List<AssemblyPath> paths = graph.getPaths();
		log.info("Merging reads from "+paths.size()+" diploid paths");
		double averageHaploidRd = calculateAverageHaploidRD(pathBlocks);
		log.info("Global haploid read depth: "+averageHaploidRd);
		List<Set<Integer>> inputClusters = new ArrayList<Set<Integer>>();
    	//Reads of each cluster
		Map<Integer,Integer> readsClusters = new HashMap<Integer, Integer>();
		List<AssemblyVertex> verticesPaths = new ArrayList<AssemblyVertex>();
    	//Clusters that can not go in the same supercluster
    	Map<Integer,List<Integer>> clusterRestrictions = new HashMap<Integer, List<Integer>>();
    	Set<Integer> readIdsForAllHaps = new HashSet<>();
    	int inputClusterId = 0;
    	for(int i=0;i<paths.size();i++) {
    		
    		int firstIdCluster = inputClusterId;
    		AssemblyPath path = paths.get(i);
    		int pathId = path.getPathId();
    		List<PathReadsCluster> haplotypeBlocks = pathBlocks.get(pathId);
    		if(haplotypeBlocks==null) {
    			log.warning("Haplotype blocks not found for path: "+pathId);
    			continue;
    		}
    		log.info("Calculated "+haplotypeBlocks.size()+" blocks from haplotyping for path: "+pathId);
    		for(PathReadsCluster block:haplotypeBlocks) {
    			int firstIdClusterBlock = inputClusterId;
    			List<Set<Integer>> inputClustersBlock = block.getPhasedReadIds();
    			double rd = block.getReadDepth();
    			double proportion = block.calculateProportion(); 
    			boolean addToAll = inputClustersBlock.size()==1 && rd>1.4*averageHaploidRd;
    			boolean merge = inputClustersBlock.size()>1 && rd<1.4*averageHaploidRd /*&& proportion < 0.3*/;
    			System.out.println("Path: "+pathId+ " next block with "+inputClustersBlock.size() +" clusters. Total read depth: "+rd+" proportion: "+proportion+" addToAll: "+addToAll+" merge: "+merge);
    			Set<Integer> mergedCluster = new HashSet<>();
    			for(Set<Integer> inputCluster:inputClustersBlock) {
    				if(addToAll) {
    					readIdsForAllHaps.addAll(inputCluster);
    					for(int readId:inputCluster) { 
            				if(pathId==debugIdx)  System.out.println("Add to all read: "+readId+" "+graph.getSequence(readId).getName());
            			}
    				} else {
    					for(int readId:inputCluster) { 
            				readsClusters.put(readId, inputClusterId);
            				if(pathId==debugIdx) System.out.println("Cluster: "+inputClusterId+" read: "+readId+" "+graph.getSequence(readId).getName());
            			}
    					if(!merge) {
    						inputClusters.add(inputCluster);
            				clusterRestrictions.put(inputClusterId, new ArrayList<Integer>());
                			inputClusterId++;
    					} else {
    						mergedCluster.addAll(inputCluster);
    					}
    				}
    			}
    			if(merge) {
    				inputClusters.add(mergedCluster);
    				clusterRestrictions.put(inputClusterId, new ArrayList<Integer>());
    				inputClusterId++;
    				continue;
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
    		
    		List<AssemblyVertex> verticesPath = path.extractVertices();
    		System.out.println("Extracted "+verticesPath.size()+" vertices for path: "+pathId);
    		verticesPaths.addAll(verticesPath);
    	}
    	//Build connections graph
    	Map<String,ReadsClusterEdge> clusterEdgesMap = new HashMap<String, ReadsClusterEdge>();
    	for (AssemblyVertex vertex:verticesPaths) {
    		int readId = vertex.getSequenceIndex();
    		Integer clusterId = readsClusters.get(readId);
    		if(clusterId ==null) continue;
    		List<AssemblyEdge> edges = graph.getEdges(vertex);
    		for(AssemblyEdge edge:edges) {
    			if(edge.isSameSequenceEdge()) continue;
    			AssemblyVertex v2 = edge.getConnectingVertex(vertex);
    			Integer id2 = readsClusters.get(v2.getSequenceIndex());
    			//if(clusterId == 20) System.out.println("MergeHaplotypeClusters. Connecting vertex: "+v2+" cluster id: "+id2+" restrictions: "+clusterRestrictions.get(clusterId));
    			if(id2!=null && clusterId!=id2 && !clusterRestrictions.get(clusterId).contains(id2) && ! clusterRestrictions.get(id2).contains(clusterId)) {
    				int idMin = Math.min(clusterId, id2);
    				int idMax = Math.max(clusterId, id2);
    				String key = ReadsClusterEdge.getKey(idMin,idMax);
    				ReadsClusterEdge clusterEdge = clusterEdgesMap.computeIfAbsent(key, (v)->new ReadsClusterEdge(idMin, idMax));
    				clusterEdge.addAssemblyEdge(edge);
    				//if((idMax==94 || idMax == 95) && (idMin == 0 || idMin==1)) System.out.println("Using edge "+edge+" for clusters joining. Current cluster edge:  "+clusterEdge);
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
    	for(int i=0;i<inputClusters.size();i++) {
    		Integer assignment = inputClustersAssignment.get(i);
    		if(assignment==null) {
    			System.out.println("Joining cluster "+i+" without edges. assignment: 0");
    			assignCluster(inputClustersAssignment, i, 0, clusterRestrictions.get(i), ploidy);
    		}
    	}
    	log.info("Finished haplotype reads clustering. Assigned "+inputClustersAssignment.size()+" input clusters");
    	//add back to input clusters unassigned and unmapped reads
    	for(int i=0;i<paths.size();i++) {
    		recoverNotClusteredReads(graph, paths.get(i), readsClusters, inputClusters, readIdsForAllHaps);
    	}
		
    	//Build answer from clusters
    	List<Set<Integer>> answer = new ArrayList<Set<Integer>>(ploidy);
    	for(int i=0;i<ploidy;i++) {
    		Set<Integer> hapSuperCluster = new HashSet<Integer>();
    		hapSuperCluster.addAll(readIdsForAllHaps);
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

	private double calculateAverageHaploidRD(Map<Integer, List<PathReadsCluster>> pathBlocks) {
		double sumHaploidRD = 0;
		double sumPhasedClusters = 0;
		double totalRD = 0;
		double totalClusters = 0;
		for(List<PathReadsCluster> pathCLusters:pathBlocks.values()) {
			for(PathReadsCluster cluster:pathCLusters) {
				totalRD+=cluster.getReadDepth();
				totalClusters++;
				double ploidy = cluster.getPhasedReadIds().size();
				if(ploidy<=1) continue;
				double devProportion = 0.5-Math.abs(0.5-cluster.calculateProportion());
				double weight = devProportion*cluster.getNumVariants()+0.1;
				sumHaploidRD+=cluster.getReadDepth()*weight/ploidy;
				sumPhasedClusters+=weight;
			}
		}
		if(sumPhasedClusters>0) return sumHaploidRD/sumPhasedClusters;
		return totalRD/(totalClusters);
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
	
	private void recoverNotClusteredReads(AssemblyGraph graph, AssemblyPath path, Map<Integer,Integer> readsClusters, List<Set<Integer>> inputClusters, Set<Integer> readIdsForAllHaps) {
		List<AssemblyEdge> edges = path.getEdges();
		for(int i=0;i<edges.size();i++) {
			AssemblyEdge edge = edges.get(i);
			if(!edge.isSameSequenceEdge()) continue;
			int pathSequenceId = edge.getVertex1().getSequenceIndex();
			if(readIdsForAllHaps.contains(pathSequenceId)) continue;
			List<AssemblyEmbedded> embeddedList = graph.getAllEmbedded(edge.getVertex1().getSequenceIndex());
			Integer clusterId = readsClusters.get(pathSequenceId);
			Map<Integer,Integer> clusterVotes = new HashMap<Integer, Integer>();
			if(clusterId == null) {
				Integer clusterIdP = (i>1)?readsClusters.get(edges.get(i-2).getVertex1().getSequenceIndex()):null;
				if(clusterIdP!=null) clusterVotes.compute(clusterIdP, (k,v)->(v==null?1:v+1));
				Integer clusterIdN = (i<edges.size()-2)?readsClusters.get(edges.get(i+2).getVertex1().getSequenceIndex()):null;
				if(clusterIdN!=null) clusterVotes.compute(clusterIdN, (k,v)->(v==null?1:v+1));
				for(AssemblyEmbedded embedded:embeddedList) {
					int readId = embedded.getSequenceId();
					Integer embeddedClusterId = readsClusters.get(readId);
					if(embeddedClusterId!=null) clusterVotes.compute(embeddedClusterId, (k,v)->(v==null?1:v+1));
				}
				int maxVotes = -1;
				for(Map.Entry<Integer, Integer> entry:clusterVotes.entrySet()) {
					if(clusterId == null || maxVotes<entry.getValue()) {
						clusterId = entry.getKey();
						maxVotes = entry.getValue();
					}
				}
				if(clusterId == null || clusterId<0 || clusterId>=inputClusters.size()) continue;
				inputClusters.get(clusterId).add(pathSequenceId);
				readsClusters.put(pathSequenceId, clusterId);
			}
			for(AssemblyEmbedded embedded:embeddedList) {
				int readId = embedded.getSequenceId();
				if(readIdsForAllHaps.contains(readId)) continue;
				Integer embeddedClusterId = readsClusters.get(readId);
				if(embeddedClusterId==null) {
					inputClusters.get(clusterId).add(readId);
					readsClusters.put(readId, clusterId);
				}
			}
		}
	}

}
class SimpleHeterozygousVariantsDetectorPileupListener implements PileupListener {
	private static final int MIN_DEPTH_ALLELE = 3;
	private String consensus;
	private List<CalledGenomicVariant> heterozygousVariants = new ArrayList<CalledGenomicVariant>();
	private boolean callIndels = false;

	public SimpleHeterozygousVariantsDetectorPileupListener(String consensus) {
		super();
		this.consensus = consensus;
	}
	
	public boolean isCallIndels() {
		return callIndels;
	}



	public void setCallIndels(boolean callIndels) {
		this.callIndels = callIndels;
	}



	public List<CalledGenomicVariant> getHeterozygousVariants() {
		return heterozygousVariants;
	}
	@Override
	public void onPileup(PileupRecord pileup) {
		int idxDebug = -1;
		int pos = pileup.getPosition();
		char refBase = consensus.charAt(pos-1);
		if(!DNASequence.isInAlphabeth(refBase)) return;
		int n = DNASequence.BASES_STRING.length();
		int [] acgtCounts = new int [n];
		List<ReadAlignment> alns = pileup.getAlignments();
		if(pos==idxDebug) System.out.println("SimpleHetVars. Sequence name: "+pileup.getSequenceName()+" Alignments: "+alns.size());
		Map<String,Integer> alleleCounts = new HashMap<String, Integer>();
		String indelAllele = null;
		String refAlleleIndel = null;
		int lastRefIndel = -1;
		for(ReadAlignment aln:alns) {
			CharSequence alleleCall = aln.getAlleleCall(pos);
			if(pos==idxDebug) System.out.println("SimpleHetVars. Sequence name "+pileup.getSequenceName()+". Next allele: "+alleleCall);
			if(alleleCall==null || alleleCall.length()==0) continue;
			String alleleStr = alleleCall.toString();
			//Counts for SNVs
			char calledBase = alleleStr.charAt(0);
			int idxBase = DNASequence.BASES_STRING.indexOf(calledBase);
			if(idxBase>=0) acgtCounts[idxBase]++;
			//Counts for biallelic indels
			if(alleleStr.length()>1) {
				GenomicRegion r = aln.getIndelCall(pos);
				if(r==null) continue;
				indelAllele = alleleStr;
				lastRefIndel = r.getLast();
			} else refAlleleIndel = alleleStr;
			alleleCounts.compute(alleleStr, (k,v)->(v==null)?1:v+1);
		}
		if(indelAllele!=null) {
			//Possible indel
			if(!callIndels) return;
			if(refAlleleIndel==null) return;
			int countNonIndelAllele = alleleCounts.get(refAlleleIndel);
			int countIndelAllele = alleleCounts.get(indelAllele);
			if(pos==idxDebug) System.out.println("SimpleHetVars. Sequence name: "+pileup.getSequenceName()+" Ref count: "+countNonIndelAllele+" indel Count: "+countIndelAllele+" total alns: "+alns.size()+" ref: "+refBase+" nonindel allele: "+refAlleleIndel+" indel allele: "+indelAllele+" last ref indel "+lastRefIndel);
			if(alleleCounts.size()>2) return;
			if(refBase!=refAlleleIndel.charAt(0)) return;
			if(refBase!=indelAllele.charAt(0)) return;
			
			if(countNonIndelAllele+countIndelAllele<alns.size()-1) return;
			if(countIndelAllele<MIN_DEPTH_ALLELE) return;
			if(countNonIndelAllele<MIN_DEPTH_ALLELE) return;
			if(countIndelAllele < 0.1*countNonIndelAllele) return;
			List<String> alleles = new ArrayList<String>(2);
			//Fix refrence allele for indel call
			if(indelAllele.length()>2) {
				//Insertion
				if(pos == consensus.length()) return;
				if(isMononucleotide(indelAllele)) return;
				refAlleleIndel = consensus.substring(pos-1,pos+1);
			} else {
				//Deletion. Recover reference end
				if(lastRefIndel <= pos+1 || lastRefIndel>=consensus.length()) return;
				refAlleleIndel = consensus.substring(pos-1,lastRefIndel);
				if(isMononucleotide(refAlleleIndel)) return;
			}
			alleles.add(refAlleleIndel);
			alleles.add(indelAllele);
			if(pos==idxDebug) System.out.println("SimpleHetVars. Adding indel. Sequence name: "+pileup.getSequenceName()+" Ref count: "+countNonIndelAllele+" indel Count: "+countIndelAllele+" total alns: "+alns.size()+" ref: "+refBase+" nonindel allele: "+refAlleleIndel+" indel allele: "+indelAllele);
			heterozygousVariants.add(new CalledGenomicVariantImpl(new GenomicVariantImpl(pileup.getSequenceName(), pileup.getPosition(), alleles ), CalledGenomicVariant.GENOTYPE_HETERO));
			
			return;
		}
		int maxIdx = NumberArrays.getIndexMaximum(acgtCounts);
		int secondMaxIdx = NumberArrays.getIndexMaximum(acgtCounts, maxIdx);
		int maxCount = acgtCounts[maxIdx];
		int secondCount = acgtCounts[secondMaxIdx];
		char maxBp = DNASequence.BASES_STRING.charAt(maxIdx);
		char secondBp = DNASequence.BASES_STRING.charAt(secondMaxIdx);
		
		char altBase = (maxBp==refBase)?secondBp:maxBp;
		if(pos==idxDebug) System.out.println("SimpleHetVars. Max count: "+maxCount+" secondCount: "+secondCount+" total alns: "+alns.size()+" ref: "+refBase+" maxBp: "+maxBp+" secondBp: "+secondBp);
		if(maxCount+secondCount<0.9*alns.size()) return;
		if(secondCount < MIN_DEPTH_ALLELE) return;
		//TODO: Define better
		if(secondCount < 0.2*maxCount) return;
		if(refBase!=maxBp && refBase != secondBp) return;
		CalledSNV calledSNV = new CalledSNV(new SNV(pileup.getSequenceName(), pileup.getPosition(), refBase, altBase), CalledGenomicVariant.GENOTYPE_HETERO);
		calledSNV.setAllBaseCounts(acgtCounts);
		
		heterozygousVariants.add(calledSNV);
		
	}

	private boolean isMononucleotide(String allele) {
		int n = allele.length();
		char firstC = allele.charAt(0);
		int firstCCount = 0;
		char lastC = allele.charAt(n-1);
		int lastCCount = 0;
		for(int i=0;i<n;i++) {
			char c = allele.charAt(i);
			if(c==firstC) firstCCount++;
			if(c==lastC) lastCCount++;
		}
		if(firstC==lastC) return firstCCount==n;
		return firstCCount==n-1 || lastCCount==n-1;
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
	private List<Integer> scores = new ArrayList<Integer>();
	public ReadsClusterEdge(int clusterId1, int clusterId2) {
		super();
		this.clusterId1 = clusterId1;
		this.clusterId2 = clusterId2;
	}
	public void addAssemblyEdge(AssemblyEdge edge) {
		numEdges++;
		int score = (int) (edge.getScore()/(edge.getIndelsPerKbp()+1));
		scores.add(score);
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
		Collections.sort(scores);
		int numMax = Math.min(3, 1+scores.size()/10);
		double score = 0;
		int n = 0;
		for(int j=scores.size()-1;j>=0 && n<numMax;j--) {
			score+=scores.get(j);
			n++;
		}
		score/=n;
		return (int) score;
	}
	public static String getKey(int id1, int id2) {
		return ""+id1+" "+id2;
	}
	public String toString() {
		return ""+clusterId1+" "+clusterId2+" edges: "+numEdges+" score: "+getScore();
	}
	
}
class PathReadsCluster {
	private List<Set<Integer>> readIds = new ArrayList<Set<Integer>>();
	private int numVariants = 0;
	private int numReads = 0;
	private int minReads = 0;
	private double readDepth = 0;
	
	public PathReadsCluster(double readDepth, int numVariants) {
		super();
		this.readDepth = readDepth;
		this.numVariants = numVariants;
	}
	public void addReadIds(Set<Integer> readIds ) {
		this.readIds.add(readIds);
		numReads+=readIds.size();
		minReads = minReads>0?Math.min(minReads, readIds.size()):readIds.size();
	}
	public List<Set<Integer>> getPhasedReadIds() {
		return readIds;
	}
	public double getReadDepth() {
		return readDepth;
	}
	
	public int getNumVariants() {
		return numVariants;
	}
	public boolean isPhased() {
		return readIds.size()>1;
	}
	
	public int getNumReads() {
		return numReads;
	}
	public int getMinReads() {
		return minReads;
	}
	public double calculateProportion () {
		if(readIds.size()<=1) return 1;
		double n1 = readIds.get(0).size();
		double n2 = readIds.get(1).size();
		return Math.min(n1, n2)/(n1+n2);
	}
}
