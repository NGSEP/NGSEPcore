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
		Map<Integer,List<PhasedPathBlock>> pathBlocks = new TreeMap<Integer, List<PhasedPathBlock>>();
		for(int i = 0; i < paths.size(); i++)
		{
			AssemblyPath path = paths.get(i);
			int pathId = i+1;
			path.setPathId(pathId);
			List<PhasedPathBlock> blocks = clusterReadsPath(graph,path, ploidy);
			pathBlocks.put(pathId,blocks);
		}
    	return mergePathClusters(graph,pathBlocks,ploidy);
	}
	
	List<PhasedPathBlock> clusterReadsPath(AssemblyGraph graph, AssemblyPath path, int ploidy) {
		int pathIdx = path.getPathId();
		AssemblyPathReadsAligner aligner = new AssemblyPathReadsAligner();
		aligner.setLog(log);
		aligner.setAlignEmbedded(true);
		aligner.setNumThreads(numThreads);
		aligner.alignPathReads(graph, path);
		StringBuilder rawConsensus = aligner.getConsensus();
		List<ReadAlignment> alignments = aligner.getAlignedReads();
		Set<Integer> unalignedReadIds = aligner.getUnalignedReadIds();
		List<List<ReadAlignment>> clusters = null;
		String sequenceName = "diploidPath_"+pathIdx;
		int countHetVars = 0;
		if(alignments.size()>0) {
			
			for(ReadAlignment aln:alignments) aln.setSequenceName(sequenceName);
			Collections.sort(alignments, GenomicRegionPositionComparator.getInstance());
			List<CalledGenomicVariant> hetVars = findHeterozygousVariants(rawConsensus, alignments, sequenceName);
			countHetVars = hetVars.size();
			if(pathIdx == debugIdx) savePathFiles("debug_"+pathIdx,sequenceName, rawConsensus.toString(),alignments,hetVars);
			if(countHetVars>0) {
				SingleIndividualHaplotyper sih = new SingleIndividualHaplotyper();
				sih.setAlgorithmName(SingleIndividualHaplotyper.ALGORITHM_NAME_REFHAP);
				//sih.setAlgorithmName(SingleIndividualHaplotyper.ALGORITHM_NAME_DGS);
				try {
					clusters = sih.phaseSequenceVariants(sequenceName, hetVars, alignments);
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
		log.info("Path: "+sequenceName+". het vars: "+countHetVars+" alignments: "+alignments.size()+" clusters from haplotyping: "+clusters.size());
		
		//Collect phased blocks first
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
		//Build unphased (haploid) blocks with reads not embedded in phased blocks
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
	
	private List<CalledGenomicVariant> findHeterozygousVariants(StringBuilder consensus, List<ReadAlignment> alignments, String sequenceName) {
		AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
		generator.setLog(log);
		QualifiedSequenceList metadata = new QualifiedSequenceList();
		metadata.add(new QualifiedSequence(sequenceName,consensus.length()));
		generator.setSequencesMetadata(metadata);
		generator.setMaxAlnsPerStartPos(0);
		SimpleHeterozygousVariantsDetectorPileupListener hetVarsListener = new SimpleHeterozygousVariantsDetectorPileupListener(consensus);
		hetVarsListener.setCallIndels(true);
		generator.addListener(hetVarsListener);
		
		int count = 0;
		for(ReadAlignment aln:alignments) {
			generator.processAlignment(aln);
			count++;
			if(count%1000==0) log.info("Sequence: "+sequenceName+". identified heterozygous SNVs from "+count+" alignments"); 
		}
		generator.notifyEndOfAlignments();
		List<CalledGenomicVariant> hetVars = hetVarsListener.getHeterozygousVariants();
		List<CalledGenomicVariant> filteredVars = new ArrayList<CalledGenomicVariant>();
		int countSNVs = 0;
		int n = hetVars.size();
		for(int i=0;i<n;i++) {
			CalledGenomicVariant hetVar = hetVars.get(i);
			int posBefore = (i==0)?-20:hetVars.get(i-1).getLast();
			int posAfter = (i==n-1)?consensus.length()+20:hetVars.get(i+1).getFirst();
			if(hetVar.getFirst()-20>posBefore && hetVar.getLast()+20<posAfter) {
				filteredVars.add(hetVar);
				if(hetVar instanceof CalledSNV) countSNVs++;
			}
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
				vcfWriter.printVCFRecord(new VCFRecord(call, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, call, header), out);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		}
	}
	
	private List<Set<Integer>> mergePathClusters(AssemblyGraph graph, Map<Integer,List<PhasedPathBlock>> pathBlocks, int ploidy) {
		List<AssemblyPath> paths = graph.getPaths();
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
    		AssemblyPath path = paths.get(i);
    		int pathId = path.getPathId();
    		List<PhasedPathBlock> haplotypeBlocks = pathBlocks.get(pathId);
    		if(haplotypeBlocks==null) {
    			log.warning("Haplotype blocks not found for path: "+pathId);
    			continue;
    		}
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
    				if((idMax==16 || idMax == 17) && (idMin == 4 || idMin==3)) System.out.println("Using edge "+edge+" for clusters joining. Current cluster edge:  "+clusterEdge);
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
    	//add back to input clusters unassigned and unmapped reads
    	for(int i=0;i<paths.size();i++) {
    		recoverNotClusteredReads(graph, paths.get(i), readsClusters, inputClusters);
    	}
		
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
	
	private void recoverNotClusteredReads(AssemblyGraph graph, AssemblyPath path, Map<Integer,Integer> readsClusters, List<Set<Integer>> inputClusters) {
		List<AssemblyEdge> edges = path.getEdges();
		for(int i=0;i<edges.size();i++) {
			AssemblyEdge edge = edges.get(i);
			if(!edge.isSameSequenceEdge()) continue;
			int pathSequenceId = edge.getVertex1().getSequenceIndex();
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
	private StringBuilder consensus;
	private List<CalledGenomicVariant> heterozygousVariants = new ArrayList<CalledGenomicVariant>();
	private boolean callIndels = false;

	public SimpleHeterozygousVariantsDetectorPileupListener(StringBuilder consensus) {
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
		if(maxCount+secondCount<alns.size()-1) return;
		if(secondCount < MIN_DEPTH_ALLELE) return;
		if(secondCount < 0.1*maxCount) return;
		if(refBase!=maxBp && refBase != secondBp) return;
		heterozygousVariants.add(new CalledSNV(new SNV(pileup.getSequenceName(), pileup.getPosition(), refBase, altBase), CalledGenomicVariant.GENOTYPE_HETERO));
		
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
