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
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.haplotyping.HaplotypeBlock;
import ngsep.haplotyping.SingleIndividualHaplotyper;
import ngsep.math.Distribution;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
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
	
	private int globalPloidy = 2;
	
	
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
	
	

	public int getGlobalPloidy() {
		return globalPloidy;
	}

	public void setGlobalPloidy(int globalPloidy) {
		this.globalPloidy = globalPloidy;
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
				double rd = block.calculateReadDepth();
				double proportion = block.calculateProportion();
				List<Set<Integer>> phasedReadIdsBlock = block.getPhasedReadIds();
				boolean homozygousBlock = phasedReadIdsBlock.size()==1 && rd>1.4*averageHaploidRd;
    			boolean falsePhasedBlock = phasedReadIdsBlock.size()>1 && rd<1.5*averageHaploidRd && (rd<averageHaploidRd || Math.abs(0.5-proportion)>0.25 || block.getNumVariants()<10 );
    			//if(globalPloidy==1) falsePhasedBlock = phasedReadIdsBlock.size()>1 && rd<1.8*averageHaploidRd;
    			GenomicRegion region = block.getBlockRegion();
    			System.out.println("Path: "+pathId+" next block from "+region.getFirst()+" to "+region.getLast()+" vars: "+block.getNumVariants()+" phasedGroups: "+phasedReadIdsBlock.size() +" reads: "+block.getNumReads()+ " basepairs "+block.getTotalBasePairs()+" rd: "+rd+" proportion: "+proportion+" homozygousBlock: "+homozygousBlock+" incorrectlyPhasedBlock: "+falsePhasedBlock);
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
		path.setSequenceName(sequenceName);
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
			
			List<CalledGenomicVariant> hetVars = findHeterozygousVariants(path, alignments);
			countHetVars = hetVars.size();
			if(pathIdx == debugIdx) savePathFiles(path,alignments,hetVars);
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
			int totalBasePairs = 0;
			String seqName = null;
			int firstRegion = 0;
			int lastRegion = 0;
			for(ReadAlignment aln:alignments) {
				sequenceIds.add(aln.getReadNumber());
				totalBasePairs+=aln.getReadLength();
				seqName = aln.getSequenceName();
				firstRegion = firstRegion>0?Math.min(firstRegion, aln.getFirst()):aln.getFirst();
				lastRegion = Math.max(lastRegion, aln.getLast());
			}
			if(seqName!=null) {
				PathReadsCluster cluster = new PathReadsCluster(totalBasePairs,0);
				cluster.addReadIds(sequenceIds);
				cluster.setBlockRegion(new GenomicRegionImpl(seqName, firstRegion, lastRegion));
				answer.add(cluster);
			}
			return answer;
		}
		log.info("Path: "+sequenceName+". het vars: "+countHetVars+" alignments: "+alignments.size()+" haplotype blocks: "+haplotypeBlocks.size());
		
		//Collect phased blocks first
		Set<Integer> readIdsInPhasedBlocks = new HashSet<Integer>();
		List<GenomicRegion> phasedPathBoundaries = new ArrayList<GenomicRegion>(haplotypeBlocks.size());
		for(HaplotypeBlock block:haplotypeBlocks) {
			List<GenomicVariant> hetVars = block.getAllVariants();
			List<List<Integer>> readIdsClusters = block.getClusteredFragmentIds();
			//Not real phasing
			if(readIdsClusters.size()<2) continue;
			//Not enough variants for a correction block
			if(globalPloidy==1 && hetVars.size()<5) continue;
			List<Integer> readIdsHap0 = readIdsClusters.get(0);
			List<Integer> readIdsHap1 = readIdsClusters.get(1);
			//Not real phasing
			if(readIdsHap0.size()==0 || readIdsHap1.size()==0) continue;
			int totalBasePairs = 0;
			Set<Integer> sequenceIdsHap0 = new HashSet<Integer>();
			for(int readId:readIdsHap0) {
				ReadAlignment aln = alnsByReadId.get(readId);
				sequenceIdsHap0.add(aln.getReadNumber());
				totalBasePairs+=aln.getReadLength();
				if(aln.getReadNumber()==975) System.out.println("Adding read "+aln.getReadName()+" to first haplotype of block with "+block.getNumVariants()+" variants");
			}
			Set<Integer> sequenceIdsHap1 = new HashSet<Integer>();
			for(int readId:readIdsHap1) {
				ReadAlignment aln = alnsByReadId.get(readId);
				sequenceIdsHap1.add(aln.getReadNumber());
				totalBasePairs+=aln.getReadLength();
				//if(aln.getReadNumber()==61) System.out.println("Adding read "+aln.getReadName()+" to second haplotype of block with "+block.getNumVariants()+" variants");
			}
			GenomicRegionImpl region = new GenomicRegionImpl(hetVars.get(0).getSequenceName(), hetVars.get(0).getFirst(), hetVars.get(hetVars.size()-1).getLast());
			PathReadsCluster cluster = new PathReadsCluster(totalBasePairs,block.getCallsLenght());
			cluster.addReadIds(sequenceIdsHap0);
			cluster.addReadIds(sequenceIdsHap1);
			answer.add(cluster);
			if(pathIdx == debugIdx) System.out.println("Path: "+pathIdx+" Adding phased block from "+region.getFirst()+" to "+region.getLast()+" with "+sequenceIdsHap0.size()+" and "+sequenceIdsHap1.size()+" reads");
			readIdsInPhasedBlocks.addAll(sequenceIdsHap0);
			readIdsInPhasedBlocks.addAll(sequenceIdsHap1);
			phasedPathBoundaries.add(region);
			cluster.setBlockRegion(region);
		}
		Collections.sort(phasedPathBoundaries, GenomicRegionPositionComparator.getInstance());
		//Build unphased (haploid) blocks with reads not embedded in phased blocks
		int i=0;
		int totalBasePairs = 0;
		GenomicRegionImpl unphasedRegion = null;
		Set<Integer> nextUnphasedBlock = new HashSet<Integer>();
		for(ReadAlignment aln:alignments) {
			if(pathIdx == debugIdx) System.out.println("Next aligned read "+aln.getReadNumber() +" "+aln.getReadName()+" in phased block: "+readIdsInPhasedBlocks.contains(aln.getReadNumber()));
			if(readIdsInPhasedBlocks.contains(aln.getReadNumber())) continue;
			boolean addRead = true;
			//Find next phased block
			while(i<phasedPathBoundaries.size()) {
				GenomicRegion phasedRegion = phasedPathBoundaries.get(i);
				if(phasedRegion.getLast()>aln.getLast()) {
					if(aln.getFirst()>=phasedRegion.getFirst()) {
						addRead = false;
						if(pathIdx == debugIdx) System.out.println("Unphased aln within phased region "+phasedRegion.getFirst()+"-"+phasedRegion.getLast()+": "+aln);
						//Alignment contained in block but not assigned
						/*PhasedPathBlock block = phasedBlocksByFirst.get(region.getFirst());
						if(block!=null) {
							block.addHomozygousReadId(aln.getReadNumber());
						}*/
					}
					break;
				}
				if(nextUnphasedBlock.size()>0) {
					PathReadsCluster cluster = new PathReadsCluster(totalBasePairs,0);
					cluster.addReadIds(nextUnphasedBlock);
					cluster.setBlockRegion(unphasedRegion);
					answer.add(cluster);
					if(pathIdx == debugIdx) System.out.println("Path: "+pathIdx+" Adding unphased block from "+unphasedRegion.getFirst()+" to "+unphasedRegion.getLast()+" with "+nextUnphasedBlock.size()+" reads");
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
			PathReadsCluster cluster = new PathReadsCluster(totalBasePairs,0);
			cluster.addReadIds(nextUnphasedBlock);
			cluster.setBlockRegion(unphasedRegion);
			answer.add(cluster);
			if(pathIdx == debugIdx) System.out.println("Path: "+pathIdx+" Adding unphased block with "+nextUnphasedBlock.size()+" reads");
		}
		return answer;
	}
	
	private List<CalledGenomicVariant> findHeterozygousVariants(AssemblyPath path, List<ReadAlignment> alignments) {
		String consensus = path.getConsensus();
		String sequenceName = path.getSequenceName();
		double readDepth = 0;
		int count = 0;
		for(ReadAlignment aln:alignments) {
			aln.updateAlleleCallsInfo();
			readDepth+=aln.getReadLength();
			count++;
			if(count%1000==0) log.info("Sequence: "+sequenceName+". identified heterozygous SNVs from "+count+" alignments"); 
		}
		
		readDepth/=consensus.length();
		List<SimpleVariantsDetectorPileupListener> callers = new ArrayList<>();
		for(int i=0;i<numThreads;i++) {
			SimpleVariantsDetectorPileupListener varsListener = new SimpleVariantsDetectorPileupListener(consensus);
			varsListener.setCallIndels(true);
			callers.add(varsListener);
		}
		
		QualifiedSequenceList metadata = new QualifiedSequenceList();
		metadata.add(new QualifiedSequence(sequenceName,consensus.length()));
		List<AlignmentsPileupGenerator> generators = createGenerators(callers, metadata, numThreads,log );
		ThreadPoolExecutor pool = new ThreadPoolExecutor(numThreads, numThreads, consensus.length(), TimeUnit.SECONDS, new LinkedBlockingQueue<Runnable>());
		try {
			for(AlignmentsPileupGenerator generator:generators) {
				pool.execute(()->generator.processAlignments(alignments));
			}
			pool.shutdown();
	    	pool.awaitTermination(consensus.length(), TimeUnit.SECONDS);
	    	if(!pool.isShutdown()) {
				throw new InterruptedException("The ThreadPoolExecutor was not shutdown after an await Termination call");
			}
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
		
		List<CalledGenomicVariant> allVars = new ArrayList<>();
		List<CalledGenomicVariant> hetVars = new ArrayList<>();
		for (SimpleVariantsDetectorPileupListener caller:callers) {
			List<CalledGenomicVariant> calls = caller.getCalls();
			if(path.getPathId()==debugIdx) System.out.println("Collecting het calls for path: "+path.getPathId()+" Heterozygous Calls: "+calls.size()+" Limits: "+ (calls.size()>0?calls.get(0).getFirst():0)+" "+(calls.size()>0?calls.get(calls.size()-1).getLast():0));
			allVars.addAll(calls);
			for(CalledGenomicVariant call:calls) {
				if(call.isHeterozygous()) hetVars.add(call);
			}
			
		}
		GenomicRegionSortedCollection<GenomicRegion> denseVariantRegions = calculateDenseRegions(path,allVars,consensus.length());
		
		Collections.sort(hetVars,GenomicRegionPositionComparator.getInstance());
		List<CalledGenomicVariant> filteredVars = new ArrayList<CalledGenomicVariant>();
		int countSNVs = 0;
		double snvsRD = 0;
		List<Double> alleleDosages = new ArrayList<>();
		int n = hetVars.size();
		int distance = 20;
		for(int i=0;i<n;i++) {
			CalledGenomicVariant hetVar = hetVars.get(i);
			int varDepth = hetVar.getTotalReadDepth();
			if(path.getPathId() == debugIdx) System.out.println("Next raw heterozygous variant at "+hetVar.getSequenceName()+":"+hetVar.getFirst()+" alleles: "+hetVar.getAlleles()[0]+" "+hetVar.getAlleles()[1]);
			if(varDepth<0.5*readDepth || (varDepth<0.7*readDepth && varDepth<15) ) {
				if(path.getPathId() == debugIdx) System.out.println("Removing heterozygous variant at "+hetVar.getSequenceName()+":"+hetVar.getFirst()+" with small depth: "+varDepth+" average: "+readDepth);
				continue;
			}
			GenomicRegionSortedCollection< GenomicRegion> spanning = denseVariantRegions.findSpanningRegions(hetVar);
			if (spanning.size()>0) {
				if(path.getPathId() == debugIdx) System.out.println("Removing heterozygous variant at "+hetVar.getSequenceName()+":"+hetVar.getFirst()+" within dense region. ");
				continue;
			}
			int posBefore = (i==0)?-distance:hetVars.get(i-1).getLast();
			int posAfter = (i==n-1)?consensus.length()+distance:hetVars.get(i+1).getFirst();
			if(hetVar.getFirst()-distance>posBefore && hetVar.getLast()+distance<posAfter) {
				if(hetVar instanceof CalledSNV) {
					countSNVs++;
					filteredVars.add(hetVar);
					snvsRD+=hetVar.getTotalReadDepth();
					CalledSNV csnv = (CalledSNV)hetVar;
					double ad = (double)Math.min(csnv.getCountReference(), csnv.getCountAlternative())/(csnv.getTotalReadDepth());
					alleleDosages.add(ad);
				}
			} else if(path.getPathId() == debugIdx) System.out.println("Removing heterozygous variant at "+hetVar.getSequenceName()+":"+hetVar.getFirst()+" close to other variant. Limits next: "+posBefore+" "+posAfter);
		}
		if(countSNVs==0) return filteredVars;
		snvsRD/=countSNVs;
		Collections.sort(alleleDosages);
		double medianAD = alleleDosages.get(countSNVs/2);
		if(path.getPathId() == debugIdx) System.out.println("Median allele dosage: "+medianAD+" average RD: "+readDepth+" SNVs average rd: "+snvsRD);
		//Temporary filter while multigroup haplotyping is implemented
		List<CalledGenomicVariant> filteredVars2 = new ArrayList<CalledGenomicVariant>();
		for(CalledGenomicVariant hetVar:filteredVars) {
			CalledSNV csnv = (CalledSNV)hetVar;
			if(hetVar.getTotalReadDepth()>1.6*snvsRD) {
				if(path.getPathId() == debugIdx) System.out.println("Removing heterozygous variant at "+hetVar.getSequenceName()+":"+hetVar.getFirst()+" with abnormal depth. Average by SNVs: "+snvsRD+" variant values "+csnv.getTotalReadDepth());
				continue;
			}
			double ad = (double)Math.min(csnv.getCountReference(), csnv.getCountAlternative())/(csnv.getTotalReadDepth());
			if(hetVar.getTotalReadDepth()>1.5*readDepth && Math.abs(ad-medianAD)>0.2) {
				if(path.getPathId() == debugIdx) System.out.println("Removing heterozygous variant at "+hetVar.getSequenceName()+":"+hetVar.getFirst()+" with abnormal depth and allele dosage. Averages: "+readDepth+" "+medianAD+" variant values "+csnv.getTotalReadDepth()+" "+ad);
				continue;
			}
			filteredVars2.add(hetVar);
		}
		
		//if(countSNVs==0) return new ArrayList<CalledGenomicVariant>();
		log.info("Called variants in sequence: "+sequenceName+". Total heterozygous variants: "+hetVars.size()+" alignments: "+count+" filtered variants: "+filteredVars2.size()+" SNVs: "+countSNVs);
		return filteredVars2;
	}

	public static List<AlignmentsPileupGenerator> createGenerators(List<? extends PileupListener> callers, QualifiedSequenceList metadata, int numThreads, Logger log ) {
		List<AlignmentsPileupGenerator> generators = new ArrayList<AlignmentsPileupGenerator>();
		QualifiedSequence contigM = metadata.get(0);
		int intervalLength = contigM.getLength() / numThreads;
		int nextStart = 1;
		int nextEnd = intervalLength+1;
		for(int i=0;i<numThreads;i++) {
			AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
			generator.setLog(log);
			generator.setSequencesMetadata(metadata);
			generator.setMaxAlnsPerStartPos(0);
			generator.addListener(callers.get(i));
			generator.setQuerySeq(contigM.getName());
			generator.setQueryFirst(nextStart);
			if(i<numThreads-1) generator.setQueryLast(nextEnd);
			generators.add(generator);
			nextStart+=intervalLength;
			nextEnd +=intervalLength;
		}
		return generators;
	}
	
	private GenomicRegionSortedCollection<GenomicRegion> calculateDenseRegions(AssemblyPath path,  List<CalledGenomicVariant> allVars,int consensusLength) {
		//100bp windows with overlap 50
		GenomicRegionSortedCollection<GenomicRegion> regions = new GenomicRegionSortedCollection<>();
		int numW = consensusLength/50+2;
		int [] windows = new int [numW];
		for(CalledGenomicVariant var:allVars) {
			int first = var.getFirst();
			int d1 = first / 100;
			int w1 = d1*2;
			int w2 = w1+1;
			if(first < d1*100+50) w2 = w1-1;
			windows[w1]++;
			if(w2>=0 && w2<numW) windows[w2]++;
		}
		Distribution distValues = new Distribution(0, 100,1);
		for(int value:windows) distValues.processDatapoint(value);
		if(path.getPathId()== debugIdx) {
			System.out.println("Distribution of variants density");
			distValues.printDistributionInt(System.out);
		}
		double average = distValues.getAverage();
		double stdev = Math.sqrt(distValues.getVariance());
		double limit = average+10*Math.max(1,stdev);
		for(int i=0;i<windows.length;i++) {
			if(windows[i]>limit) {
				int first = 50*i+1;
				int last = first+99;
				if(path.getPathId()== debugIdx) System.out.println("Next dense region: "+first+" - "+last+" variants: "+windows[i]+" limit: "+limit);
				regions.add(new GenomicRegionImpl(path.getSequenceName(), first, last));
			}
		}
		return regions;
	}

	private void savePathFiles(AssemblyPath path, List<ReadAlignment> alignments, List<CalledGenomicVariant> hetSNVs) {
		String outPrefix = "debug_"+path.getPathId();
		
		/*FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList sequences = new QualifiedSequenceList();
		String sequenceName = path.getSequenceName();
		sequences.add(new QualifiedSequence(sequenceName,path.getConsensus()));
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
		}*/
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

	private double calculateAverageHaploidRD(Map<Integer, List<PathReadsCluster>> pathBlocks) {
		double sumHaploidRD = 0;
		double sumPhasedClusters = 0;
		double totalRD = 0;
		double totalClusters = 0;
		for(List<PathReadsCluster> pathCLusters:pathBlocks.values()) {
			for(PathReadsCluster cluster:pathCLusters) {
				double rd = cluster.calculateReadDepth(); 
				totalRD+=rd;
				totalClusters++;
				double blockPloidy = cluster.getPhasedReadIds().size();
				if(globalPloidy>=2) {
					if(blockPloidy<=1 || cluster.getNumVariants()<10) continue;
					double devProportion = 0.5-Math.abs(0.5-cluster.calculateProportion());
					double weight = devProportion*cluster.getNumVariants()+0.1;
					sumHaploidRD+=rd*weight/blockPloidy;
					sumPhasedClusters+=weight;
				} else {
					if(blockPloidy>1 || cluster.getNumReads()<20) continue;
					double weight = cluster.getBlockRegion().length()/10000+0.1;
					sumHaploidRD+=rd*weight;
					sumPhasedClusters+=weight;
				}
				
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
class SimpleVariantsDetectorPileupListener implements PileupListener {
	private static final int MIN_DEPTH_ALLELE = 3;
	private String consensus;
	private List<CalledGenomicVariant> calls = new ArrayList<CalledGenomicVariant>();
	private boolean callIndels = false;

	public SimpleVariantsDetectorPileupListener(String consensus) {
		super();
		this.consensus = consensus;
	}
	
	

	public boolean isCallIndels() {
		return callIndels;
	}



	public void setCallIndels(boolean callIndels) {
		this.callIndels = callIndels;
	}



	public List<CalledGenomicVariant> getCalls() {
		return calls;
	}

	@Override
	public void onPileup(PileupRecord pileup) {
		int idxDebug = -1;
		int pos = pileup.getPosition();
		char refBase = consensus.charAt(pos-1);
		if(!DNASequence.EMPTY_DNA_SEQUENCE.isInAlphabet(refBase)) return;
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
			if(pos==idxDebug) System.out.println("SimpleHetVars. Sequence name "+pileup.getSequenceName()+". Next allele: "+alleleCall+" read: "+aln.getReadName());
			if(alleleCall==null || alleleCall.length()==0) continue;
			String alleleStr = alleleCall.toString();
			if(alleleStr.length()==1) {
				//Counts for SNVs
				char calledBase = alleleStr.charAt(0);
				int idxBase = DNASequence.BASES_STRING.indexOf(calledBase);
				if(idxBase>=0) acgtCounts[idxBase]++;
			}
			//Counts for biallelic indels
			if(alleleStr.length()>1) {
				GenomicRegion r = aln.getIndelCall(pos);
				if(r==null) continue;
				indelAllele = alleleStr;
				lastRefIndel = r.getLast();
			} else refAlleleIndel = alleleStr;
			alleleCounts.compute(alleleStr, (k,v)->(v==null)?1:v+1);
		}
		if(indelAllele!=null && alleleCounts.get(indelAllele)>2) {
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
			//Fix reference allele for indel call
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
			calls.add(new CalledGenomicVariantImpl(new GenomicVariantImpl(pileup.getSequenceName(), pileup.getPosition(), alleles ), CalledGenomicVariant.GENOTYPE_HETERO));
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
		if(maxCount<MIN_DEPTH_ALLELE) return;
		if(refBase!=maxBp && refBase != secondBp) return;
		//TODO: Define better
		boolean secondCountLow = secondCount < MIN_DEPTH_ALLELE || secondCount<0.2*maxCount;
		CalledSNV calledSNV;
		if(maxBp!=refBase) {
			if(secondCountLow) {
				calledSNV = new CalledSNV(new SNV(pileup.getSequenceName(), pileup.getPosition(), refBase, altBase), CalledGenomicVariant.GENOTYPE_HOMOALT);
			} else {
				calledSNV = new CalledSNV(new SNV(pileup.getSequenceName(), pileup.getPosition(), refBase, altBase), CalledGenomicVariant.GENOTYPE_HETERO);
			}
		} else if(secondCountLow) return;
		else {
			calledSNV = new CalledSNV(new SNV(pileup.getSequenceName(), pileup.getPosition(), refBase, altBase), CalledGenomicVariant.GENOTYPE_HETERO);
		}
		calledSNV.setAllBaseCounts(acgtCounts);
		calls.add(calledSNV);
		
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
	private int totalBasePairs = 0;
	private GenomicRegion blockRegion;
	
	public PathReadsCluster(int totalBasePairs, int numVariants) {
		super();
		this.totalBasePairs = totalBasePairs;
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
	public double calculateReadDepth() {
		double length = blockRegion.length();
		double averageReadLength = totalBasePairs/numReads;
		if(numVariants>0) {
			//Add basepairs corresponding to the boundaries of the region
			length+=(0.5+0.5/numVariants)*averageReadLength;
		} else {
			length-=Math.min(0.75*length, averageReadLength);
		}
		return (double)totalBasePairs/length;
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
	
	
	public int getTotalBasePairs() {
		return totalBasePairs;
	}
	
	public void setTotalBasePairs(int totalBasePairs) {
		this.totalBasePairs = totalBasePairs;
	}
	public GenomicRegion getBlockRegion() {
		return blockRegion;
	}
	public void setBlockRegion(GenomicRegion blockRegion) {
		this.blockRegion = blockRegion;
	}
	public double calculateProportion () {
		if(readIds.size()<=1) return 1;
		double n1 = readIds.get(0).size();
		double n2 = readIds.get(1).size();
		return Math.min(n1, n2)/(n1+n2);
	}
}
