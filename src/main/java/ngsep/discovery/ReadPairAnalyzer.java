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
package ngsep.discovery;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import JSci.maths.statistics.PoissonDistribution;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.GenomicRegionPositionComparator;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.graphs.CliquesFinder;
import ngsep.math.Distribution;
import ngsep.math.PhredScoreHelper;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.variants.CalledCNV;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.ReadPairCalledGenomicVariant;

public class ReadPairAnalyzer {
	public static final String DEF_READGROUP = "";
	public static final int DEF_MAX_LEN_DELETION = 1000000;
	public static final int DEF_SPLIT_READ_SEED = 8;
	
	private Logger log = Logger.getLogger(ReadPairAnalyzer.class.getName());
	private int maxLengthDeletion = DEF_MAX_LEN_DELETION;
	private boolean ignoreProperPairFlag = false;
	private GenomicRegionSortedCollection<CalledCNV> duplications = new GenomicRegionSortedCollection<CalledCNV>();
	private int minMQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	private int seedSize = DEF_SPLIT_READ_SEED;
	
	private ReferenceGenome reference;
	private List<String> seqNames;
	//private GenomicRegionComparator comparator;
	
	//Objects with infomation of abnormal alignments
	private long coveredGenome=0;
	private int maxAvgInsertLength=0;
	private double weightedAvgInsertLength = 0;
	private Map<String, Distribution> insertLengthDistributions;
	private Map<String, Integer> insertLengthModes;
	private Map<String, Double> insertLengthStdevs;
	private Map<String,List<SameChromosomeAbnormalLengthAln>> deletionAlns = new TreeMap<String, List<SameChromosomeAbnormalLengthAln>>();
	private Map<String,List<SameChromosomeAbnormalLengthAln>> insertionAlns = new TreeMap<String, List<SameChromosomeAbnormalLengthAln>>();
	private Map<String,List<SameChromosomeAbnormalLengthAln>> inversionAlns = new TreeMap<String, List<SameChromosomeAbnormalLengthAln>>();
	
	
	
	public ReadPairAnalyzer() {
		super();
	}

	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	
	private void dispose() {
		deletionAlns.clear();
		insertionAlns.clear();
		inversionAlns.clear();
	}

	public boolean isIgnoreProperPairFlag() {
		return ignoreProperPairFlag;
	}

	public void setIgnoreProperPairFlag(boolean ignoreProperPairFlag) {
		this.ignoreProperPairFlag = ignoreProperPairFlag;
	}

	/**
	 * @return the minMQ
	 */
	public int getMinMQ() {
		return minMQ;
	}

	/**
	 * @param minMQ the minMQ to set
	 */
	public void setMinMQ(int minMQ) {
		this.minMQ = minMQ;
	}

	public ReferenceGenome getReference() {
		return reference;
	}

	public void setReference(ReferenceGenome reference) {
		this.reference = reference;
		this.seqNames = reference.getSequenceNamesStringList();
	}
	
	
	
	public int getMaxLengthDeletion() {
		return maxLengthDeletion;
	}

	public void setMaxLengthDeletion(int maxLengthDeletion) {
		this.maxLengthDeletion = maxLengthDeletion;
	}

	public int getSeedSize() {
		return seedSize;
	}

	public void setSeedSize(int seedSize) {
		this.seedSize = seedSize;
	}

	public GenomicRegionSortedCollection<CalledCNV> getDuplications() {
		return duplications;
	}

	public void setDuplications(GenomicRegionSortedCollection<CalledCNV> duplications) {
		this.duplications = duplications;
	}

	public List<CalledGenomicVariant> findVariants(String filename) throws IOException {
		try {
			List<CalledGenomicVariant> calls = new ArrayList<CalledGenomicVariant>();
			log.info("Calculating insert length distributions");
			calculateInsertLengthDistributions(filename);
			log.info("Calculated insert length distributions for "+(insertLengthDistributions.size()-1)+" read groups. Distributing abnormally aligned reads");
			
			resetDuplicationCounts();
			distributeReadsNonProperPair(filename);
			
			log.info("Finding deletions");
			List<? extends CalledGenomicVariant> deletions = findDeletions();
			log.info("Found "+deletions.size()+" deletion candidates");
			calls.addAll(deletions);
			
			log.info("Finding insertions");
			List<? extends CalledGenomicVariant> insertions = findInsertions();
			log.info("Found "+insertions.size()+" insertion candidates");
			calls.addAll(insertions);
			
			log.info("Finding breakpoints for identified indel candidates and identifiying new indels only based on split reads");
			List<? extends CalledGenomicVariant> splitReadIndels = analyzeSplitReads(calls,filename);
			log.info("Identified "+splitReadIndels.size()+" indel candidates only based on split reads");
			calls.addAll(splitReadIndels);
			
			log.info("Finding inversions");
			List<? extends CalledGenomicVariant> inversions = findInversions();
			log.info("Found "+inversions.size()+" inversion candidates");
			calls.addAll(inversions);
			
			//log.info("Finding traslocations");
			//List<ImpreciseCalledGenomicVariant> traslocations = findTraslocations(filename);
			//log.info("Found "+traslocations.size()+" traslocations");
			//answer.addAll(traslocations);
			log.info("Sorting list with "+calls.size()+" events");
			GenomicRegionComparator comparator = new GenomicRegionComparator(reference.getSequencesMetadata());
			Collections.sort(calls,comparator);
			return calls;
		} finally {
			dispose();
		}
	}	

	private void resetDuplicationCounts() {
		for(CalledCNV cnv:duplications) {
			cnv.setTandemFragments(0);
			cnv.setTransDupFragments(0);
		}
	}

	private void calculateInsertLengthDistributions(String filename) throws IOException {
		insertLengthDistributions = new TreeMap<String, Distribution>();
		//Default distribution for alignments without read group
		insertLengthDistributions.put(DEF_READGROUP, new Distribution(1, 200000, 1));
		
		ReadAlignmentFileReader reader = null;
		int numPairedUniqueAlnReads = 0;
		int firstPos = 0;
		int lastPos = 0;
		coveredGenome = 0;
		try {
			reader = new ReadAlignmentFileReader(filename);
			reader.setLoadMode (ReadAlignmentFileReader.LOAD_MODE_MINIMAL);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			filterFlags += ReadAlignment.FLAG_MATE_UNMAPPED;
			filterFlags += ReadAlignment.FLAG_MULTIPLE_ALN;
			filterFlags += ReadAlignment.FLAG_MATE_DIFFERENT_SEQUENCE;
			reader.setFilterFlags(filterFlags);
			reader.setRequiredFlags(ReadAlignment.FLAG_PAIRED);
			reader.setMinMQ(minMQ);
			String currentSeqName = null;
			List<String> readGroups = reader.getReadGroups();
			for(String rg:readGroups) insertLengthDistributions.put(rg, new Distribution(1, 200000, 1));
			
			Iterator<ReadAlignment> it = reader.iterator();
			while(it.hasNext()) {
				ReadAlignment aln = it.next();
				//Updating covered genome
				boolean sequenceChange = !aln.getSequenceName().equals(currentSeqName);
				if(sequenceChange) {
					if(currentSeqName!=null) {
						coveredGenome+=(lastPos-firstPos+1);
					}
					currentSeqName = aln.getSequenceName();
					firstPos = aln.getFirst();
					lastPos = aln.getLast();
				} else if (aln.getFirst()>lastPos) {
					coveredGenome+=(lastPos-firstPos+1);
					firstPos = aln.getFirst();
					lastPos = aln.getLast();
				} else if(aln.getLast()>lastPos) {
					lastPos = aln.getLast();
				}
				if(!ignoreProperPairFlag && !aln.isProperPair()) {
					continue;
				}
				if(aln.getInferredInsertSize()<=0) {
					continue;
				}
				
				Distribution dist  = getDistribution(aln);
				dist.processDatapoint(aln.getInferredInsertSize());
				numPairedUniqueAlnReads++;
				if(numPairedUniqueAlnReads%1000000==0) log.info("Processed "+numPairedUniqueAlnReads+" uniquely aligned paired-end reads with consistent reference sequence");
				//if(numPairedUniqueAlnReads%1000000==0) log.info("Last processed name: "+aln.getReadName()+" Located at: "+aln.getSequenceName()+":"+aln.getFirst()+". Flags: "+aln.getFlags()+". Insert length: "+aln.getInferredInsertSize());
			}
		} finally {
			if(reader!=null) reader.close();
		}
		if(numPairedUniqueAlnReads==0) throw new IOException("BAM file does not have paired-end reads with unique alignments. Please skip read pair analysis for this dataset");
		insertLengthModes = new TreeMap<String, Integer>();
		insertLengthStdevs = new TreeMap<String, Double>();
		maxAvgInsertLength = 0;
		weightedAvgInsertLength = 0;
		double sumWeights = 0;
		for(String rg:insertLengthDistributions.keySet()) {
			Distribution dist = insertLengthDistributions.get(rg);
			if(dist.getCount()>0) {
				int mode = (int)Math.round(dist.getMaximumBinStart());
				//double stdev = dist.getEstimatedStandardDeviationPeak(mode);
				double stdev;
				if(!ignoreProperPairFlag) stdev = Math.sqrt(dist.getVariance());
				else stdev = estimateStdevPeak(dist,mode);
				if(stdev <50) stdev = 50;
				if(stdev > mode) {
					log.warning("Estimated standard deviation of insert length "+stdev+" is larger than average insert length "+mode);
					stdev = mode;
				}
				log.info("Found "+((int)dist.getCount())+" uniquely aligned paired end reads for read group "+rg+". Estimated insert length: "+mode+ " estimated standard deviation: "+stdev);
				insertLengthModes.put(rg,mode);
				insertLengthStdevs.put(rg,stdev);
				if(mode>maxAvgInsertLength) maxAvgInsertLength = mode;
				weightedAvgInsertLength+=mode*dist.getCount();
				sumWeights+=dist.getCount();
			}
		}
		weightedAvgInsertLength/=sumWeights;
		
	}

	private double estimateStdevPeak(Distribution dist, int mode) {
		double [] distN = dist.getDistribution();
		int start = mode/2;
		if(start <0) start = 0;
		int end = mode+start;
		double sum = 0;
		double sum2 = 0;
		double n = 0;
		for(int i=start;i<distN.length && i<=end;i++) {
			sum += distN[i]*i;
			sum2 += distN[i]*i*i;
			n+=distN[i];
			//System.out.println("Next value: "+i+" next frequency "+distN[i]+" sum squares: "+sum2+" sum: "+sum+" datapoints: "+n);
			
		}
		if(n<2) return 0;
		double var = (sum2-sum*sum/n)/(n-1);
		if(var<0) return 0;
		//System.out.println("Estimated variance: "+var+" sum squares: "+sum2+" sum: "+sum+" datapoints: "+n);
		return Math.sqrt(var);
	}

	private Distribution getDistribution(ReadAlignment aln) {
		Distribution dist=insertLengthDistributions.get(getDistributionRG(aln));
		if(dist != null) return dist;
		return insertLengthDistributions.get(DEF_READGROUP);
	}
	private int getDistributionMode(ReadAlignment aln) {
		Integer mode = insertLengthModes.get(getDistributionRG(aln));
		if(mode != null) return mode;
		return insertLengthModes.get(DEF_READGROUP);
	}
	private double getDistributionSD(ReadAlignment aln) {
		Double distSD=insertLengthStdevs.get(getDistributionRG(aln));
		if(distSD != null) return distSD;
		return insertLengthStdevs.get(DEF_READGROUP);
	}
	private String getDistributionRG(ReadAlignment aln) {
		String readGroup = aln.getReadGroup();
		if(readGroup!=null) {
			return readGroup;
		}
		//Take the default read group if not found
		return DEF_READGROUP;
	}

	private void distributeReadsNonProperPair(String filename) throws IOException {
		int firstDebug = -1;
		int lastDebug = -1;
		
		ReadAlignmentFileReader reader = null;
		int numReads = 0;
		try {
			reader = new ReadAlignmentFileReader(filename);
			//reader.setLoadMode (ReadAlignmentFileReader.LOAD_MODE_FULL);
			reader.setLoadMode (ReadAlignmentFileReader.LOAD_MODE_MINIMAL);
			int filterFlags = ReadAlignment.FLAG_READ_UNMAPPED;
			filterFlags += ReadAlignment.FLAG_MULTIPLE_ALN;
			reader.setFilterFlags(filterFlags);
			reader.setRequiredFlags(ReadAlignment.FLAG_PAIRED);
			reader.setMinMQ(minMQ);
			String currentSeqName = null;
			List<SameChromosomeAbnormalLengthAln> seqDelAlns=null;
			List<SameChromosomeAbnormalLengthAln> seqInsAlns=null;
			List<SameChromosomeAbnormalLengthAln> seqInvAlns=null;
			Iterator<ReadAlignment> it = reader.iterator();
			while(it.hasNext()) {
				ReadAlignment aln = it.next();
				boolean sequenceChange = !aln.getSequenceName().equals(currentSeqName);
				if(sequenceChange) {
					if(currentSeqName!=null) {
						log.info("Finished sequence "+currentSeqName+" deletion alns: "+seqDelAlns.size()+" insertion alns: "+seqInsAlns.size()+" inversion alns: "+seqInvAlns.size());
						deletionAlns.put(currentSeqName, seqDelAlns);
						insertionAlns.put(currentSeqName, seqInsAlns);
						inversionAlns.put(currentSeqName, seqInvAlns);
					}
					currentSeqName = aln.getSequenceName();
					seqDelAlns = new ArrayList<SameChromosomeAbnormalLengthAln>();
					seqInsAlns = new ArrayList<SameChromosomeAbnormalLengthAln>();
					seqInvAlns = new ArrayList<SameChromosomeAbnormalLengthAln>();
				}
				numReads++;
				if(numReads%1000000==0) log.info("Processed "+numReads+" paired-end reads with unique alignments");
				
				int avgInsertLength = getDistributionMode(aln);
				int status = getAlignmentStatus(aln, avgInsertLength);
				if(status == 0) {
					//Proper pair
					continue;
				}
				if(aln.getFirst()>firstDebug && aln.getLast()<lastDebug) log.info("Status aln "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" is "+status+" insertLength: "+aln.getInferredInsertSize()+" avgLength: "+avgInsertLength);
				if(intersectWithDuplication (aln,avgInsertLength)) {
					if(aln.getFirst()>firstDebug && aln.getLast()<lastDebug) log.info("Aln "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" intersect with duplications");
					continue;
				}
				int length2 = aln.getReadLength()/2;
				
				if(status == 1 && aln.getInferredInsertSize()>0) {
					//Less than normal
					int predictedLength = avgInsertLength-aln.getInferredInsertSize();
					seqInsAlns.add(new SameChromosomeAbnormalLengthAln(aln.getFirst()+length2, aln.getMateFirst()+length2, predictedLength));
				} else if (status == 2 && aln.getInferredInsertSize()>0 ) {
					//More than normal
					int predictedLength = aln.getInferredInsertSize()-avgInsertLength;
					seqDelAlns.add(new SameChromosomeAbnormalLengthAln(aln.getFirst()+length2, aln.getMateFirst()+length2, predictedLength));
					if(aln.getFirst()>firstDebug && aln.getLast()<lastDebug) log.info("Predicted length deletion aln "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" is "+predictedLength+" numDelreads sequence: "+seqDelAlns.size());
				} else if (status == 3) {
					//Inversion candidate
					int invFirst = 0;
					int invLast = 0;
					
					if(!aln.isNegativeStrand() && aln.getFirst()>aln.getMateFirst()) {
						invFirst = aln.getMateFirst() + length2;
						invLast = aln.getLast() + avgInsertLength;
					} else if (aln.isNegativeStrand() && aln.getFirst()<aln.getMateFirst()) {
						invFirst = Math.max(1,aln.getFirst() - avgInsertLength);
						invLast = aln.getMateFirst();
					}
					int invLength = invLast-invFirst+1;
					if(aln.getFirst()>firstDebug && aln.getLast()<lastDebug) log.info("Aln "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" invFirst "+invFirst+" invLast: "+invLast+" length: "+invLength);
					if(invFirst>0 && invLast > 0 && invLength < maxLengthDeletion) {
						seqInvAlns.add(new SameChromosomeAbnormalLengthAln(invFirst, invLast, invLength));
					}
				}
				
			}
			if(currentSeqName!=null) {
				log.info("Finished sequence "+currentSeqName+" deletion alns: "+seqDelAlns.size()+" insertion alns: "+seqInsAlns.size()+" inversion alns: "+seqInvAlns.size());
				deletionAlns.put(currentSeqName, seqDelAlns);
				insertionAlns.put(currentSeqName, seqInsAlns);
				inversionAlns.put(currentSeqName, seqInvAlns);
			}
		} finally {
			if(reader!=null) reader.close();
			
		}
	}

	private boolean intersectWithDuplication(ReadAlignment aln, int avgInsertLength) {
		GenomicRegionSortedCollection<CalledCNV> cnvsAln = duplications.findSpanningRegions(aln);
		if(cnvsAln.size()==0) {
			cnvsAln = duplications.findSpanningRegions(aln.getMateSequenceName(), aln.getMateFirst(), aln.getMateFirst()+aln.getReadLength());
			return cnvsAln.size()>0;
		}
		//else System.out.println("Found "+cnvsAln.size()+" duplications for read "+aln.getReadName()+" ");
		for(CalledCNV cnv:cnvsAln) {
			//System.out.println("Next duplication "+cnv.getSequenceName()+":"+cnv.getFirst()+"-"+cnv.getLast()+" called by "+cnv.getSource());
			cnv.addPairedFragment(aln, avgInsertLength);
		}
		return true;
	}

	private int getAlignmentStatus(ReadAlignment aln, int avgInsertLength) {
		if(aln.isMateUnmapped()) return 4;
		int absInsert = Math.abs(aln.getInferredInsertSize());
		boolean properPair;
		if(ignoreProperPairFlag) {
			//TODO: Do this better
			double stdev = getDistributionSD(aln);
			double minInsertLength = avgInsertLength - 3*stdev;
			double maxInsertLength = avgInsertLength + 3*stdev;
			properPair = aln.isMateSameSequence();
			//Opposite strands
			properPair = properPair && aln.isNegativeStrand()!=aln.isMateNegativeStrand();
			//Ends facing each other
			properPair = properPair && (aln.isNegativeStrand() == (aln.getFirst()>aln.getMateFirst()));
			//Insert length within limits
			properPair = properPair && (minInsertLength <= absInsert && absInsert <= maxInsertLength);
			
		} else properPair = aln.isProperPair();
		if(properPair) return 0;
		if(aln.isMateDifferentSequence()) return 5;
		if(aln.isNegativeStrand()==aln.isMateNegativeStrand()) return 3;
		if(aln.isPositiveStrand() && aln.getFirst()>aln.getMateFirst()) return 6;
		if(aln.isNegativeStrand() && aln.getFirst()<aln.getMateFirst()) return 7;
		if(absInsert<avgInsertLength) return 1;
		if(absInsert>avgInsertLength && absInsert <maxLengthDeletion) return 2;
		return 8;
	}
	

	private List<ReadPairCalledGenomicVariant> findInsertions() {
		List<ReadPairCalledGenomicVariant> insertions = new ArrayList<ReadPairCalledGenomicVariant>();
		int totalInsAlns = 0;
		for(String seqName:seqNames) {
			List<SameChromosomeAbnormalLengthAln> seqInsAlns = insertionAlns.get(seqName);
			if(seqInsAlns==null) continue;
			Collections.sort(seqInsAlns);
			totalInsAlns+=seqInsAlns.size();
			List<List<SameChromosomeAbnormalLengthAln>> nonOverlappingAlns = distributeNonOverlappingAlns(seqInsAlns);
			for(List<SameChromosomeAbnormalLengthAln> overlappingAlns:nonOverlappingAlns){
				insertions.addAll(buildCandidateEvents(seqName,overlappingAlns,false));
			}
		}
		assignGenotypeQualities(insertions,totalInsAlns);
		return insertions;
	}

	private List<? extends CalledGenomicVariant> findDeletions() {
		List<ReadPairCalledGenomicVariant> deletions = new ArrayList<ReadPairCalledGenomicVariant>();
		int totalDelAlns = 0;
		for(String seqName:seqNames) {
			List<SameChromosomeAbnormalLengthAln> seqDelAlns = deletionAlns.get(seqName);
			if(seqDelAlns==null) continue;
			
			Collections.sort(seqDelAlns);
			totalDelAlns+=seqDelAlns.size();
			List<List<SameChromosomeAbnormalLengthAln>> nonOverlappingAlns = distributeNonOverlappingAlns(seqDelAlns);
			log.info("Finding deletions for sequence: "+seqName+" from "+seqDelAlns.size()+" alignments. NonOv clusters: "+nonOverlappingAlns.size());
			for(List<SameChromosomeAbnormalLengthAln> overlappingAlns:nonOverlappingAlns){
				deletions.addAll(buildCandidateEvents(seqName,overlappingAlns,true));
			}
		}
		assignGenotypeQualities(deletions,totalDelAlns);
		return deletions;
		//return makeSuperInterfaceList(deletions);
	}

	private List<List<SameChromosomeAbnormalLengthAln>> distributeNonOverlappingAlns(List<SameChromosomeAbnormalLengthAln> seqAlns) {
		List<List<SameChromosomeAbnormalLengthAln>> answer = new ArrayList<List<SameChromosomeAbnormalLengthAln>>();
		List<SameChromosomeAbnormalLengthAln> lastOverlappingAlns=null;
		int lastEnd = -1;
		for(SameChromosomeAbnormalLengthAln aln:seqAlns) {
			if(aln.getFirst()>lastEnd) {
				if(lastOverlappingAlns!=null) answer.add(lastOverlappingAlns);
				lastOverlappingAlns = new ArrayList<SameChromosomeAbnormalLengthAln>();
			}
			lastOverlappingAlns.add(aln);
			lastEnd = aln.getLast();
		}
		if(lastOverlappingAlns!=null) answer.add(lastOverlappingAlns);
		return answer;
	}

	private List<ReadPairCalledGenomicVariant> buildCandidateEvents(String seqName, List<SameChromosomeAbnormalLengthAln> overlappingAlns, boolean deletions) {
		List<ReadPairCalledGenomicVariant> answer = new ArrayList<ReadPairCalledGenomicVariant>();
		if(overlappingAlns.size()<=1) return answer;
		List<List<SameChromosomeAbnormalLengthAln>> minClusters = findMinimumClusters (seqName,overlappingAlns,deletions);
		for(List<SameChromosomeAbnormalLengthAln> cluster:minClusters) {
			answer.add(buildIndel(seqName,cluster,deletions));
		}
		return answer;
	}

	private List<List<SameChromosomeAbnormalLengthAln>> findMinimumClusters(String seqName, List<SameChromosomeAbnormalLengthAln> overlappingAlns, boolean deletions) {
		boolean [][] consistencyMatrix = buildConsistencyMatrix(overlappingAlns,deletions);
		int first = overlappingAlns.get(0).getFirst();
		int firstDebug = -1;
		int lastDebug = -1;
		if(first > firstDebug && first < lastDebug) System.out.println("Calculating potential SVs from "+overlappingAlns.size()+" overlapping alignments starting at "+seqName+":"+first);
		if(first > firstDebug && first < lastDebug) printAlignments(overlappingAlns);
		if(first > firstDebug && first < lastDebug) printMatrix(consistencyMatrix);
		List<List<Integer>> components = CliquesFinder.findCliques(consistencyMatrix);
		if(first > firstDebug && first < lastDebug) System.out.println("Number of connected components: "+components.size());
		List<List<SameChromosomeAbnormalLengthAln>> answer = new ArrayList<List<SameChromosomeAbnormalLengthAln>>();
		for(List<Integer> idxs:components) {
			List<SameChromosomeAbnormalLengthAln> nextList = new ArrayList<SameChromosomeAbnormalLengthAln>();
			int j=0;
			for(int i=0;i<overlappingAlns.size();i++) {
				SameChromosomeAbnormalLengthAln aln = overlappingAlns.get(i);
				if(j<idxs.size() && idxs.get(j)==i) {
					nextList.add(aln);
					j++;
				}
			}
			if(first > firstDebug && first < lastDebug) System.out.println("Next connected component contains "+nextList.size()+" alignments");
			answer.add(nextList);
			
		}
		
		return answer;
	}

	public void printAlignments( List<SameChromosomeAbnormalLengthAln> alns) {
		System.out.println("Overlapping alns");
		for(SameChromosomeAbnormalLengthAln aln:alns) {
			System.out.println("First: "+aln.getFirst()+" last: "+aln.getLast()+" event length: "+aln.getEventLength());
		}
	}

	public void printMatrix(boolean[][] consistencyMatrix) {
		System.out.println("Consistency matrix: ");
		for(int i=0;i<consistencyMatrix.length;i++) {
			for(int j=0;j<consistencyMatrix.length;j++) {
				System.out.print(" "+consistencyMatrix[i][j]);
			}
			System.out.println();
		}
	}

	private boolean[][] buildConsistencyMatrix(List<SameChromosomeAbnormalLengthAln> overlappingAlns, boolean deletions) {
		int nAlns = overlappingAlns.size();
		boolean matrix [] [] = new boolean [nAlns][nAlns];
		for(int i=0;i<nAlns;i++) {
			matrix[i][i] = true;
			for(int j=i+1;j<nAlns;j++) {
				 matrix[i][j] = matrix[j][i] = areConsistent(overlappingAlns.get(i),overlappingAlns.get(j),deletions);  
			}
		}
		return matrix;
	}

	private boolean areConsistent(SameChromosomeAbnormalLengthAln aln1, SameChromosomeAbnormalLengthAln aln2,boolean deletions) {
		//Consistent alignments must overlap
		int overlap = GenomicRegionSpanComparator.getInstance().getSpanLength(aln1.getFirst(), aln1.getLast(), aln2.getFirst(), aln2.getLast());
		
		if(overlap<=0) return false;
		if(deletions) {
			int span1 = aln1.getLast()-aln1.getFirst()+1;
			int span2 = aln2.getLast()-aln2.getFirst()+1;
			int avgInsert = span1-aln1.getEventLength();
			//Span + avg insert length is approximately the same as event length + 2*avgInsertLength 
			if(aln1.getEventLength()>span2+avgInsert) return false;
			if(aln2.getEventLength()>span1+avgInsert) return false;
			//Check that the overlap can accommodate the potential event
			int minEventLength = Math.min(aln1.getEventLength(),aln2.getEventLength());
			if(overlap < minEventLength) return false;
		}
		return true;
	}

	

	private ReadPairCalledGenomicVariant buildIndel(String seqName, List<SameChromosomeAbnormalLengthAln> consistentAlns, boolean deletion) {
		int firstDebug = -1;
		int lastDebug = -1;
		int first = -1;
		int last = -1;
		double sumPredictedLength = 0;
		for(SameChromosomeAbnormalLengthAln aln:consistentAlns) {
			if(first == -1 ||  aln.getFirst()>first) first = aln.getFirst();
			if(last==-1 || aln.getLast()<last) last = aln.getLast();
			sumPredictedLength+=aln.getEventLength();
			if(first > firstDebug && first < lastDebug) log.info("Creating indel from "+consistentAlns.size()+" consistent alignments at "+seqName+". Next aln limits "+aln.getFirst()+"-"+aln.getLast()+" new event coordinates "+first+"-"+last);
		}
		int avgPredictedLength = (int)Math.round(sumPredictedLength/consistentAlns.size());
		int span = last-first+1;
		if(span<0) {
			//Some reads are still inconsistent
			System.err.println("Negative span "+span+" found for predicted indel at "+seqName+":"+last+"-"+first+". Deletion? "+deletion);
			int tmp = last;
			last = first;
			first = tmp;
			span = last-first+1;
		}
		if(deletion && span < avgPredictedLength ) {
			System.err.println("Span "+span+" smaller than predicted length "+avgPredictedLength+" for deletion at "+seqName+":"+first+"-"+last);
			int r = avgPredictedLength - span;
			first -=r;
			last += r;
		}
		byte type = GenomicVariant.TYPE_LARGEINS;
		if(deletion) type = GenomicVariant.TYPE_LARGEDEL;
		GenomicVariantImpl var = new GenomicVariantImpl(seqName, first, last, type);
		var.setLength(avgPredictedLength);
		//TODO: Calculate total depth and check if heterozygous
		byte genotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
		ReadPairCalledGenomicVariant answer = new ReadPairCalledGenomicVariant(var, genotype, avgPredictedLength);
		answer.setSupportingFragments(consistentAlns.size());
		return answer;
	}
	
	/**
	 * Assigns qualities finding deviations from the average of the proportion totalReads/genomeSize.
	 * @param events to evaluate
	 * @param totalAlns Total alignments found in the data that support events to be evaluated
	 */
	private void assignGenotypeQualities(List<ReadPairCalledGenomicVariant> events, int totalAlns) {
		double avgFrags = Math.max(0.5, weightedAvgInsertLength*((double)totalAlns)/coveredGenome);
		log.info("Assigning qualities to "+events.size()+" events found from "+totalAlns+" total alignments. Average fragments null hypothesis: "+avgFrags);
		PoissonDistribution dist = new PoissonDistribution(avgFrags);
		for(ReadPairCalledGenomicVariant e:events) {
			int eventFrags = e.getSupportingFragments();
			//This is useful only for variants called from split reads 
			eventFrags+=e.getNumSplitReads();
			double cumulative = dist.cumulative(eventFrags);
			if(cumulative<0) cumulative=0;
			if(cumulative>1) cumulative=1;
			e.setGenotypeQuality(PhredScoreHelper.calculatePhredScore(1-cumulative));
			
		}
	}

	private List<ReadPairCalledGenomicVariant> analyzeSplitReads(List<CalledGenomicVariant> events, String filename) throws IOException {
		List<ReadPairCalledGenomicVariant> newEvents = new ArrayList<ReadPairCalledGenomicVariant>();
		Map<String,List<ReadPairCalledGenomicVariant>> sortedEventsPerSeq = buildSortedEventsMap (events);
		int totalSplitReads = 0;
		int totalAlns = 0;
		ReadAlignmentFileReader reader = null;
		try {
			reader = new ReadAlignmentFileReader(filename);
			reader.setFilterFlags(ReadAlignment.FLAG_SECONDARY);
			reader.setRequiredFlags(ReadAlignment.FLAG_PAIRED);
			reader.setMinMQ(minMQ);
			int i = 0;
			String currentSeqName = null;
			Map<String, ReadAlignment> unmappedPairs = new TreeMap<String, ReadAlignment>();
			List<ReadPairCalledGenomicVariant> callsSeq = new ArrayList<ReadPairCalledGenomicVariant>();
			List<List<ReadAlignment>> alnsEventsSeq = new ArrayList<List<ReadAlignment>>();
			GenomicRegionSortedCollection<GenomicRegionWithAlignment> alnsForSplitRead = new GenomicRegionSortedCollection<GenomicRegionWithAlignment>();
			Iterator<ReadAlignment> it = reader.iterator();
			while(it.hasNext()) {
				ReadAlignment aln = it.next();
				
				if(aln.isReadUnmapped() && aln.isMateUnmapped()) continue;
				totalAlns++;
				if(totalAlns%1000000==0) log.info("Processed "+totalAlns+" primary alignments with at least one end mapped");
				ReadAlignment mate = null;
				if(aln.isReadUnmapped() || aln.isMateUnmapped()) {
					mate = unmappedPairs.remove(aln.getReadName());
					if(mate == null) {
						unmappedPairs.put(aln.getReadName(), aln);
						continue;
					}
					//Switch aln with mate if the unmapped read is found first 
					if(mate.isReadUnmapped()) {
						ReadAlignment tmp = mate;
						mate = aln;
						aln = tmp;
					}
					//Fix mate orientation info in unmapped alignment
					fixMateInfo(aln,mate);
				}
				if(!aln.isReadUnmapped() && !aln.isUnique()) {
					continue;
				}
				String seqName = aln.getSequenceName();
				if(seqName == null) seqName = mate.getSequenceName();
				
				boolean seqChange = !seqName.equals(currentSeqName);
				if(seqChange) {
					for(;i<callsSeq.size();i++) findBreakpoint(callsSeq.get(i),alnsEventsSeq.get(i));
					if(alnsForSplitRead.size()>0) newEvents.addAll(buildSplitReadIndels(alnsForSplitRead,null));
					currentSeqName = seqName;
					callsSeq = sortedEventsPerSeq.get(currentSeqName);
					//If no events are found within the given sequence
					if(callsSeq == null) callsSeq = new ArrayList<ReadPairCalledGenomicVariant>();
					alnsEventsSeq = new ArrayList<List<ReadAlignment>>();
					for(int j=0;j<callsSeq.size();j++) {
						alnsEventsSeq.add(new ArrayList<ReadAlignment>());
					}
					i = 0;
				}
				
				int avgInsertLength = getDistributionMode(aln);
				GenomicRegion alnRegion = predictGenomicRegion(aln, avgInsertLength);
				if(alnRegion==null) {
					continue;
				}
				//if("ERR607633.84206".equals(alnRecord.getReadName())) System.out.println("Predicted region "+alnRegion.getSequenceName()+":"+alnRegion.getFirst()+"-"+alnRegion.getLast()+" Avg insert length: "+avgInsertLength);
				int minPossiblePosNextAlns = alnRegion.getFirst() - 3*maxAvgInsertLength;
				while(i<callsSeq.size()) {
					ReadPairCalledGenomicVariant call = callsSeq.get(i);
					int cmp = GenomicRegionPositionComparator.getInstance().compare(alnRegion, call);
					if(cmp >0 &&  call.getLast()<minPossiblePosNextAlns) {
						//No other alignments will fall within the event
						findBreakpoint(call,alnsEventsSeq.get(i));
						alnsEventsSeq.get(i).clear();
						i++;
					} else {
						break;
					}
				}
				boolean readInEvent = false;
				for(int j=i;j<callsSeq.size();j++) {
					CalledGenomicVariant call = callsSeq.get(j);
					if(GenomicRegionSpanComparator.getInstance().span(alnRegion, call)) {
						alnsEventsSeq.get(j).add(aln);
						readInEvent = true;
					} else if (alnRegion.getFirst()<call.getFirst()) {
						break;
					}
				}
				if(!readInEvent) {
					int absPredictedInsertLength = Math.abs(aln.getInferredInsertSize());
					if(aln.isReadUnmapped() || (aln.isPartialAlignment(2*seedSize+1) && absPredictedInsertLength>0 && absPredictedInsertLength<2*avgInsertLength)) {
						alnsForSplitRead.add(new GenomicRegionWithAlignment(alnRegion, aln));
						totalSplitReads++;
						if(alnsForSplitRead.size()%1000==0) newEvents.addAll(buildSplitReadIndels(alnsForSplitRead, minPossiblePosNextAlns));
					}
				}
			}
			for(;i<callsSeq.size();i++) findBreakpoint(callsSeq.get(i),alnsEventsSeq.get(i));
			if(alnsForSplitRead.size()>0) newEvents.addAll(buildSplitReadIndels(alnsForSplitRead,null));
		} finally {
			if(reader!=null) reader.close();
		}
		assignGenotypeQualities(newEvents, totalSplitReads);
		return newEvents;
	}

	private Map<String, List<ReadPairCalledGenomicVariant>> buildSortedEventsMap(List<CalledGenomicVariant> events) {
		Map<String, List<ReadPairCalledGenomicVariant>> answer = new TreeMap<String, List<ReadPairCalledGenomicVariant>>();
		for(CalledGenomicVariant call:events) {
			List<ReadPairCalledGenomicVariant> callsSeq = answer.get(call.getSequenceName());
			if(callsSeq==null) {
				callsSeq = new ArrayList<ReadPairCalledGenomicVariant>();
				answer.put(call.getSequenceName(), callsSeq);
			}
			//This will throw ClassCastException if invalid events are provided to this method
			callsSeq.add((ReadPairCalledGenomicVariant)call);
		}
		for(List<ReadPairCalledGenomicVariant> list:answer.values()) Collections.sort(list,GenomicRegionPositionComparator.getInstance());
		return answer;
	}

	private void fixMateInfo(ReadAlignment aln, ReadAlignment mate) {
		aln.setMateSequenceName(mate.getSequenceName());
		if(aln.getMateFirst()!=mate.getFirst()) {
			//log.warning("Mate start inconsistent for unmapped read "+alnRecord.getReadName()+" Mate start: "+mate.getAlignmentStart()+" start in unmapped read: "+alnRecord.getMateAlignmentStart());
			aln.setMateFirst(mate.getFirst());
		}
		if(aln.isMateNegativeStrand()!=mate.isNegativeStrand()) {
			//log.warning("Mate orientation inconsistent for unmapped read "+alnRecord.getReadName()+" Mate orientation: "+mate.getReadNegativeStrandFlag()+" orientation in unmapped read: "+alnRecord.getMateNegativeStrandFlag());
			aln.setMateNegativeStrand(mate.isNegativeStrand());
		}
	}

	private GenomicRegion predictGenomicRegion(ReadAlignment aln, int avgInsertLength) {
		GenomicRegionImpl alnRegion = null;
		int readLength = aln.getReadLength();
		int quarterRL = readLength/4;
		if(aln.isReadUnmapped()) {
			int predictedFirst;
			if(aln.isMateNegativeStrand()) {
				predictedFirst = aln.getMateFirst() - avgInsertLength;
			} else {
				predictedFirst = aln.getMateFirst() + avgInsertLength;
			}
			//if("ERR607635.192675".equals(alnRecord.getReadName())) System.out.println("Predicted first "+predictedFirst+" avgInsertLength "+avgInsertLength+" Mate start "+alnRecord.getMateAlignmentStart());
			int stdev = (int)Math.round(getDistributionSD(aln));
			alnRegion = new GenomicRegionImpl(aln.getMateSequenceName(), predictedFirst - 2*stdev, predictedFirst+readLength+2*stdev);
		} else if (/*isPartialAlignment (aln) &&*/ isWellOrientedMate(aln)) { 
			//Added 100 on each side to increase the chances of finding mid-size indels with splitreads
			alnRegion = new GenomicRegionImpl(aln.getSequenceName(), aln.getFirst()-quarterRL-100, aln.getLast()+quarterRL+100);
		}
		return alnRegion;
	}

	private boolean isWellOrientedMate(ReadAlignment aln) {
		return !aln.isMateUnmapped() && (!aln.isMateNegativeStrand() && aln.getMateFirst()<aln.getFirst()) || (aln.isMateNegativeStrand() && aln.getMateFirst()>aln.getFirst());
	}

	

	private void findBreakpoint(ReadPairCalledGenomicVariant event, List<ReadAlignment> alns) {
		//TODO: Use parameter
		int refFirst = event.getFirst()-100;
		int refLast = event.getLast()+100;
		int debugStart = -1;
		int debugEnd = -1;
		if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Finding breakpoints for indel starting at "+event.getSequenceName()+":"+event.getFirst()+"-"+event.getLast()+" using "+alns.size()+" alignments");
		if(alns.size()==0) return;
		
		CharSequence seq = reference.getReference(event.getSequenceName(), refFirst, refLast);
		if(seq==null) return;
		String referenceSeq = seq.toString().toUpperCase();
		if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Reference");
		if(refFirst > debugStart && refFirst<debugEnd) System.out.println(referenceSeq);
		int relativeFirst = 0;
		int relativeLast = referenceSeq.length()-1;
		int numSplitReads = 0;
		
		for(ReadAlignment aln: alns) {
			if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Next read: "+aln.getReadName()+". Unmapped flag: "+aln.isReadUnmapped()+". Mate unmapped flag: "+aln.isMateUnmapped()+" location "+aln.getSequenceName()+":"+aln.getFirst()+"-"+aln.getLast()+" CIGAR: "+aln.getCigarString()+" Negative strand: "+aln.isNegativeStrand()+" Mate location "+aln.getMateSequenceName()+":"+aln.getMateFirst()+" mate negative strand: "+aln.isMateNegativeStrand());
			boolean readCandidate = aln.isReadUnmapped();
			if(!readCandidate && !aln.isMateUnmapped()) { 
				int avgInsertLength = getDistributionMode(aln);
				double stdevInsertLength = getDistributionSD(aln);
				int predictedInsertLength = aln.getInferredInsertSize();
				readCandidate = predictedInsertLength>0 && predictedInsertLength<event.length()+avgInsertLength+2*stdevInsertLength;
			}
			if(readCandidate) {
				String read = aln.getReadCharacters().toString().toUpperCase();
				if(aln.isReadUnmapped() && !aln.isMateNegativeStrand()) {
					read = DNAMaskedSequence.getReverseComplement(read).toString();
				}
				if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Next split-read candidate: "+aln.getReadName()+". seqence: "+read);		
				int s = 100;
				boolean deletion = event.getType() == GenomicVariant.TYPE_LARGEDEL; 
				if(deletion) {
					s = Math.min(referenceSeq.length()-100, event.getPredictedLength());
					int minSpan = (int)Math.round(0.7*referenceSeq.length());
					if(s<minSpan) s = minSpan;
				}
				
				//Object with the last coordinate of the left side and the first coordinate of the right side 
				SameChromosomeAbnormalLengthAln[]  localAln = align(referenceSeq,read,s);
				if(localAln[0] != null) {
					int localFirst = localAln[0].getFirst();
					int localLast = localAln[0].getLast();
					boolean validAlignment = deletion && localFirst>=0 && localLast>=0;
					validAlignment = validAlignment || (!deletion && localAln[1].getFirst()>=0 && localAln[1].getLast()==-1 && localAln[1].getFirst()<read.length()-10);
					validAlignment = validAlignment || (!deletion && localAln[1].getFirst()<0 && localAln[1].getLast()>=0 && localAln[1].getLast()>10);
					if(validAlignment) {
						if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Local alignment found. Reference limits: "+localFirst+"-"+localLast);
						numSplitReads++;
						if(relativeFirst<localFirst && localFirst<relativeLast) relativeFirst = localFirst;
						if(relativeFirst<localLast && localLast<relativeLast) relativeLast = localLast;
					}
				}
				
			}
		}
		if(relativeFirst>0) {
			event.setFirst(refFirst+relativeFirst);
			event.setLast(refFirst+relativeLast);
			event.setNumSplitReads(numSplitReads);
			
		}
		if(refFirst > debugStart && refFirst<debugEnd) System.out.println("New limits for indel "+event.getSequenceName()+":"+event.getFirst()+"-"+event.getLast());
	}
	/**
	 * Returns the coordinates of a split alignment between the reference sequence and the read
	 * @param referenceSeq Template to align
	 * @param read Sequence to perform split alignment
	 * @param span Predicted length of the span. PRE: span < referenceSeq.length
	 * @return SameChromosomeAbnormalLengthAln [] Array with two positions. The first has the coordinates of the split alignment relative to the template.
	 * The second has the coordinates of the split alignment relative to the read 
	 */
	private SameChromosomeAbnormalLengthAln [] align(String referenceSeq, String read, int span) {
		SameChromosomeAbnormalLengthAln [] answer = new SameChromosomeAbnormalLengthAln[2];
		answer[0] = null;
		answer[1] = null;
		//Try three seeds
		int seed2 = seedSize/2;
		int l = read.length();
		//Look for left seed at the beginning of the sequence
		int firstS = -1;
		int firstR = -1;
		String refFirst = referenceSeq.substring(0,referenceSeq.length()-span);
		for(int i=5;(i<=5+seedSize && i+seedSize<l && firstS <0);i+=seed2) {
			String tagLeft = read.substring(i,i+seedSize);
			firstS = refFirst.indexOf(tagLeft);
			if(firstS>=0) firstR = i;
		}
		
		//Look for right seed at the end of the sequence
		int lastS = -1;
		int lastR = -1;
		String refLast = referenceSeq.substring(span);
		for(int i=l-5-seedSize;(i>=l-5-2*seedSize && i>=0 && lastS <0);i-=seed2) {
			String tagRight = read.substring(i,i+seedSize);
			lastS = refLast.lastIndexOf(tagRight);
			if(lastS>=0) {
				lastS+=span;
				lastR = i;
			}
		}
		if(firstR==-1 && lastR==-1) return answer;
		if(lastR!=-1 ) {
			//Just in case the two tags found overlap
			while (firstR >= lastR) {
				lastR++;
				lastS++;
			}
			if(lastR>=l) {
				lastR=-1;
				lastS=-1;
			}
		}
		//Get left tag as close as possible to the breakpoint
		if(firstR>=0) {
			int diffs = 0;
			boolean lastDifferent = false;
			while(lastR==-1 || (firstS+1<lastS && firstR+1<lastR)) {
				boolean diff = referenceSeq.charAt(firstS)!=read.charAt(firstR); 
				if(diff) diffs++;
				if(diffs>1) {
					//Go back to the last position of agreement
					firstS--;
					firstR--;
					if(lastDifferent) {
						//Go back two positions if the last two are mismatches
						firstS--;
						firstR--;
					}
					break;
				}
				lastDifferent = diff;
				if(firstS+1==referenceSeq.length()) break;
				if(firstR+1==l) break;
				firstS++;
				firstR++;
			}
		}
		
		//Get right tag as close as possible to the breakpoint
		if(lastR>=0) {
			int diffs = 0;
			boolean lastDifferent = false;
			while(lastR==-1 || (firstS+1<lastS && firstR+1<lastR)) {
				boolean diff = referenceSeq.charAt(lastS)!=read.charAt(lastR); 
				if(diff) diffs++;
				if(diffs>1) {
					//Go forward to the last position of agreement
					lastS++;
					lastR++;
					if(lastDifferent) {
						//Go forward two positions if the last two are mismatches
						lastS++;
						lastR++;
					}
					break;
				}
				lastDifferent = diff;
				if(lastS==0) break;
				if(lastR==0) break;
				lastS--;
				lastR--;
			}
			
		}
		answer[0] = new SameChromosomeAbnormalLengthAln(firstS, lastS, lastS-firstS+1);
		answer[1] = new SameChromosomeAbnormalLengthAln(firstR, lastR, lastR-firstR+1);
		return answer;
	}
	
	private List<ReadPairCalledGenomicVariant> buildSplitReadIndels(GenomicRegionSortedCollection<GenomicRegionWithAlignment> alnsForSplitRead, Integer limitPos) {
		List<ReadPairCalledGenomicVariant> answer = new ArrayList<ReadPairCalledGenomicVariant>();
		List<GenomicRegionWithAlignment> alnsList = alnsForSplitRead.asList();
		if(alnsList.size()==0) return answer;
		List<ReadAlignment> overlappingAlns = new ArrayList<ReadAlignment>();
		String seqName = null;
		int firstPos = -1;
		int lastPos = -1;
		for(GenomicRegionWithAlignment aln:alnsList) {
			if(firstPos==-1 || GenomicRegionSpanComparator.getInstance().span(firstPos, lastPos, aln.getFirst(), aln.getLast())) {
				overlappingAlns.add(aln.getAlignment());
				if(firstPos == -1) {
					seqName = aln.getSequenceName();
					firstPos = aln.getFirst();
				}
				if(lastPos < aln.getLast()) lastPos = aln.getLast();
			} else break;
		}
		GenomicRegionImpl indelRegion = new GenomicRegionImpl(seqName, firstPos, lastPos);
		if(limitPos ==null || indelRegion.getLast() < limitPos) {
			ReadPairCalledGenomicVariant indel = buildSplitReadIndel (indelRegion,overlappingAlns); 
			if(indel!=null) answer.add (indel);
			alnsForSplitRead.removeFirst(overlappingAlns.size());
			answer.addAll(buildSplitReadIndels(alnsForSplitRead, limitPos));
		}
		return answer;
	}

	

	

	private ReadPairCalledGenomicVariant buildSplitReadIndel(GenomicRegionImpl indelRegion, List<ReadAlignment> alns) {
		
		int refFirst = indelRegion.getFirst();
		int refLast = indelRegion.getLast();
		int debugStart = -1;
		int debugEnd = -1;
		if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Finding indel with split reads starting at "+indelRegion.getSequenceName()+":"+indelRegion.getFirst()+"-"+indelRegion.getLast()+" using "+alns.size()+" alignments");
		if(alns.size()<=1) return null;
		
		CharSequence seq = reference.getReference(indelRegion.getSequenceName(), refFirst, refLast);
		if(seq==null) return null;
		String referenceSeq = seq.toString().toUpperCase();
		if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Reference");
		if(refFirst > debugStart && refFirst<debugEnd) System.out.println(referenceSeq);
		int relativeFirstDel = 0;
		int relativeLastDel = referenceSeq.length()-1;
		int numSplitReadsDel = 0;
		int sumLengthDel = 0;
		int relativeFirstIns = 0;
		int relativeLastIns = referenceSeq.length()-1;
		int numSplitReadsIns = 0;
		int sumLengthIns = 0;
		
		for(ReadAlignment aln: alns) {
			if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Next read: "+aln.getReadName()+". Unmapped flag: "+aln.isReadUnmapped()+". Mate unmapped flag: "+aln.isMateUnmapped()+" location "+aln.getSequenceName()+":"+aln.getFirst()+"-"+aln.getLast()+" CIGAR: "+aln.getCigarString()+" Negative strand: "+aln.isNegativeStrand()+" Mate location "+aln.getMateSequenceName()+":"+aln.getMateFirst()+" mate negative strand: "+aln.isMateNegativeStrand());
			
			String read = aln.getReadCharacters().toString().toUpperCase();
			if(aln.isReadUnmapped() && !aln.isMateNegativeStrand()) {
				read = DNAMaskedSequence.getReverseComplement(read).toString();
			}
			if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Next split-read candidate: "+aln.getReadName()+". sequence: "+read);		
			//Object with the last coordinate of the left side and the first coordinate of the right side 
			SameChromosomeAbnormalLengthAln[] localAln = align(referenceSeq,read,0);
			if(localAln[0] != null && localAln[1]!=null) {
				int localRefFirst = localAln[0].getFirst();
				int localRefLast = localAln[0].getLast();
				int localRefLength = localAln[0].getEventLength();
				int localReadFirst = localAln[1].getFirst();
				int localReadLast = localAln[1].getLast();
				int localReadLength = localAln[1].getEventLength();
				
				if(localRefFirst>=0 && localRefLast>=0 && localReadFirst>=0 && localReadLast>=0 ) {
					if(refFirst > debugStart && refFirst<debugEnd) System.out.println("Local alignment found. Limits: "+localRefFirst+"-"+localRefLast+" read limits: "+localReadFirst+" "+localReadLast);
					int minEventLength = read.length()/4;
					int maxEventLengthIns = read.length()-3*seedSize;
					int diff = localReadLength-localRefLength;
					if(diff>=minEventLength && diff<=maxEventLengthIns) {
						//Insertion
						numSplitReadsIns++;
						sumLengthIns += diff;
						if(relativeFirstIns<localRefFirst && localRefFirst<relativeLastIns) relativeFirstIns = localRefFirst;
						if(relativeFirstIns<localRefLast && localRefLast<relativeLastIns) relativeLastIns = localRefLast;
					} else if (diff <= -minEventLength) {
						//Deletion
						numSplitReadsDel++;
						sumLengthDel += Math.abs(diff);
						if(relativeFirstDel<localRefFirst && localRefFirst<relativeLastDel) relativeFirstDel = localRefFirst;
						if(relativeFirstDel<localRefLast && localRefLast<relativeLastDel) relativeLastDel = localRefLast;
					}
					
				}
			}
		}
		//TODO: Calculate total depth and check if heterozygous
		byte genotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
		GenomicVariantImpl var = null;
		ReadPairCalledGenomicVariant answer = null;
		if(numSplitReadsDel>=numSplitReadsIns && numSplitReadsDel>1) {
			int predictedLength = sumLengthDel/numSplitReadsDel;
			var = new GenomicVariantImpl(indelRegion.getSequenceName(), refFirst+relativeFirstDel, refFirst+relativeLastDel, GenomicVariant.TYPE_LARGEDEL);
			var.setLength(predictedLength);
			answer = new ReadPairCalledGenomicVariant(var, genotype, predictedLength);
			answer.setNumSplitReads(numSplitReadsDel);
		}
		else if (numSplitReadsIns>1) {
			int predictedLength = sumLengthIns/numSplitReadsIns;
			var = new GenomicVariantImpl(indelRegion.getSequenceName(), refFirst+relativeFirstIns, refFirst+relativeLastIns, GenomicVariant.TYPE_LARGEINS);
			var.setLength(predictedLength);
			answer = new ReadPairCalledGenomicVariant(var, genotype, predictedLength);
			answer.setNumSplitReads(numSplitReadsIns);
		}
		return answer;
	}

	private List<ReadPairCalledGenomicVariant> findInversions() {
		List<ReadPairCalledGenomicVariant> inversions = new ArrayList<ReadPairCalledGenomicVariant>();
		int totalInvAlns = 0;
		for(String seqName:seqNames) {
			List<SameChromosomeAbnormalLengthAln> seqInvAlns = inversionAlns.get(seqName);
			if(seqInvAlns==null) continue;
			Collections.sort(seqInvAlns);
			totalInvAlns+=seqInvAlns.size();
			List<List<SameChromosomeAbnormalLengthAln>> nonOverlappingAlns = distributeNonOverlappingAlns(seqInvAlns);
			for(List<SameChromosomeAbnormalLengthAln> overlappingAlns:nonOverlappingAlns){
				ReadPairCalledGenomicVariant inv = buildInversion(seqName,overlappingAlns); 
				if (inv!=null) inversions.add(inv);
			}
		}
		
		assignGenotypeQualities(inversions,totalInvAlns);
		return inversions;
	}

	private ReadPairCalledGenomicVariant buildInversion(String seqName, List<SameChromosomeAbnormalLengthAln> overlappingAlns) {
		int first = -1;
		int last = -1;
		double sumPredictedLength = 0;
		for(SameChromosomeAbnormalLengthAln aln:overlappingAlns) {
			if(first == -1 ||  aln.getFirst()>first) first = aln.getFirst();
			if(last == -1 || aln.getLast()<last) last = aln.getLast();
			sumPredictedLength+=aln.getEventLength();
		}
		int n = overlappingAlns.size();
		int avgPredictedLength = (int)Math.round(sumPredictedLength/n);
		int span = last - first + 1;
		
		if(first < 0 || last < 0 || span < 0.5*avgPredictedLength) return null;
		byte type = GenomicVariant.TYPE_INVERSION;
		GenomicVariantImpl var = new GenomicVariantImpl(seqName, first, last, type);
		var.setLength(avgPredictedLength);
		//TODO: Calculate total depth and check if heterozygous
		byte genotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
		ReadPairCalledGenomicVariant answer = new ReadPairCalledGenomicVariant(var, genotype, avgPredictedLength);
		answer.setSupportingFragments(n);
		return answer;
	}
	
}
class SameChromosomeAbnormalLengthAln implements Comparable<SameChromosomeAbnormalLengthAln> {

	private int first;
	private int last;
	private int eventLength;
	public SameChromosomeAbnormalLengthAln(int first, int last, int eventLength) {
		this.first = first;
		this.last = last;
		this.eventLength = eventLength;
	}
	
	public int getFirst() {
		return first;
	}

	public void setFirst(int first) {
		this.first = first;
	}

	public int getLast() {
		return last;
	}

	public void setLast(int last) {
		this.last = last;
	}

	public int getEventLength() {
		return eventLength;
	}
	public void setEventLength(int eventLength) {
		this.eventLength = eventLength;
	}

	@Override
	public int compareTo(SameChromosomeAbnormalLengthAln o) {
		if(this.first!=o.first) return this.first-o.first;
		if(this.last!=o.last) return this.last-o.last;
		return this.eventLength-o.eventLength;
	}
}
class GenomicRegionWithAlignment extends GenomicRegionImpl {

	private ReadAlignment rAln;
	public GenomicRegionWithAlignment(GenomicRegion gr, ReadAlignment rAln) {
		super(gr.getSequenceName(), gr.getFirst(), gr.getLast());
		this.rAln = rAln;
	}
	public ReadAlignment getAlignment() {
		return rAln;
	}
	
	
}
