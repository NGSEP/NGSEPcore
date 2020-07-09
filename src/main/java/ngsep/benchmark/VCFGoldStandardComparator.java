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
package ngsep.benchmark;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.math.Distribution;
import ngsep.math.LogMath;
import ngsep.math.PhredScoreHelper;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.VariantCallReport;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

/**
 * Compares a gold standard VCF or gVCF with a VCF file with predicted calls for the same sample  
 * @author Jorge Duitama
 *
 */
public class VCFGoldStandardComparator {
	
	// Constants for default values
	public static final int DEF_MIN_CLUSTER_DISTANCE = 5;
	public static final int DEF_MAX_CLUSTER_DISTANCE = 30;
	
	// Logging and progress
	private Logger log = Logger.getLogger(VCFGoldStandardComparator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private String gsFile = null;
	private ReferenceGenome genome;
	private String outputFile = null;
	private Map<String, List<GenomicRegion>> complexRegions;
	private Map<String, List<GenomicRegion>> confidenceRegions;
	private boolean genomicVCF = false;
	
	// Model attributes
	private int mode = 0;
	private short minQuality = 0;
	private long confidenceRegionsLength = 0;
	private Map<Byte, GoldStandardComparisonCounts> countsPerType = new HashMap<>();
	private Distribution distClusterSizeGS = new Distribution(0, 10, 1);
	private Distribution distClusterTestHet = new Distribution(0, 10, 1);
	private Distribution distClusterSpan = new Distribution(0, 1000, 100);
	VCFFileWriter writer = new VCFFileWriter();
	
	// Get and set methods
	
	
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		if (log == null) throw new NullPointerException("Log can not be null");
		this.log = log;
	}

	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
	
	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	
	public String getGsFile() {
		return gsFile;
	}
	public void setGsFile(String gsFile) {
		this.gsFile = gsFile;
	}
	
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	public void setGenome(String genomeFile) throws IOException {
		setGenome(OptionValuesDecoder.loadGenome(genomeFile,log));
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public boolean isGenomicVCF() {
		return genomicVCF;
	}
	public void setGenomicVCF(boolean genomicVCF) {
		this.genomicVCF = genomicVCF;
	}
	public void setGenomicVCF(Boolean genomicVCF) {
		this.setGenomicVCF(genomicVCF.booleanValue());
	}
	
	public List<GenomicRegion> getComplexRegions(String sequenceName) {
		return complexRegions.get(sequenceName);
	}
	public void setComplexRegions(String complexRegionsFile) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		this.complexRegions = handler.loadRegionsAsMap(complexRegionsFile);
	}

	public List<GenomicRegion> getConfidenceRegions(String sequenceName) {
		return confidenceRegions.get(sequenceName);
	}
	public void setConfidenceRegions(String regionsFile) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		this.confidenceRegions = handler.loadRegionsAsMap(regionsFile);
		confidenceRegionsLength = 0;
		for(List<GenomicRegion> regions: confidenceRegions.values()) {
			for(GenomicRegion g:regions) {
				confidenceRegionsLength+=g.length();
			}
		}
	}
	
	public static void main(String[] args) throws Exception {
		VCFGoldStandardComparator instance = new VCFGoldStandardComparator();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		if(args.length>=i+2) {
			instance.mode = Integer.parseInt(args[i++]);
			instance.minQuality = Short.parseShort(args[i++]);
		}
		instance.run();
	}
	
	public void run () throws IOException {
		if (inputFile == null) throw new IOException("The input (test) VCF file is a required parameter");
		if (gsFile == null) throw new IOException("The gold standard VCF file is a required parameter");
		if (genome == null) throw new IOException("The file with the reference genome is a required parameter");
		compareFiles(gsFile,inputFile);
		if (outputFile==null) printStatistics (System.out);
		else {
			try(PrintStream out = new PrintStream(outputFile)) {
				printStatistics (out);
			}
		}
	}

	public void compareFiles(String vcfGS, String vcfTest) throws IOException {
		log.info("Comparing gold standard variant file "+vcfGS+ " with test file "+vcfTest);
		QualifiedSequenceList sequenceNames = genome.getSequencesMetadata();
		if(confidenceRegions==null && genomicVCF) {
			log.info("Loading confidence regions from gold standard VCF");
			loadConfidenceRegionsFromVCF(vcfGS,sequenceNames);
		}
		//Assume that the gold standard is perfect over the entire genome
		boolean countNonGSAsFP = false;
		if(confidenceRegions==null) {
			log.info("No confidence regions were provided. Assuming that the gold standard is complete over the entire genome");
			confidenceRegionsLength = genome.getTotalLength();
			countNonGSAsFP = true;
		}
		initCounts(countNonGSAsFP);
		try (VCFFileReader inGS = new VCFFileReader(vcfGS);
			 VCFFileReader inTest = new VCFFileReader(vcfTest)) {
			inGS.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			inGS.setSequences(sequenceNames);
			inTest.setSequences(sequenceNames);
			LinkedList<GenomicRegion> confidenceRegionsSeq=new LinkedList<>();
			List<GenomicRegion> complexRegionsSeq=null;
			int pIdx = 0;
			List<CalledGenomicVariant> gsCalls = new ArrayList<>();
			
			int sequenceIdx = -1;
			int seqLen = 0;
			int clusterFirst = 0;
			int clusterLast = 0;
			byte clusterType = GenomicVariant.TYPE_UNDETERMINED;
			
			
			Iterator<VCFRecord> itGS = inGS.iterator();
			VCFRecord recordGS = loadNextRecord(itGS, false);
			//System.out.println("First GS record: "+recordGS.getFirst());
			
			Iterator<VCFRecord> itTest = inTest.iterator();
			LinkedList<CalledGenomicVariant> testCallsSequence = new LinkedList<>();
			VCFRecord firstRecordNextSequence = loadNextRecord(itTest, false);
			
			int countProcessedClusters = 0;
			while(recordGS!=null) {
				int nextSeqIdxGS = sequenceNames.indexOf(recordGS.getSequenceName());
				if(nextSeqIdxGS>sequenceIdx) {
					if(sequenceIdx>=0) processClusterCalls (gsCalls, clusterFirst, clusterLast, clusterType, testCallsSequence, confidenceRegionsSeq, seqLen);
					countProcessedClusters++;
					sequenceIdx = nextSeqIdxGS;
					QualifiedSequence seqObj = sequenceNames.get(sequenceIdx); 
					String sequenceName = seqObj.getName();
					seqLen = seqObj.getLength();
					log.info("Starting sequence "+sequenceName);
					
					testCallsSequence.clear();
					if(firstRecordNextSequence!=null && firstRecordNextSequence.getSequenceName().equals(sequenceName)) {
						CalledGenomicVariant callTest = firstRecordNextSequence.getCalls().get(0);
						testCallsSequence.add(callTest);
						firstRecordNextSequence = loadTestCallsSequence(itTest, sequenceName, testCallsSequence);
					}
					
					if(confidenceRegions!=null) {
						confidenceRegionsSeq.clear(); 
						List<GenomicRegion> crsl = confidenceRegions.get(sequenceName);
						if(crsl!=null) confidenceRegionsSeq.addAll(crsl);
					}
					if(complexRegions!=null) complexRegionsSeq = complexRegions.get(sequenceName);
					pIdx = 0;
					gsCalls.clear();
					clusterFirst=clusterLast=0;
					clusterType = GenomicVariant.TYPE_UNDETERMINED;
				} else if (nextSeqIdxGS<sequenceIdx) {
					log.severe("Disorder detected in gold standard after sequence name: "+sequenceNames.get(sequenceIdx).getName());
					break;
				}
				if(complexRegionsSeq!=null) {
					for(;pIdx<complexRegionsSeq.size();pIdx++) {
						GenomicRegion region = complexRegionsSeq.get(pIdx);
						if(clusterLast+DEF_MIN_CLUSTER_DISTANCE<=region.getFirst()) break;
						if(clusterFirst<region.getLast()+DEF_MIN_CLUSTER_DISTANCE) {
							clusterType = GenomicVariant.TYPE_STR;
							clusterLast = Math.max(clusterLast, region.getLast());
							clusterFirst = Math.min(clusterFirst, region.getFirst());
						}
					}
				}
				int nextClusterFirst = clusterLast+Math.max(DEF_MIN_CLUSTER_DISTANCE, clusterLast-clusterFirst+1);
				nextClusterFirst = Math.min(nextClusterFirst, clusterLast+DEF_MAX_CLUSTER_DISTANCE);
				
				boolean gsClose = nextClusterFirst>recordGS.getFirst();
				if (!gsClose) {
					processClusterCalls (gsCalls, clusterFirst, clusterLast, clusterType, testCallsSequence, confidenceRegionsSeq, seqLen);
					countProcessedClusters++;
					gsCalls.clear();
					clusterFirst = clusterLast = 0;
					clusterType = GenomicVariant.TYPE_UNDETERMINED;
				}
				CalledGenomicVariant callGS = recordGS.getCalls().get(0);
				gsCalls.add(callGS);
				if(clusterType==GenomicVariant.TYPE_UNDETERMINED) clusterType = loadType(callGS);
				//Default type for cluster calls
				if(gsCalls.size()>1) clusterType= GenomicVariant.TYPE_STR;
				
				if(clusterFirst == 0 ) clusterFirst = recordGS.getFirst();
				else clusterFirst = Math.min(clusterFirst, recordGS.getFirst());
				clusterLast = Math.max(clusterLast, recordGS.getLast());
				recordGS = loadNextRecord(itGS, false);
				if(countProcessedClusters%10000==0) log.info("Processed "+countProcessedClusters+" clusters. Current cluster coordinates "+sequenceNames.get(sequenceIdx).getName()+": "+clusterFirst+"-"+clusterLast);
				if (progressNotifier!=null && countProcessedClusters%1000==0) {
					int progress = countProcessedClusters/1000;
					if (!progressNotifier.keepRunning(progress)) {
						log.info("Process canceled");
						return;
					}
				}
			}
			processClusterCalls (gsCalls, clusterFirst, sequenceNames.get(sequenceIdx).getLength(), clusterType, testCallsSequence, confidenceRegionsSeq, seqLen );
		}
	}
	
	private void loadConfidenceRegionsFromVCF(String filename, QualifiedSequenceList sequenceNames) throws IOException {
		confidenceRegions = new HashMap<>();
		confidenceRegionsLength = 0;
		try (VCFFileReader in = new VCFFileReader(filename)) {
			in.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
			in.setSequences(sequenceNames);
			String sequenceName = null;
			List<GenomicRegion> regionsSeq = null;
			GenomicRegionImpl r = null;
			Iterator<VCFRecord> it = in.iterator();
			while(it.hasNext()) {
				VCFRecord record = it.next();
				CalledGenomicVariant call = record.getCalls().get(0); 
				if(call.isUndecided()) continue;
				if(sequenceName==null || !sequenceName.equals(call.getSequenceName())) {
					if(r!=null) {
						regionsSeq.add(r);
						confidenceRegionsLength+=r.length();
					}
					sequenceName = call.getSequenceName();
					regionsSeq =  new ArrayList<>();
					confidenceRegions.put(sequenceName, regionsSeq);
					r = new GenomicRegionImpl(sequenceName, record.getFirst(), record.getLast());
				} else if( r.getLast()<record.getFirst()-1) {
					regionsSeq.add(r);
					confidenceRegionsLength+=r.length();
					r = new GenomicRegionImpl(sequenceName, record.getFirst(), record.getLast());
				} else {
					r.setLast(Math.max(r.getLast(), record.getLast()));
				}
			}
			if(r!=null) {
				regionsSeq.add(r);
				confidenceRegionsLength+=r.length();
			}
		}
		
	}

	private VCFRecord loadTestCallsSequence(Iterator<VCFRecord> itTest, String seqName, LinkedList<CalledGenomicVariant> testCallsSequence) {
		VCFRecord answer = loadNextRecord(itTest, false);
		while(answer!=null && answer.getSequenceName().equals(seqName)) {
			CalledGenomicVariant call = answer.getCalls().get(0);
			testCallsSequence.add(call);
			answer = loadNextRecord(itTest, false);
		}
		return answer;
	}

	private void initCounts(boolean countNonGSAsFP) {
		countsPerType.put(GenomicVariant.TYPE_BIALLELIC_SNV, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_INDEL, new GoldStandardComparisonCounts());
		countsPerType.put(GenomicVariant.TYPE_STR, new GoldStandardComparisonCounts());
		double confidentMbp = (double)confidenceRegionsLength/1000000.0;
		log.info("Confident MBP: "+confidentMbp);
		for(GoldStandardComparisonCounts counts:countsPerType.values()) {
			counts.setCountNonGSAsFP(countNonGSAsFP);
			counts.setConfidentMbp(confidentMbp);
		}
	}
	
	
	

	private VCFRecord loadNextRecord(Iterator<VCFRecord> it, boolean loadHomozygousReference) {
		VCFRecord record = null;
		while(it.hasNext()) {
			VCFRecord record2 = it.next();
			CalledGenomicVariant call = record2.getCalls().get(0); 
			if (call.isUndecided()) continue;
			if(!loadHomozygousReference && call.isHomozygousReference()) continue;
			record = record2;
			break;
		}
		return record;
	}

	private void processClusterCalls(List<CalledGenomicVariant> gsCalls, int clusterFirst, int clusterLast, byte clusterGSType, LinkedList<CalledGenomicVariant> testCallsSeq, LinkedList<GenomicRegion> confidenceRegionsSeq, int sequenceLength) {	
		int lastRowCounts = GoldStandardComparisonCounts.NUM_ROWS_COUNTS-1;
		if(gsCalls.size()==0) return;
		List<CalledGenomicVariant> testCallsCluster = new ArrayList<>();
		int heterozygous=0;
		int regionFirst = clusterFirst-DEF_MIN_CLUSTER_DISTANCE;
		int regionLast = clusterLast+DEF_MIN_CLUSTER_DISTANCE;
		// Process individually potential false positives
		while(testCallsSeq.size()>0) {
			CalledGenomicVariant testCall = testCallsSeq.peek();
			if(clusterFirst == 567239) log.info("Cluster limits: "+clusterFirst+"-"+clusterLast+" GS size: "+gsCalls.size()+ " cluster type: "+clusterGSType+" test call: "+testCall.getFirst()+"-"+testCall.getLast());
			if(testCall.getLast()+DEF_MIN_CLUSTER_DISTANCE<=clusterFirst) {
				processPossibleFalsePositive(testCall, confidenceRegionsSeq, lastRowCounts);
				testCallsSeq.removeFirst();
				continue;
			} else if(testCall.getFirst()>=clusterLast+DEF_MIN_CLUSTER_DISTANCE) break;
			testCallsCluster.add(testCall);
			regionFirst = Math.min(regionFirst, testCall.getFirst());
			regionLast = Math.max(regionLast, testCall.getLast());
			if(testCall.isHeterozygous()) heterozygous++;
			testCallsSeq.removeFirst();
		}
			
		
		//if(gsCalls.size()>0 && testCalls.size()>0) System.out.println("GS Calls: "+gsCalls.size()+" test calls: "+testCalls.size()+" first: "+clusterFirst+" first gs call: "+gsCalls.get(0).getFirst()+" first test call: "+testCalls.get(0).getFirst());
		while(confidenceRegionsSeq.size()>0) {
			GenomicRegion nextConfident = confidenceRegionsSeq.peek();
			if(nextConfident.getLast()<clusterFirst) confidenceRegionsSeq.removeFirst();
			else break;
		}
		boolean clusterInConfidenceRegion = false;
		for(GenomicRegion r:confidenceRegionsSeq) {
			if(r.getFirst()<=clusterFirst && r.getLast()>=clusterLast) clusterInConfidenceRegion = true;
			else if (r.getFirst()>clusterFirst) break;
		}
		if (confidenceRegions!=null && !clusterInConfidenceRegion) return;
		distClusterSizeGS.processDatapoint(gsCalls.size());
		distClusterTestHet.processDatapoint(heterozygous);
		distClusterSpan.processDatapoint(clusterLast-clusterFirst+1);
		if(testCallsCluster.size()==0) {
			processFalseNegativeCluster(gsCalls, clusterGSType, lastRowCounts);
		} else  {
			CalledGenomicVariant firstGS = gsCalls.get(0);
			CalledGenomicVariant firstTest = testCallsCluster.get(0);
			int genotypeFirstGS = getGenotypeNumber(firstGS);
			int genotypeFirstTest = getGenotypeNumber(firstTest);
			short qualFirstTest = loadGenotypeQuality(firstTest);
			GoldStandardComparisonCounts counts = countsPerType.get(clusterGSType);
			if(clusterGSType == GenomicVariant.TYPE_BIALLELIC_SNV) {
				//Process SNVs
				if(firstGS.getFirst()==firstTest.getFirst() && firstGS.getAlleles()[1].equals(firstTest.getAlleles()[1])) {
					// Match between gold standard and test SNV
					int k = 3*genotypeFirstGS+genotypeFirstTest;
					int row = Math.min(qualFirstTest/10, lastRowCounts); 
					counts.update(0,row,k);
					counts.update(row+1,lastRowCounts,9+genotypeFirstGS);
					if(mode == 3 && genotypeFirstGS!=genotypeFirstTest &&  qualFirstTest>=minQuality) {
						System.out.println("Variant "+firstGS.getSequenceName()+": "+firstGS.getFirst()+"genotypeGS: "+genotypeFirstGS+" genotypeTest: "+genotypeFirstTest+" alternativeGS: "+firstGS.getAlleles()[1]+" alternative Test: "+firstTest.getAlleles()[1]);
					}
				} else {
					//Isolated SNV calls in different close positions or with different alternative alleles
					processPossibleFalsePositive(firstTest, confidenceRegionsSeq, lastRowCounts);
					processFalseNegativeCluster(gsCalls, clusterGSType, lastRowCounts);
				}
			} else {
				regionFirst = Math.max(1, regionFirst);
				regionLast = Math.min(sequenceLength, regionLast);
				CharSequence refS = genome.getReference(firstGS.getSequenceName(), regionFirst, regionLast);
				if(refS==null) {
					System.err.println("WARN: Null reference sequence at "+firstGS.getSequenceName()+": "+regionFirst+"-"+regionLast+" cluster coordinates: "+clusterFirst+"-"+clusterLast+"cluster type: "+clusterGSType+" cluster size: "+gsCalls.size());
					return;
				}
				String reference = refS.toString().toUpperCase();
				GoldStandardHaplotypeReconstruction gsHaps = new GoldStandardHaplotypeReconstruction(reference, gsCalls, regionFirst);
				String [] gsHaplotypes = gsHaps.getPhasedAlleles();
				int genotypeGS = gsHaps.getGenotypeNumber();
				String [] matchingHaplotypes = gsHaps.buildMatchingHaplotypes(testCallsCluster);
				int genotypeTest = CalledGenomicVariant.GENOTYPE_HOMOALT;
				if(!matchingHaplotypes[0].equals(matchingHaplotypes[1])) genotypeTest = CalledGenomicVariant.GENOTYPE_HETERO;
				boolean sequenceMismatch = false;
				if(genotypeGS==genotypeTest) {
					//Record genotype errors if haplotype sequences do not match even if the genotypes match 
					if(genotypeGS==CalledGenomicVariant.GENOTYPE_HOMOALT && !gsHaplotypes[0].equals(matchingHaplotypes[0])) {
						genotypeTest = CalledGenomicVariant.GENOTYPE_HETERO;
						sequenceMismatch = true;
						//System.out.println("Changing genotype test to heterozygous due to error in sequence. GS hap: "+gsHaplotypes[0]+" test hap: "+matchingHaplotypes[0]);
					} else if (genotypeGS==CalledGenomicVariant.GENOTYPE_HETERO && (!gsHaplotypes[0].equals(matchingHaplotypes[0]) || (!gsHaplotypes[1].equals(matchingHaplotypes[1])))) {
						genotypeTest = CalledGenomicVariant.GENOTYPE_HOMOALT;
						sequenceMismatch = true;
					}
				}
				short qualTest = calculateQuality(testCallsCluster);
				int k = 3*genotypeGS+genotypeTest;
				int row = Math.min(qualTest/10, lastRowCounts); 
				counts.update(0,row,k);
				counts.update(row+1,lastRowCounts,9+genotypeGS);
				if(mode == 3 && genotypeGS!=genotypeTest &&  qualTest>=minQuality) {
					System.out.println("Variant "+gsCalls.get(0).getSequenceName()+": "+gsHaps.getFirst()+" type: "+clusterGSType+" genotypeGS: "+genotypeGS+" genotypeTest: "+genotypeTest+" sequenceMismatch: "+sequenceMismatch);
					System.out.println(gsHaplotypes[0]);
					System.out.println(gsHaplotypes[1]);
					System.out.println(matchingHaplotypes[0]);
					System.out.println(matchingHaplotypes[1]);
				}
			}
		}
	}

	

	private void processPossibleFalsePositive(CalledGenomicVariant call, LinkedList<GenomicRegion> confidenceRegionsSeq, int lastRowCounts) {
		boolean clusterInConfidenceRegion = false;
		for(GenomicRegion r:confidenceRegionsSeq) {
			if(r.getFirst()<=call.getFirst() && r.getLast()>=call.getLast()) clusterInConfidenceRegion = true;
			else if (r.getFirst()>call.getLast()) break;
		}
		int n = 0;
		if (confidenceRegions!=null && !clusterInConfidenceRegion) n=12;
		int genotypeTest = getGenotypeNumber(call);
		short qualTest = loadGenotypeQuality(call);
		byte typeTest = loadType(call); 
		
		countsPerType.get(typeTest).update(0,Math.min(qualTest/10, lastRowCounts),n+genotypeTest);
		if(mode == 2 && n==0 && qualTest>=minQuality) {
			System.out.println("Variant "+call.getSequenceName()+": "+call.getFirst()+" genotype: "+genotypeTest+" GQ: "+qualTest+" type: "+typeTest);
		}
	}

	private void processFalseNegativeCluster(List<CalledGenomicVariant> gsCalls, byte clusterType, int lastRowCounts) {
		int clusterGenotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
		for(CalledGenomicVariant call:gsCalls) {
			int genotypeGS = getGenotypeNumber(call);
			if(genotypeGS==CalledGenomicVariant.GENOTYPE_HETERO) {
				clusterGenotype = genotypeGS;
				break;
			}
		}
		countsPerType.get(clusterType).update(0,lastRowCounts,9+clusterGenotype);
		if(mode == 1) {
			System.out.println("False negative cluster with "+gsCalls.size()+" gold standard calls");
			for(CalledGenomicVariant call:gsCalls) {
				System.out.println("Variant "+call.getSequenceName()+": "+call.getFirst()+" genotype: "+getGenotypeNumber(call)+" type: "+loadType(call)+" alt allele: "+call.getAlleles()[1]);
			}
		}
		
	}
	
	public void printStatistics(PrintStream out) {
		out.println("Distribution of number of gold standard variants per cluster: ");
		distClusterSizeGS.printDistributionInt(out);
		out.println();
		out.println("Distribution of number of heterozygoys test variants per cluster: ");
		distClusterTestHet.printDistributionInt(out);
		out.println();
		out.println("Cluster span distribution: ");
		distClusterSpan.printDistributionInt(out);
		out.println();
		out.println("SNVs");
		countsPerType.get(GenomicVariant.TYPE_BIALLELIC_SNV).print(out);
		out.println("Indels");
		countsPerType.get(GenomicVariant.TYPE_INDEL).print(out);
		out.println("STRs/OTHER");
		countsPerType.get(GenomicVariant.TYPE_STR).print(out);
		
	}	
	private byte loadType(CalledGenomicVariant call) {
		byte type = call.getType();
		if(type == GenomicVariant.TYPE_BIALLELIC_SNV) return type;
		if(type == GenomicVariant.TYPE_INDEL) return type;
		if(type == GenomicVariant.TYPE_STR) return type;
		String [] alleles = call.getAlleles();
		if(alleles.length==2 && alleles[0].length()!=alleles[1].length()) return GenomicVariant.TYPE_INDEL;
		else if(alleles.length==2 && alleles[0].length()==1) return GenomicVariant.TYPE_BIALLELIC_SNV;
		return GenomicVariant.TYPE_STR;
	}

	private static int getGenotypeNumber(CalledGenomicVariant call) {
		if(call.isHeterozygous()) return CalledGenomicVariant.GENOTYPE_HETERO;
		else if (!call.isHomozygousReference()) return CalledGenomicVariant.GENOTYPE_HOMOALT;
		return CalledGenomicVariant.GENOTYPE_HOMOREF;
	}
	private short loadGenotypeQuality(CalledGenomicVariant call) {
		short q = call.getGenotypeQuality();
		if(q>0) return q;
		String [] alleles = call.getAlleles();
		String [] calledAlleles = call.getCalledAlleles();
		VariantCallReport report = call.getCallReport();
		if(report==null || calledAlleles.length==0 || !report.logConditionalsPresent()) return 0;
		double logP;
		Double sum=null;
		if(calledAlleles.length==1) logP = report.getLogConditionalProbability(calledAlleles[0], calledAlleles[0]);
		else logP = report.getLogConditionalProbability(calledAlleles[0], calledAlleles[1]);
		for(int i=0;i<alleles.length;i++) {
			for(int j=i;j<alleles.length;j++) {
				sum = LogMath.logSum(sum, report.getLogConditionalProbability(alleles[i], alleles[j]));
			}
		}
		double logPos = LogMath.logProduct(logP, -sum);
		double pos = LogMath.power10(logPos);
		q = PhredScoreHelper.calculatePhredScore(1-pos);
		call.setGenotypeQuality(q);
		//if(call.getFirst()==376) System.out.println("Conditional P: "+logP+" sum: "+sum+" logpos: "+logPos+" posterior: "+pos+" q: "+q);
		return q;
	}
	private short calculateQuality(List<CalledGenomicVariant> testCalls) {
		short quality = 255;
		for(CalledGenomicVariant call:testCalls) {
			short q = loadGenotypeQuality(call);
			if(q<quality) quality = q;
		}
		return quality;
	}
}
class GoldStandardComparisonCounts {
	public static final int NUM_ROWS_COUNTS = 10;
	private int [][] counts;
	private static final DecimalFormat DF = new DecimalFormat("0.0000");
	private boolean countNonGSAsFP = false;
	private double confidentMbp=3000;
	public GoldStandardComparisonCounts () {
		counts = new int [NUM_ROWS_COUNTS][15];
	}
	
	
	/**
	 * @return the countNonGSAsFP
	 */
	public boolean isCountNonGSAsFP() {
		return countNonGSAsFP;
	}


	/**
	 * @param countNonGSAsFP the countNonGSAsFP to set
	 */
	public void setCountNonGSAsFP(boolean countNonGSAsFP) {
		this.countNonGSAsFP = countNonGSAsFP;
	}
	
	/**
	 * @return the confidentMbp
	 */
	public double getConfidentMbp() {
		return confidentMbp;
	}


	/**
	 * @param confidentMbp the confidentMbp to set
	 */
	public void setConfidentMbp(double confidentMbp) {
		this.confidentMbp = confidentMbp;
	}


	public void update(int firstRow, int lastRow, int column) {
		for(int i=firstRow;i<=lastRow;i++) {
			counts[i][column]++;
		}	
	}
	
	public void print(PrintStream out) {
		for (int i=0;i<counts.length; i++) {
			out.print(""+(i*10));
			for(int j=0;j<counts[0].length;j++) {
				out.print("\t"+counts[i][j]);
			}
			printComparisonStats(counts[i], out);
			out.println();
		}
	}

	private void printComparisonStats(int[] row, PrintStream out) {
		int gsTotal0 = row[0]+row[1]+row[2]+row[9];
		int gsTotal1 = row[3]+row[4]+row[5]+row[10];
		int gsTotal2 = row[6]+row[7]+row[8]+row[11];
		out.print("\t"+gsTotal0+"\t"+gsTotal1+"\t"+gsTotal2);
		
		int testTotal0 = row[0]+row[3]+row[6]+row[12];
		int testTotal1 = row[1]+row[4]+row[7]+row[13];
		int testTotal2 = row[2]+row[5]+row[8]+row[14];
		out.print("\t"+testTotal0+"\t"+testTotal1+"\t"+testTotal2);
		
		double recall1 = 0;
		if(gsTotal1>0) {
			recall1 = (double)row[4] / gsTotal1;
		}
		int fd1 = row[1] + row [7];
		double denomFDR1 = testTotal1;
		if(countNonGSAsFP) fd1 += row[13];
		else denomFDR1 -= row[13];
		double fppm1 = fd1/confidentMbp;
		double fdr1 = 0;
		double precision1 = 1;
		if(denomFDR1>1) {
			fdr1 = (double)fd1 / denomFDR1;
			precision1 = (double)row[4] / denomFDR1;
		}
		double f1 = 0;
		if(precision1 + recall1 > 0) {
			f1 = 2.0*precision1*recall1/(precision1+recall1);
		}
		
		out.print("\t"+DF.format(recall1)+"\t"+fd1+"\t"+DF.format(fppm1)+"\t"+DF.format(fdr1)+"\t"+DF.format(precision1)+"\t"+DF.format(f1));
		
		double recall2 = 0;
		if(gsTotal2>0) {
			recall2 = (double)row[8] / gsTotal2;
		}
		int fd2 = row[2] + row [5];
		double denomFDR2 = testTotal2;
		if(countNonGSAsFP) fd2 += row[14];
		else denomFDR2 -= row[14];
		double fppm2 = fd2/confidentMbp;
		double fdr2 = 0;
		double precision2 = 1;
		if(denomFDR2>0) {
			fdr2 = (double)fd2 / denomFDR2;
			precision2 = (double)row[8] / denomFDR2;
		}
		double f2 = 0;
		if(precision2 + recall2 > 0) {
			f2 = 2.0*precision2*recall2/(precision2+recall2);
		}
		out.print("\t"+DF.format(recall2)+"\t"+fd2+"\t"+DF.format(fppm2)+"\t"+DF.format(fdr2)+"\t"+DF.format(precision2)+"\t"+DF.format(f2));
		
	}
	
}
class GoldStandardHaplotypeReconstruction implements CalledGenomicVariant {
	private List<CalledGenomicVariant> calls;
	private String sequenceName;
	private String reference;
	private String haplotype0;
	private String haplotype1;
	private int numAlleles = 1;
	private byte genotypeN = 0;
	private int first;
	private int last;
	private Map<Integer, Integer> startsH0;
	private Map<Integer, Integer> startsH1;
	private int posPrint = -1;
	
	public GoldStandardHaplotypeReconstruction (String reference, List<CalledGenomicVariant> calls, int first) {
		this.calls = calls;
		this.first = first;
		this.last = first+reference.length()-1;
		this.reference = reference;
		buildHaplotypes();
		boolean b0 = haplotype0.equals(reference);
		boolean b1 = haplotype1.equals(reference);
		boolean b3 = haplotype0.equals(haplotype1);
		
		if(!b0) numAlleles++;
		if(!b1) {
			if(!b3) numAlleles++;
		}
		genotypeN = CalledGenomicVariant.GENOTYPE_HOMOREF;
		if(!b3) genotypeN = CalledGenomicVariant.GENOTYPE_HETERO;
		else if (!b0 || !b1) genotypeN = CalledGenomicVariant.GENOTYPE_HOMOALT; 
		
	}
	
	private void buildHaplotypes() {
		StringBuilder hap0 = new StringBuilder();
		StringBuilder hap1 = new StringBuilder();
		startsH0 = new HashMap<>();
		startsH1 = new HashMap<>();
		if(first==posPrint) System.out.println("Cluster first: "+first+" Gold standard calls: "+calls.size());
		int pos = first;
		for(CalledGenomicVariant call:calls) {
			if(first==posPrint) System.out.println("Cluster first: "+first+" Gold standard call: "+call.getSequenceName()+": "+call.getFirst()+" heterozygous: "+call.isHeterozygous()+" phased: "+call.isPhased());
			String [] calledAlleles = call.getCalledAlleles();
			String [] phasedAlleles = call.getPhasedAlleles();
			//System.out.println("Pos: "+pos+" first call: "+call.getFirst()+" first this: "+first);
			if(call.getFirst()>pos) {
				String refString = reference.substring(pos-first, call.getFirst()-first); 
				hap0.append(refString);
				hap1.append(refString);
			}
			startsH0.put(call.getFirst(), hap0.length());
			startsH1.put(call.getFirst(), hap1.length());
			if(calledAlleles.length==1) {
				hap0.append(calledAlleles[0]);
				hap1.append(calledAlleles[0]);
			} else if(phasedAlleles!=null && phasedAlleles.length==2) {
				hap0.append(phasedAlleles[0]);
				hap1.append(phasedAlleles[1]);
			} else {
				System.err.println("WARN: Unphased heterozygous gold standard call at "+call.getSequenceName()+": "+call.getFirst());
				hap0.append(calledAlleles[0]);
				hap1.append(calledAlleles[1]);
			}
			pos = call.getLast()+1;
		}
		if(pos<=last) {
			String refEnd = reference.substring(pos-first);
			hap0.append(refEnd);
			hap1.append(refEnd);
		}
		haplotype0 = hap0.toString();
		haplotype1 = hap1.toString();	
	}

	public String [] buildMatchingHaplotypes(List<CalledGenomicVariant> testCalls) {
		
		//Count heterozygous calls
		int hetCalls = 0;
		for(CalledGenomicVariant call:testCalls) {
			if(call.isHeterozygous()) hetCalls++;
		}
		if(hetCalls<=8) return exhaustiveMatchingHaplotypes (testCalls,hetCalls);
		return greedyMatchingHaplotypes (testCalls);
		
	}
	
	private String[] greedyMatchingHaplotypes(List<CalledGenomicVariant> testCalls) {
		StringBuilder hap0 = new StringBuilder();
		StringBuilder hap1 = new StringBuilder();
		int pos = first;
		for(CalledGenomicVariant call:testCalls) {
			String [] calledAlleles = call.getCalledAlleles();
			if(call.getFirst()>pos) {
				String refString = reference.substring(pos-first, call.getFirst()-first); 
				hap0.append(refString);
				hap1.append(refString);
			}
			if(calledAlleles.length==1) {
				//Homozygous alternative
				hap0.append(calledAlleles[0]);
				hap1.append(calledAlleles[0]);
			} else {
				HammingSequenceDistanceMeasure measure = new HammingSequenceDistanceMeasure();
				String search00 = hap0+calledAlleles[0];
				String search11 = hap1+calledAlleles[1];
				
				String search01 = hap0+calledAlleles[1];
				String search10 = hap1+calledAlleles[0];
				
				double d0 = 0;
				if(haplotype0.length()>=search00.length()) {
					String subject00 = haplotype0.substring(0,search00.length());
					d0 += measure.calculateDistance(search00, subject00);
				} else {
					d0+=2*(search00.length()-haplotype0.length());
				}
				
				if(haplotype1.length()>=search11.length()) {
					String subject11 = haplotype1.substring(0,search11.length());
					d0+=measure.calculateDistance(search11, subject11);
				} else {
					d0+=2*(search11.length()-haplotype1.length());
				}
				double d1=0;
				if(haplotype0.length()>=search01.length()) {
					String subject01 = haplotype0.substring(0,search01.length());
					d1 += measure.calculateDistance(search01, subject01);
				} else {
					d1+=2*(search01.length()-haplotype0.length());
				}
				
				if(haplotype1.length()>=search10.length()) {
					String subject10 = haplotype1.substring(0,search10.length());
					d1+=measure.calculateDistance(search10, subject10);
				} else {
					d1+=2*(search10.length()-haplotype1.length());
				}
				if(first == posPrint) {
					System.out.println("Comparing haplotypes: ");
					System.out.println(haplotype0);
					System.out.println(search00);
					System.out.println(search01);
					System.out.println();
					System.out.println(haplotype1);
					System.out.println(search11);
					System.out.println(search10);
				}
				if(d0<d1) {
					hap0.append(calledAlleles[0]);
					hap1.append(calledAlleles[1]);
				} else {
					hap0.append(calledAlleles[1]);
					hap1.append(calledAlleles[0]);
				}
			}
			pos = call.getLast()+1;
		}
		if(pos<=last) {
			String refEnd = reference.substring(pos-first);
			hap0.append(refEnd);
			hap1.append(refEnd);
		}
		String [] answer = {hap0.toString(),hap1.toString()}; 
		return answer;
	}

	private String[] exhaustiveMatchingHaplotypes(List<CalledGenomicVariant> testCalls, int hetCalls) {
		if(first == posPrint) System.out.println("Exhaustive matching test haplotypes. Total calls: "+testCalls.size()+" heterozygous calls: "+hetCalls);
		StringBuilder [][] haplotypePairs = new StringBuilder[(int) Math.pow(2, hetCalls)][2];
		for(int i=0;i<haplotypePairs.length;i++) {
			haplotypePairs[i][0] = new StringBuilder(2*haplotype0.length());
			haplotypePairs[i][1] = new StringBuilder(2*haplotype1.length());
		}
		int pos = first;
		int power = 1;
		for(CalledGenomicVariant call:testCalls) {
			String [] calledAlleles = call.getCalledAlleles();
			if(call.getFirst()>pos) {
				if(call.getFirst() <first) throw new RuntimeException("Test call at pos: "+call.getFirst()+" before cluster start: "+first);
				String refString = reference.substring(pos-first, call.getFirst()-first);
				for(int i=0;i<haplotypePairs.length;i++) {
					haplotypePairs[i][0].append(refString);
					haplotypePairs[i][1].append(refString);
				}
				
			}
			if(calledAlleles.length==1) {
				//Homozygous
				for(int i=0;i<haplotypePairs.length;i++) {
					haplotypePairs[i][0].append(calledAlleles[0]);
					haplotypePairs[i][1].append(calledAlleles[0]);
				}
			} else {
				for(int i=0;i<haplotypePairs.length;i++) {
					boolean allele1First = i%(2*power)>=power;
					if(allele1First) {
						haplotypePairs[i][0].append(calledAlleles[1]);
						haplotypePairs[i][1].append(calledAlleles[0]);
					} else {
						haplotypePairs[i][0].append(calledAlleles[0]);
						haplotypePairs[i][1].append(calledAlleles[1]);
					}
				}
				power*=2;
			}
			pos = call.getLast()+1;
		}
		if(pos<=last) {
			String refEnd = reference.substring(pos-first);
			for(int i=0;i<haplotypePairs.length;i++) {
				haplotypePairs[i][0].append(refEnd);
				haplotypePairs[i][1].append(refEnd);
			}
		}
		//Pick best haplotype pair
		String [] answer = new String [2];
		double minD = 0;
		if(first == posPrint) {
			System.out.println("Find best match");
			System.out.println(haplotype0);
			System.out.println(haplotype1);
		}
		for(int i=0;i<haplotypePairs.length;i++) {
			HammingSequenceDistanceMeasure measure = new HammingSequenceDistanceMeasure();
			double d=0;
			String hapTest0 = haplotypePairs[i][0].toString();
			String hapTest1 = haplotypePairs[i][1].toString();
			if(haplotype0.length()==hapTest0.length()) {
				d += measure.calculateDistance(hapTest0, haplotype0);
			} else if(haplotype0.length()>hapTest0.length()) {
				String subject = haplotype0.substring(0,hapTest0.length());
				d += measure.calculateDistance(hapTest0, subject)+2*(haplotype0.length()-hapTest0.length());
			} else {
				String query = hapTest0.substring(0,haplotype0.length());
				d += measure.calculateDistance(query, haplotype0)+2*(hapTest0.length()-haplotype0.length());
			}
			if(haplotype1.length()==hapTest1.length()) {
				d += measure.calculateDistance(hapTest1, haplotype1);
			} else if(haplotype1.length()>hapTest1.length()) {
				String subject = haplotype1.substring(0,hapTest1.length());
				d += measure.calculateDistance(hapTest1, subject)+2*(haplotype1.length()-hapTest1.length());
			} else {
				String query = hapTest1.substring(0,haplotype1.length());
				d += measure.calculateDistance(query, haplotype1)+2*(hapTest1.length()-haplotype1.length());
			}
			if (first == posPrint) System.out.println(hapTest0+" "+d);
			if(answer[0]==null || minD>d) {
				answer[0] = hapTest0;
				answer[1] = hapTest1;
				minD = d;
			}
		}
		
		
		return answer;
	}

	/**
	 * @return the haplotype0
	 */
	public String getHaplotype0() {
		return haplotype0;
	}

	/**
	 * @return the haplotype1
	 */
	public String getHaplotype1() {
		return haplotype1;
	}
	
	public int getGenotypeNumber () {
		return genotypeN;
	}

	/**
	 * @return the startsH0
	 */
	public Integer getStartH0(int referencePos) {
		return startsH0.get(referencePos);
	}

	/**
	 * @return the startsH1
	 */
	public Integer getStartsH1(int referencePos) {
		return startsH1.get(referencePos);
	}

	@Override
	public String getSequenceName() {
		return sequenceName;
	}

	@Override
	public int getFirst() {
		return first;
	}

	@Override
	public int getLast() {
		return last;
	}

	@Override
	public int length() {
		return reference.length();
	}

	@Override
	public boolean isPositiveStrand() {
		return true;
	}

	@Override
	public boolean isNegativeStrand() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public String[] getAlleles() {
		Set <String> allelesSet = new HashSet<>();
		allelesSet.add(reference);
		allelesSet.add(haplotype0);
		allelesSet.add(haplotype1);
		return allelesSet.toArray(new String [0]);
	}

	@Override
	public String getReference() {
		return reference;
	}

	@Override
	public String getId() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setId(String id) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public short getVariantQS() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setVariantQS(short qualityScore) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean isCompatible(GenomicVariant variant) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isBiallelic() {
		return numAlleles==2;
	}

	@Override
	public boolean isSNV() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public byte getType() {
		// TODO Auto-generated method stub
		return GenomicVariant.TYPE_STR;
	}

	@Override
	public void setType(byte type) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String getSampleId() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setSampleId(String sampleId) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String[] getCalledAlleles() {
		Set <String> allelesSet = new HashSet<>();
		allelesSet.add(haplotype0);
		allelesSet.add(haplotype1);
		return allelesSet.toArray(new String [0]);
	}

	@Override
	public byte[] getIndexesCalledAlleles() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public short[] getAllelesCopyNumber() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public short getCopyNumber() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getTotalReadDepth() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setTotalReadDepth(int depth) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public short getGenotypeQuality() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setGenotypeQuality(short genotypeQuality) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void makeUndecided() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean isUndecided() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isHeterozygous() {
		return genotypeN == 1;
	}

	@Override
	public boolean isHomozygous() {
		return genotypeN != 1;
	}

	@Override
	public boolean isHomozygousReference() {
		return genotypeN == 0;
	}

	@Override
	public void updateAllelesCopyNumberFromCounts(short totalCopyNumber) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setAllelesCopyNumber(short[] allelesCN) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public VariantCallReport getCallReport() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int[] getAllCounts() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isPhased() {
		return true;
	}

	@Override
	public String[] getPhasedAlleles() {
		String [] haplotypes = new String[2];
		haplotypes[0] = haplotype0;
		haplotypes[1] = haplotype1;
		return haplotypes;
	}

	@Override
	public byte[] getIndexesPhasedAlleles() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public byte getStrandBiasScore() {
		// TODO Auto-generated method stub
		return 0;
	}
	
}
