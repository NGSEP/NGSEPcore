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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.discovery.rd.ReadDepthBin;
import ngsep.discovery.rd.ReadDepthDistribution;
import ngsep.discovery.rd.SingleSampleReadDepthAlgorithm;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.GenomicRegionSpanComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.AbstractLimitedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledCNV;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.Sample;
import ngsep.variants.io.GFFVariantsFileHandler;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;


public class SingleSampleVariantsDetector implements PileupListener {
	
	// Constants for default values
	public static final int DEF_MAX_ALNS_PER_START_POS = AlignmentsPileupGenerator.DEF_MAX_ALNS_PER_START_POS;
	public static final double DEF_MIN_ALLELE_FREQUENCY = 0;
	public static final double DEF_HETEROZYGOSITY_RATE_DIPLOID = SingleSampleVariantPileupListener.DEF_HETEROZYGOSITY_RATE_DIPLOID;
	public static final short DEF_MIN_QUALITY = SingleSampleVariantPileupListener.DEF_MIN_QUALITY;
	public static final short DEF_MIN_MQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	public static final byte DEF_MAX_BASE_QS = SingleSampleVariantPileupListener.DEF_MAX_BASE_QS;
	public static final byte DEF_PLOIDY = GenomicVariant.DEFAULT_PLOIDY;
	public static final String DEF_SAMPLE_ID = SingleSampleVariantPileupListener.DEF_SAMPLE_ID;
	public static final short DEF_MIN_SV_QUALITY = 20;
	public static final short DEF_BIN_SIZE = ReadDepthDistribution.DEFAULT_BIN_SIZE;
	public static final String DEF_ALGORITHM_CNV = "CNVnator";
	public static final short DEF_MAX_PCT_OVERLAP_CNVS = 100;
	public static final int DEF_MAX_LEN_DELETION = ReadPairAnalyzer.DEF_MAX_LEN_DELETION;
	public static final int DEF_SPLIT_READ_SEED = ReadPairAnalyzer.DEF_SPLIT_READ_SEED;
	
	// Logging and progress
	private Logger log = Logger.getLogger(SingleSampleVariantsDetector.class.getName());
	private ProgressNotifier progressNotifier = null;
	private double coveredGenomeSize = 0;
	private long referenceGenomeSize = 0;
	
	// Parameters
	private String inputFile=null;
	private ReferenceGenome genome;
	private String outputPrefix=null;
	private String knownVariantsFile=null;
	private String knownSTRsFile=null;
	private short normalPloidy = DEF_PLOIDY;
	private boolean printSamplePloidy = false;
	private String sampleId = DEF_SAMPLE_ID;
	private String knownSVsFile=null;
	private long inputGenomeSize = 0;
	private int binSize = ReadDepthDistribution.DEFAULT_BIN_SIZE;
	private String algCNV = DEF_ALGORITHM_CNV;
	private short minSVQuality = DEF_MIN_SV_QUALITY;
	private int maxPCTOverlapCNVs = DEF_MAX_PCT_OVERLAP_CNVS;
	private boolean findRepeats = false;
	private boolean runRDAnalysis = false;
	private boolean findSNVs = true;
	private boolean runRPAnalysis = false;
	private boolean findNewCNVs = true;
	// Classes implementing the algorithms for structural variants detection
	private MultipleMappingRegionsCalculator mmRegsCalc = new MultipleMappingRegionsCalculator();
	private ReadPairAnalyzer rpAnalyzer = new ReadPairAnalyzer();
	//Listeners
	private IndelRealignerPileupListener indelRealigner = new IndelRealignerPileupListener();
	private SingleSampleVariantPileupListener varListener = new SingleSampleVariantPileupListener();
	private boolean hetRateModified = false; 
	
	// Model attributes
	private AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator(); 
	private GenomicRegionSortedCollection<CalledGenomicVariant> calledSVs;
	
	// File handlers
	private VCFFileWriter varsFW = new VCFFileWriter();
	private VCFFileHeader header;
	private GFFVariantsFileHandler svsFH = new GFFVariantsFileHandler();
	
	
	
	//Objects for output files
	private PrintStream outVars = null;

	// Get and set methods
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
		generator.setLog(log);
		rpAnalyzer.setLog(log);
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
	
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	public void setGenome(String genomeFile) throws IOException {
		setGenome(OptionValuesDecoder.loadGenome(genomeFile,log));
	}

	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	
	public String getSampleId() {
		return sampleId;
	}
	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}
	
	public short getNormalPloidy() {
		return normalPloidy;
	}
	public void setNormalPloidy(short normalPloidy) {
		this.normalPloidy = normalPloidy;
	}
	public void setNormalPloidy(String value) {
		setNormalPloidy((short)OptionValuesDecoder.decode(value, Short.class));
	}
	
	public boolean isPrintSamplePloidy() {
		return printSamplePloidy;
	}
	public void setPrintSamplePloidy(boolean printSamplePloidy) {
		this.printSamplePloidy = printSamplePloidy;
	}
	public void setPrintSamplePloidy(Boolean printSamplePloidy) {
		this.setPrintSamplePloidy(printSamplePloidy.booleanValue());
	}
	
	/**
	 * @return int
	 * @see ngsep.discovery.AlignmentsPileupGenerator#getMinMQ()
	 */
	public int getMinMQ() {
		return generator.getMinMQ();
	}
	/**
	 * @param minMQ
	 * @see ngsep.discovery.AlignmentsPileupGenerator#setMinMQ(int)
	 */
	public void setMinMQ(int minMQ) {
		mmRegsCalc.setMinMQ(minMQ);
		generator.setMinMQ(minMQ);
		rpAnalyzer.setMinMQ(minMQ);
	}
	public void setMinMQ(String minMQ) {
		setMinMQ((int)OptionValuesDecoder.decode(minMQ, Integer.class));
	}
	
	public String getKnownVariantsFile() {
		return knownVariantsFile;
	}
	public void setKnownVariantsFile(String knownVariantsFile) {
		this.knownVariantsFile = knownVariantsFile;
	}
	
	/**
	 * @return String
	 * @see ngsep.discovery.AlignmentsPileupGenerator#getQuerySeq()
	 */
	public String getQuerySeq() {
		return generator.getQuerySeq();
	}
	/**
	 * @param querySeq
	 * @see ngsep.discovery.AlignmentsPileupGenerator#setQuerySeq(java.lang.String)
	 */
	public void setQuerySeq(String querySeq) {
		generator.setQuerySeq(querySeq);
	}

	/**
	 * @return int
	 * @see ngsep.discovery.AlignmentsPileupGenerator#getQueryFirst()
	 */
	public int getQueryFirst() {
		return generator.getQueryFirst();
	}

	/**
	 * @param queryFirst
	 * @see ngsep.discovery.AlignmentsPileupGenerator#setQueryFirst(int)
	 */
	public void setQueryFirst(int queryFirst) {
		generator.setQueryFirst(queryFirst);
	}
	public void setQueryFirst(String value) {
		setQueryFirst((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	/**
	 * @return int
	 * @see ngsep.discovery.AlignmentsPileupGenerator#getQueryLast()
	 */
	public int getQueryLast() {
		return generator.getQueryLast();
	}
	/**
	 * @param queryLast
	 * @see ngsep.discovery.AlignmentsPileupGenerator#setQueryLast(int)
	 */
	public void setQueryLast(int queryLast) {
		generator.setQueryLast(queryLast);
	}
	public void setQueryLast(String value) {
		setQueryLast((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public boolean isIgnoreLowerCaseRef() {
		return varListener.isIgnoreLowerCaseRef();
	}
	public void setIgnoreLowerCaseRef(boolean ignoreLowerCaseRef) {
		varListener.setIgnoreLowerCaseRef(ignoreLowerCaseRef);
	}
	public void setIgnoreLowerCaseRef(Boolean ignoreLowerCaseRef) {
		setIgnoreLowerCaseRef(ignoreLowerCaseRef.booleanValue());
	}
	
	/**
	 * @return int
	 * @see ngsep.discovery.AlignmentsPileupGenerator#getMaxAlnsPerStartPos()
	 */
	public int getMaxAlnsPerStartPos() {
		return generator.getMaxAlnsPerStartPos();
	}
	/**
	 * @param maxAlnsPerStartPos
	 * @see ngsep.discovery.AlignmentsPileupGenerator#setMaxAlnsPerStartPos(int)
	 */
	public void setMaxAlnsPerStartPos(int maxAlnsPerStartPos) {
		generator.setMaxAlnsPerStartPos(maxAlnsPerStartPos);
	}
	public void setMaxAlnsPerStartPos(String value) {
		setMaxAlnsPerStartPos((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	/**
	 * @return boolean
	 * @see ngsep.discovery.AlignmentsPileupGenerator#isProcessNonUniquePrimaryAlignments()
	 */
	public boolean isProcessNonUniquePrimaryAlignments() {
		return generator.isProcessNonUniquePrimaryAlignments();
	}
	/**
	 * @param processNonUniquePrimaryAlignments
	 * @see ngsep.discovery.AlignmentsPileupGenerator#setProcessNonUniquePrimaryAlignments(boolean)
	 */
	public void setProcessNonUniquePrimaryAlignments(boolean processNonUniquePrimaryAlignments) {
		generator.setProcessNonUniquePrimaryAlignments(processNonUniquePrimaryAlignments);
	}
	public void setProcessNonUniquePrimaryAlignments(Boolean processNonUniquePrimaryAlignments) {
		setProcessNonUniquePrimaryAlignments(processNonUniquePrimaryAlignments.booleanValue());
	}
	
	/**
	 * @return boolean
	 * @see ngsep.discovery.AlignmentsPileupGenerator#isProcessSecondaryAlignments()
	 */
	public boolean isProcessSecondaryAlignments() {
		return generator.isProcessSecondaryAlignments();
	}
	/**
	 * @param processSecondaryAlignments
	 * @see ngsep.discovery.AlignmentsPileupGenerator#setProcessSecondaryAlignments(boolean)
	 */
	public void setProcessSecondaryAlignments(boolean processSecondaryAlignments) {
		generator.setProcessSecondaryAlignments(processSecondaryAlignments);
	}
	public void setProcessSecondaryAlignments(Boolean processSecondaryAlignments) {
		setProcessSecondaryAlignments(processSecondaryAlignments.booleanValue());
	}
	/**
	 * @return
	 * @see ngsep.discovery.AlignmentsPileupGenerator#getBasesToIgnore5P()
	 */
	public byte getBasesToIgnore5P() {
		return generator.getBasesToIgnore5P();
	}
	/**
	 * @param basesToIgnore5P
	 * @see ngsep.discovery.AlignmentsPileupGenerator#setBasesToIgnore5P(byte)
	 */
	public void setBasesToIgnore5P(byte basesToIgnore5P) {
		generator.setBasesToIgnore5P(basesToIgnore5P);
	}
	public void setBasesToIgnore5P(String value) {
		setBasesToIgnore5P((byte)OptionValuesDecoder.decode(value, Byte.class));
	}
	/**
	 * @return
	 * @see ngsep.discovery.AlignmentsPileupGenerator#getBasesToIgnore3P()
	 */
	public byte getBasesToIgnore3P() {
		return generator.getBasesToIgnore3P();
	}
	/**
	 * @param basesToIgnore3P
	 * @see ngsep.discovery.AlignmentsPileupGenerator#setBasesToIgnore3P(byte)
	 */
	public void setBasesToIgnore3P(byte basesToIgnore3P) {
		generator.setBasesToIgnore3P(basesToIgnore3P);
	}
	public void setBasesToIgnore3P(String value) {
		setBasesToIgnore3P((byte)OptionValuesDecoder.decode(value, Byte.class));
	}

	public double getHeterozygosityRate() {
		return varListener.getHeterozygosityRate();
	}
	public void setHeterozygosityRate(double heterozygosityRate) {
		varListener.setHeterozygosityRate(heterozygosityRate);
		hetRateModified = true;
	}	
	public void setHeterozygosityRate(String value) {
		setHeterozygosityRate((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	public byte getMaxBaseQS() {
		return varListener.getMaxBaseQS();
	}
	public void setMaxBaseQS(byte maxBaseQS) {
		varListener.setMaxBaseQS(maxBaseQS);
	}
	public void setMaxBaseQS(String value) {
		setMaxBaseQS((byte)OptionValuesDecoder.decode(value, Byte.class));
	}
	
	public String getKnownSTRsFile() {
		return knownSTRsFile;
	}
	public void setKnownSTRsFile(String knownSTRsFile) {
		this.knownSTRsFile = knownSTRsFile;
	}
	
	public short getMinQuality() {
		return varListener.getMinQuality();
	}
	public void setMinQuality(short minQuality) {
		varListener.setMinQuality(minQuality);
	}	
	public void setMinQuality(String value) {
		setMinQuality((short)OptionValuesDecoder.decode(value, Short.class));
	}

	public boolean isCallEmbeddedSNVs() {
		return varListener.isCallEmbeddedSNVs();
	}
	public void setCallEmbeddedSNVs(boolean callEmbeddedSNVs) {
		varListener.setCallEmbeddedSNVs(callEmbeddedSNVs);
	}	
	public void setCallEmbeddedSNVs(Boolean callEmbeddedSNVs) {
		setCallEmbeddedSNVs(callEmbeddedSNVs.booleanValue());
	}
	
	public boolean isCalcStrandBias () {
		return varListener.isCalcStrandBias();
	}
	public void setCalcStrandBias(boolean calcStrandBias) {
		varListener.setCalcStrandBias(calcStrandBias);
		
	}
	public void setCalcStrandBias(Boolean calcStrandBias) {
		setCalcStrandBias(calcStrandBias.booleanValue());
	}

	public String getKnownSVsFile() {
		return knownSVsFile;
	}
	public void setKnownSVsFile(String knownSVsFile) {
		this.knownSVsFile = knownSVsFile;
	}
	
	public short getMinSVQuality() {
		return minSVQuality;
	}
	public void setMinSVQuality(short minSVQuality) {
		this.minSVQuality = minSVQuality;
	}
	public void setMinSVQuality(String value) {
		setMinSVQuality((short)OptionValuesDecoder.decode(value, Short.class));
	}
	
	public boolean isFindRepeats() {
		return findRepeats;
	}
	public void setFindRepeats(boolean findRepeats) {
		this.findRepeats = findRepeats;
	}
	public void setFindRepeats(Boolean findRepeats) {
		setFindRepeats(findRepeats.booleanValue());
	}

	public boolean isRunRDAnalysis() {
		return runRDAnalysis;
	}
	public void setRunRDAnalysis(boolean runRDAnalysis) {
		this.runRDAnalysis = runRDAnalysis;
	}
	public void setRunRDAnalysis(Boolean runRDAnalysis) {
		setRunRDAnalysis(runRDAnalysis.booleanValue());
	}
	
	public boolean isFindNewCNVs() {
		return findNewCNVs;
	}
	public void setFindNewCNVs(boolean findNewCNVs) {
		this.findNewCNVs = findNewCNVs;
	}
	public void setFindNewCNVs(Boolean findNewCNVs) {
		this.setFindNewCNVs(findNewCNVs.booleanValue());
	}

	public long getInputGenomeSize() {
		return inputGenomeSize;
	}
	public void setInputGenomeSize(long inputGenomeSize) {
		this.inputGenomeSize = inputGenomeSize;
	}
	public void setInputGenomeSize(String value) {
		setInputGenomeSize((long)OptionValuesDecoder.decode(value, Long.class));
	}
	
	public int getBinSize() {
		return binSize;
	}
	public void setBinSize(int binSize) {
		this.binSize = binSize;
	}
	public void setBinSize(String value) {
		setBinSize((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	public String getAlgCNV() {
		return algCNV;
	}
	public void setAlgCNV(String algCNV) {
		this.algCNV = algCNV;
	}

	public int getMaxPCTOverlapCNVs() {
		return maxPCTOverlapCNVs;
	}
	public void setMaxPCTOverlapCNVs(int maxPCTOverlapCNVs) {
		this.maxPCTOverlapCNVs = maxPCTOverlapCNVs;
	}
	public void setMaxPCTOverlapCNVs(String value) {
		setMaxPCTOverlapCNVs((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public boolean isRunRPAnalysis() {
		return runRPAnalysis;
	}
	public void setRunRPAnalysis(boolean runRPAnalysis) {
		this.runRPAnalysis = runRPAnalysis;
	}
	public void setRunRPAnalysis(Boolean runRPAnalysis) {
		setRunRPAnalysis(runRPAnalysis.booleanValue());
	}
	
	public int getMaxLengthDeletion() {
		return rpAnalyzer.getMaxLengthDeletion();
	}
	public void setMaxLengthDeletion(int maxLengthDeletion) {
		rpAnalyzer.setMaxLengthDeletion(maxLengthDeletion);
	}
	public void setMaxLengthDeletion(String value) {
		setMaxLengthDeletion((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getSplitReadSeed() {
		return rpAnalyzer.getSeedSize();
	}
	public void setSplitReadSeed(int seedSize) {
		rpAnalyzer.setSeedSize(seedSize);
	}
	public void setSplitReadSeed(String value) {
		setSplitReadSeed((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public boolean isIgnoreProperPairFlag() {
		return rpAnalyzer.isIgnoreProperPairFlag();
	}
	public void setIgnoreProperPairFlag(boolean ignoreProperPairFlag) {
		rpAnalyzer.setIgnoreProperPairFlag(ignoreProperPairFlag);
	}
	public void setIgnoreProperPairFlag(Boolean ignoreProperPairFlag) {
		setIgnoreProperPairFlag(ignoreProperPairFlag.booleanValue());
	}
	
	public boolean isRunOnlySVsAnalyses() {
		return !findSNVs;
	}
	public void setRunOnlySVsAnalyses(boolean runOnlySVsAnalyses) {
		this.findSNVs = !runOnlySVsAnalyses;
	}
	public void setRunOnlySVsAnalyses(Boolean runOnlySVsAnalyses) {
		setRunOnlySVsAnalyses(runOnlySVsAnalyses.booleanValue());
	}
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		SingleSampleVariantsDetector instance = new SingleSampleVariantsDetector();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}

	public void run () throws IOException {
		if(!runRDAnalysis || knownSVsFile!=null) findNewCNVs = false;
		if(!hetRateModified && normalPloidy==1) {
			setHeterozygosityRate(SingleSampleVariantPileupListener.DEF_HETEROZYGOSITY_RATE_HAPLOID);
		}
		logParameters();
		if(inputFile==null) throw new IOException("The input file with alignments is a required parameter");
		if (genome==null) throw new IOException("The reference genome file is a required parameter");
		
		referenceGenomeSize = genome.getTotalLength();
		log.info("Loaded "+genome.getNumSequences()+" sequences");
		if(progressNotifier!=null && !progressNotifier.keepRunning(1)) return;  
		calledSVs = new GenomicRegionSortedCollection<CalledGenomicVariant>(genome.getSequencesMetadata());
		if (knownSVsFile!=null) {
			calledSVs.addAll(svsFH.loadVariants(knownSVsFile));
			log.info("Loaded "+calledSVs.size()+" input SVs");
		}
		if(findRepeats) {
			log.info("Finding repeats using reads with multiple alignments");
			List<CalledCNV> multipleMCnvs = mmRegsCalc.calculateMultipleMappingRegions(inputFile);
			log.info("Found "+multipleMCnvs.size()+" repeats");
			calledSVs.addAll(multipleMCnvs);
			log.info("Number of SVs after finding repeats: "+calledSVs.size());
		}
		if(progressNotifier!=null && !progressNotifier.keepRunning(4)) return;
		//Call CNVs based on read depth
		if(runRDAnalysis) {
			log.info("Running read depth (RD) analysis to identify/genotype CNVs");
			List<CalledCNV> cnvsRD = runRDAnalysis();
			if(cnvsRD !=null) {
				log.info("Found "+cnvsRD.size()+" new CNVs running the RD analysis");
				calledSVs.addAll(cnvsRD);
			}
			log.info("Total number of SVs: "+calledSVs.size());
		}
		if(progressNotifier!=null && !progressNotifier.keepRunning(10)) return;
		if(runRPAnalysis) {
			log.info("Running read pair (RP) analysis to identify indels and inversions");
			List<CalledGenomicVariant> svsRP = runRPAnalysis(); 
			log.info("Found "+svsRP.size()+" new structural variants running the RP analysis");
			calledSVs.addAll(svsRP);
			log.info("Total number of SVs: "+calledSVs.size());
		}
		if(progressNotifier!=null && !progressNotifier.keepRunning(15)) return;
		if(findSNVs) {
			try {
				findSNVS();
			} finally {
				if(outVars!=null) outVars.close();
				dispose();
			}
		}
		if(runRDAnalysis || runRPAnalysis || findRepeats) {
			log.info("Saving structural variants");
			try (PrintStream outStructural = new PrintStream(outputPrefix+"_SV.gff")) {
				GFFVariantsFileHandler svHandler = new GFFVariantsFileHandler();
				svHandler.saveVariants(calledSVs.asList(), outStructural);
			}
		}
		log.info("Variants Detector Completed");
	}


	public void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input alignments file: "+inputFile);
		if (genome!=null) out.println("Loaded reference genome from: "+genome.getFilename());
		out.println("Prefix for output files: "+outputPrefix);
		out.println("Sample id: "+sampleId);
		out.println("Normal ploidy: "+normalPloidy);
		out.println("Print header with sample ploidy in the vcf file: "+printSamplePloidy);
		out.println("Minimum mapping quality to consider an alignment unique: "+getMinMQ());
		out.println("Find SNVs: "+findSNVs);
		if(findSNVs) {
			if(knownVariantsFile!=null) out.println("File with known variants to genotype: " + knownVariantsFile);
			if(generator.getQuerySeq()!=null) {
				out.println("Analyze only region at "+getQuerySeq()+":"+getQueryFirst()+"-"+getQueryLast());
			}
			out.println("Ignore variants in lower case reference positions: " + isIgnoreLowerCaseRef());
			out.println("Maximum number of alignments starting at the same position: " + getMaxAlnsPerStartPos());
			out.println("Process non unique primary alignments: " + isProcessNonUniquePrimaryAlignments());
			out.println("Process secondary alignments: " + isProcessSecondaryAlignments());
			out.println("Base pairs to ignore from the 5' end of each read: " + getBasesToIgnore5P());
			out.println("Base pairs to ignore from the 3' end of each read: " + getBasesToIgnore3P());
			out.println("Prior heterozygosity rate: "+getHeterozygosityRate());
			out.println("Maximum base quality score (PHRED): " + getMaxBaseQS());
			if(knownSTRsFile!=null) out.println("File with known short tandem repeats: " + knownSTRsFile);
			out.println("Minimum variant quality score (PHRED): " + getMinQuality());
			out.println("Call SNVs within STRs: " + isCallEmbeddedSNVs());
			out.println("Calculate a exact fisher test p-value for strand bias: "+isCalcStrandBias());
		}
		out.println("File with known structural variants: "+knownSVsFile);
		out.println("Min quality for structural variants (PHRED) : "+getMinSVQuality());
		out.println("Find repeats using reads with multiple alignments: "+findRepeats);
		out.println("Run RD analysis to genotype given SVs and find new CNVs: "+runRDAnalysis);
		if(runRDAnalysis) {
			out.println("Identify new CNVs using the RD data: "+findNewCNVs);
			out.println("Input genome size: "+getInputGenomeSize());
			out.println("Bin size: "+getBinSize());
			out.println("Algorithms for RD analysis: "+getAlgCNV());
			out.println("Max percentage of overlap between input CNVs and new CNVs: "+getMaxPCTOverlapCNVs());
		}
		out.println("Run RP analysis to find indels and inversions: "+runRPAnalysis);
		if(runRPAnalysis) {
			out.println("Max length of deletions found with RP analysis : "+getMaxLengthDeletion());
			out.println("Size of the seed for split-read alignments : "+getSplitReadSeed());
			out.println("Ignore proper pair flag for RP analysis : "+isIgnoreProperPairFlag());
		}
		log.info(os.toString());	
	}
	
	public List<CalledCNV> runRDAnalysis() throws IOException {
		log.info("Loading bins");
		ReadDepthDistribution rdDistribution = new ReadDepthDistribution(genome, binSize);
		log.info("Loaded bins. Assembly genome size: "+rdDistribution.getGenomeSize());
		//Pass parameters
		rdDistribution.setLog(this.getLog());
		rdDistribution.setMinMQ(generator.getMinMQ());
		
		
		log.info("Processing alignments file: "+inputFile);
		rdDistribution.processAlignments(inputFile);
		log.info("Processed alignments file: "+inputFile);
		if(progressNotifier!=null && !progressNotifier.keepRunning(7)) return new ArrayList<CalledCNV>();
		rdDistribution.correctDepthByGCContent();
		log.info("Corrected GCContent biases");
		log.info("Calculating read depth parameters");
		rdDistribution.calculateReadDepthDistParameters();
		log.info("Calculated read depth parameters. Mean read depth: "+rdDistribution.getMeanReadDepth()+". Standard deviation: "+rdDistribution.getSigmaReadDepth());
		GenomicRegionSortedCollection<CalledCNV> inputCNVs = new GenomicRegionSortedCollection<CalledCNV>();
		if(calledSVs!=null) {
			log.info("Calculating normalized read depth for input repeats and CNVs");
			inputCNVs = selectCalledCNVs(calledSVs);
			calculateNormalizedAverageDepth(rdDistribution,inputCNVs);
			updateDelDupStatus(inputCNVs.asList());
			log.info("Calculated normalized read depth for input CNVs and repeats");
		}
		if(findNewCNVs) {
			List<CalledCNV> cnvs = new ArrayList<CalledCNV>();
			
			String[] algs = algCNV.split(",");
			for (String algorithm : algs){
				SingleSampleReadDepthAlgorithm algor;
				try {
					Class<?> algorithmClass = (Class<?>) Class.forName("ngsep.discovery.rd."+algorithm+"ReadDepthAlgorithm");
					Constructor<?> constructor = algorithmClass.getDeclaredConstructors()[0];
					algor = (SingleSampleReadDepthAlgorithm) constructor.newInstance();
				} catch (Exception e) {
					throw new IOException("Unrecognized algorithm "+algorithm+" for read depth analysis", e);
				} 
				List<CalledCNV> cnvsA = executeCNValgorithm(algor, rdDistribution);
				cnvs.addAll(filterCNVs(cnvsA,inputCNVs));
			}
			updateDelDupStatus(cnvs);
			return cnvs;
		}
		return null;
	}
	
	

	private void updateDelDupStatus(List<CalledCNV> cnvs) {
		for(CalledCNV cnv:cnvs) {
			if(cnv.getNumCopies()<normalPloidy) cnv.setTextGenotype(CalledCNV.TEXT_GEN_DEL);
		}
		
	}


	private List<CalledCNV> executeCNValgorithm(SingleSampleReadDepthAlgorithm algorithm, ReadDepthDistribution rdDistribution){
		SingleSampleReadDepthAlgorithm rdAlgorithm = algorithm;
		rdAlgorithm.setLog(this.getLog());
		if(inputGenomeSize>0) rdAlgorithm.setGenomeSize(inputGenomeSize);
		else rdAlgorithm.setGenomeSize(rdDistribution.getGenomeSize());
		rdAlgorithm.setNormalPloidy((byte)normalPloidy);
		rdAlgorithm.setReadDepthDistribution(rdDistribution);
		return rdAlgorithm.callCNVs();
	}
	
	private List<CalledCNV> filterCNVs(List<CalledCNV> cnvs,GenomicRegionSortedCollection<CalledCNV> inputCNVs) {
		List<CalledCNV> answer = new ArrayList<CalledCNV>();
		
		for(CalledCNV cnv:cnvs) {
			if(cnv.getGenotypeQuality()<minSVQuality) continue;
			if(maxPCTOverlapCNVs<100) {
				GenomicRegionSortedCollection<CalledCNV> spanningCNVs = inputCNVs.findSpanningRegions(cnv);
				double maxSpan = 0;
				for(CalledCNV c2:spanningCNVs) {
					int nextSpan = GenomicRegionSpanComparator.getInstance().getSpanLength(cnv.getFirst(), cnv.getLast(), c2.getFirst(), c2.getLast());
					if(maxSpan<nextSpan) maxSpan= nextSpan;
				}
				double l = cnv.length();
				if(100.0*maxSpan/l > maxPCTOverlapCNVs) continue;
			}
			answer.add(cnv);
		}
		return answer;
	}
	private GenomicRegionSortedCollection<CalledCNV> selectCalledCNVs(GenomicRegionSortedCollection<CalledGenomicVariant> svs) {
		return selectCalledCNVs(svs,false);
	}
	private GenomicRegionSortedCollection<CalledCNV> selectCalledCNVs(GenomicRegionSortedCollection<CalledGenomicVariant> svs, boolean onlyRepeats) {
		GenomicRegionSortedCollection<CalledCNV> calledCNVs = new GenomicRegionSortedCollection<CalledCNV>(svs.getSequenceNames());
		for(CalledGenomicVariant sv:svs) {
			if(sv instanceof CalledCNV && (!onlyRepeats || sv.getType()==GenomicVariant.TYPE_REPEAT)) {
				calledCNVs.add((CalledCNV)sv);
			}
		}
		return calledCNVs;
	}

	public void calculateNormalizedAverageDepth(ReadDepthDistribution rdDistribution, GenomicRegionSortedCollection<CalledCNV> calledCNVs) {
		int binSize = rdDistribution.getBinSize();
		for(String seqName:calledCNVs.getSequenceNames().getNamesStringList()) {
			List<CalledCNV> seqCNVs = calledCNVs.getSequenceRegions(seqName).asList();
			
			List<ReadDepthBin> seqBins = rdDistribution.getBins(seqName);
			if(seqBins!=null) {
				for(CalledCNV cnv:seqCNVs) {
					int binStart = (cnv.getFirst()-1)/binSize;
					int binEnd = (cnv.getLast()-1)/binSize;
					double avg = 0;
					int sumUncorrected = 0;
					int nBins = 0;
					for(int i=binStart;i<seqBins.size()&& i<=binEnd;i++) {
						ReadDepthBin bin = seqBins.get(i);
						sumUncorrected += bin.getRawReadDepth();
						avg+=bin.getCorrectedReadDepth();
						nBins++;
					}
					//TODO: Update genotype quality
					if(nBins ==0) cnv.setNumCopies((byte)normalPloidy,true);
					else {
						cnv.setTotalReadDepth(sumUncorrected);
						avg/=nBins;
						cnv.setNumCopies((float) (avg*normalPloidy/rdDistribution.getMeanReadDepth()),true);
					}
				}
			} else {
				log.warning("Sequence name: "+seqName+" in input CNVs not found in the genome");
			}
		}
	}

	static GenomicRegionSortedCollection<GenomicVariant> makeNonRedundantSTRs( ReferenceGenome genome, List<GenomicRegion> strs) {
		QualifiedSequenceList sequences = genome.getSequencesMetadata();
		GenomicRegionSortedCollection<GenomicRegion> strsC = new GenomicRegionSortedCollection<GenomicRegion>(sequences);
		strsC.addAll(strs);
		GenomicRegionSortedCollection<GenomicVariant> answer = new GenomicRegionSortedCollection<GenomicVariant>(sequences);
		
		for(QualifiedSequence seq:sequences) {
			String seqName = seq.getName();
			int first = 0;
			int last = 0;
			GenomicRegionSortedCollection<GenomicRegion> seqSTRs = strsC.getSequenceRegions(seqName);
			for(GenomicRegion r:seqSTRs) {
				if(last == 0 || !mergeSTRs(genome, seqName,first,last,r)) {
					if(last>0) {
						GenomicVariant nextVar = makeSTRVariant(genome, seqName, Math.max(1, first-1),Math.min(last+1,seq.getLength()));
						if(nextVar!=null) answer.add(nextVar);
					}
					
					first = r.getFirst();
				}
				last = r.getLast();
			}
			if(last>0) {
				GenomicVariant nextVar = makeSTRVariant(genome, seqName, Math.max(1, first-1),Math.min(last+1,seq.getLength()));
				if(nextVar!=null) answer.add(nextVar);
			}
		}
		
		return answer;
	}
	private static boolean mergeSTRs(ReferenceGenome genome, String sequenceName, int first, int last, GenomicRegion r) {
		if(r.getFirst()-last>5) return false;
		else if (r.getFirst()-last<=2) return true;
		CharSequence ref1 = genome.getReference(sequenceName, Math.max(first, last-10), last);
		CharSequence ref2 = genome.getReference(sequenceName, r.getFirst(), r.getLast());
		if(ref1==null || ref2==null) return false;
		return AbstractLimitedSequence.getOverlapLength(ref1.toString().toUpperCase(), ref2.toString().toUpperCase())>5;
	}


	private static GenomicVariant makeSTRVariant(ReferenceGenome genome, String sequenceName, int first, int last) {
		List<String> alleles = new ArrayList<String>();
		CharSequence reference = genome.getReference(sequenceName, first, last);
		if(reference==null) {
			System.err.println("Reference not found for input STR at coordinates "+sequenceName+":"+first+"-"+last);
			return null;
		}
		alleles.add(reference.toString().toUpperCase());
		GenomicVariantImpl answer = new GenomicVariantImpl(sequenceName, first, alleles);
		answer.setType(GenomicVariant.TYPE_STR);
		return answer;
	}

	public void findSNVS() throws IOException {
		if(knownVariantsFile!=null) {
			log.info("Loading input variants");
			List<GenomicVariant> knownVariants = VCFFileReader.loadVariants(knownVariantsFile,true);
			log.info("Loaded "+knownVariants.size()+" input variants");
			GenomicRegionSortedCollection<GenomicVariant> knownVarsC = new GenomicRegionSortedCollection<GenomicVariant>(genome.getSequencesMetadata());
			knownVarsC.addAll(knownVariants);
			indelRealigner.setInputVariants(knownVarsC);
			varListener.setInputVariants(knownVarsC);
		} else if(knownSTRsFile!=null) {
			log.info("Loading input short tandem repeats");
			//TODO: Choose the best format
			SimpleGenomicRegionFileHandler rfh = new SimpleGenomicRegionFileHandler();
			List<GenomicRegion> strs = rfh.loadRegions(knownSTRsFile);
			indelRealigner.setInputVariants(makeNonRedundantSTRs(genome, strs));
			log.info("Loaded "+strs.size()+" input short tandem repeats");
		}
		log.info("Finding variants");
		header = VCFFileHeader.makeDefaultEmptyHeader();
		Sample s = new Sample(sampleId);
		s.setNormalPloidy(normalPloidy);
		header.addSample(s, printSamplePloidy);
		indelRealigner.setGenome(genome);
		generator.addListener(indelRealigner);
		varListener.clear();
		varListener.setGenome(genome);
		varListener.setSample(s);
		generator.addListener(varListener);
		generator.addListener(this);
		try (PrintStream outVars = new PrintStream(outputPrefix+".vcf")) {
			this.outVars = outVars;
			varsFW.printHeader(header,outVars);
			generator.processFile(inputFile);
		}	
	}

	private void saveSequenceVariants(String sequenceName) {
		List<CalledCNV> sequenceCNVs= selectCalledCNVs(calledSVs.getSequenceRegions(sequenceName)).asList();
		List<CalledGenomicVariant> sequenceVariants = varListener.getCalledVariants();
		boolean [] varInCNV = new boolean [sequenceVariants.size()]; 
		intersectVariantsCNVs(sequenceCNVs,sequenceVariants,varInCNV);
		for(int i=0;i<sequenceVariants.size();i++) {
			CalledGenomicVariant call = sequenceVariants.get(i);
			int [] format;
			if(call instanceof CalledSNV) format = VCFRecord.DEF_FORMAT_ARRAY_NGSEP_SNV;
			else if (call.getAllCounts()!=null) format = VCFRecord.DEF_FORMAT_ARRAY_NGSEP_SNV;
			else format = VCFRecord.DEF_FORMAT_ARRAY_NGSEP_NOSNV;
			VCFRecord record = new VCFRecord(call, format, call, header);
			if(call.getStrandBiasScore()!=CalledGenomicVariant.INVALID_STRAND_BIAS_SCORE) record.addAnnotation(new GenomicVariantAnnotation(call, GenomicVariantAnnotation.ATTRIBUTE_FISHER_STRAND_BIAS, call.getStrandBiasScore()));
			if(varInCNV[i]) record.addAnnotation(new GenomicVariantAnnotation(call, GenomicVariantAnnotation.ATTRIBUTE_IN_CNV, "1"));
			varsFW.printVCFRecord(record, outVars);
		}
		outVars.flush();
		varListener.clear();
	}
	private void intersectVariantsCNVs(List<CalledCNV> sequenceCNVs,List<CalledGenomicVariant> sequenceVars, boolean [] varInCNV) {
		int threshold = getBinSize();
		int indexFirst = 0;
		for(int i=0;i<sequenceVars.size();i++) {
			CalledGenomicVariant var = sequenceVars.get(i);
			//Clean up old estimate
			varInCNV[i] = false;
			var.updateAllelesCopyNumberFromCounts(normalPloidy);
			int first = Math.max(0, var.getFirst()-threshold);
			//Discard CNVs ending before start
			while(indexFirst<sequenceCNVs.size()) {
				CalledCNV cnv = sequenceCNVs.get(indexFirst);
				if(cnv.getLast() < first) {
					indexFirst++;
				} else {
					break;
				}
			}
			if(indexFirst==sequenceCNVs.size()) continue;
			for(int j=indexFirst;j<sequenceCNVs.size();j++) {
				CalledCNV cnv = sequenceCNVs.get(j);
				int startVicinity = cnv.getFirst()-threshold;
				int endVicinity = cnv.getLast()+threshold;
				boolean closeToCNV = startVicinity <= var.getLast() && var.getFirst() <= endVicinity;
				//boolean inCNV = cnv.getStart() <= var.getEnd() && var.getStart() <= cnv.getEnd();
				//double normalizedCoverage = (var.getCountReference()+var.getCountAlternative())/cnvListener.getHaploidAverageCoverage();
				//if(inCNV || (closeToCNV && Math.abs(normalizedCoverage-normalPloidy)>0.5)) {
				if(closeToCNV) {
					varInCNV[i] = true;
					byte numCopies = (byte)Math.round(cnv.getNumCopies());
					var.updateAllelesCopyNumberFromCounts(numCopies);
					if(var.isHeterozygous()) {
						cnv.addHeterozygousVariant();
					}
				} else if (cnv.getFirst() > var.getLast()+threshold) {
					break;
				}
			}
		}
		
	}

	@Override
	public void onPileup(PileupRecord pileup) {
		coveredGenomeSize++;
		if(progressNotifier!=null && coveredGenomeSize%10000==0) {
			int progress = 15+(int)Math.round(85.0*coveredGenomeSize/referenceGenomeSize);
			generator.setKeepRunning(progressNotifier.keepRunning(progress));
		}
	}

	
	@Override
	public void onSequenceStart(QualifiedSequence sequence) {
	}

	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {
		if(sequence==null) {
			log.warning("Null sequence");
			return;
		}
		saveSequenceVariants(sequence.getName());
	}
	
	private List<CalledGenomicVariant> runRPAnalysis() throws IOException {
		GenomicRegionSortedCollection<CalledCNV> duplications = new GenomicRegionSortedCollection<CalledCNV>(genome.getSequencesMetadata());
		GenomicRegionSortedCollection<CalledCNV> calledCNVs = selectCalledCNVs(calledSVs);
		for(CalledCNV cnv:calledCNVs) {
			if(cnv.getNumCopies()>1.25*normalPloidy) duplications.add(cnv);
		}
		log.info("Using "+duplications.size()+" duplications out of "+calledCNVs.size()+" svs in the read pair algorithm");
		rpAnalyzer.setReference(genome);
		rpAnalyzer.setDuplications(duplications);
		List<CalledGenomicVariant> svsRP = rpAnalyzer.findVariants(inputFile);
		log.info("Identified "+svsRP.size()+" candidate structural variants using the read pair algorithm. Filtering by quality score");
		svsRP = filterSVsReadPair(svsRP);
		for(CalledCNV cnv:duplications) {
			if(cnv.getTandemFragments()>1 && cnv.getTandemFragments()>3*cnv.getTransDupFragments()) cnv.setTextGenotype(CalledCNV.TEXT_GEN_TANDEMDUP);
			else if(cnv.getTransDupFragments()>1 && cnv.getTransDupFragments()>3*cnv.getTandemFragments()) cnv.setTextGenotype(CalledCNV.TEXT_GEN_TRANSDUP);
		}
		return svsRP;
	}
	
	private List<CalledGenomicVariant> filterSVsReadPair(List<CalledGenomicVariant> svs) {
		List<CalledGenomicVariant> answer = new ArrayList<CalledGenomicVariant>();
		for(CalledGenomicVariant v:svs) {
			if(v.getGenotypeQuality()>=minSVQuality) answer.add(v);
		}
		return answer;
	}

	public GenomicRegionSortedCollection<CalledGenomicVariant> getCalledSVs() {
		return calledSVs;
	}

	/**
	 * Removes heavy resources loaded in memory during the process
	 */
	private void dispose() {
		indelRealigner.setInputVariants(null);
		varListener.setInputVariants(null);
		varListener.clear();
	}
}
