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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class MultisampleVariantsDetector implements PileupListener {

	// Constants for default values
	public static final int DEF_MAX_ALNS_PER_START_POS = AlignmentsPileupGenerator.DEF_MAX_ALNS_PER_START_POS;
	public static final double DEF_MIN_ALLELE_DEPTH_FREQUENCY = 0;
	public static final double DEF_HETEROZYGOSITY_RATE_DIPLOID = CountsHelper.DEF_HETEROZYGOSITY_RATE_DIPLOID;
	public static final short DEF_MIN_QUALITY = 40;
	public static final short DEF_MIN_MQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	public static final byte DEF_MAX_BASE_QS = CountsHelper.DEF_MAX_BASE_QS;
	public static final byte DEF_PLOIDY = GenomicVariant.DEFAULT_PLOIDY;
	public static final String DEF_OUTPUT_FILE = "variants.vcf";
	
	// Logging and progress
	private Logger log = Logger.getLogger(MultisampleVariantsDetector.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	
	private List<String> inputFiles = new ArrayList<>();
	private ReferenceGenome genome = null;
	private String outFilename = DEF_OUTPUT_FILE;
	private PrintStream outFile;
	private VCFFileHeader vcfFileHeader;
	private VCFFileWriter writer = new VCFFileWriter();
	
	private double heterozygosityRate = DEF_HETEROZYGOSITY_RATE_DIPLOID;
	private AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
	 
	private boolean ignoreLowerCaseRef = false;
	private boolean callEmbeddedSNVs = false;
	private double minAlleleDepthFrequency = DEF_MIN_ALLELE_DEPTH_FREQUENCY;
	private short minQuality = DEF_MIN_QUALITY;
	private byte maxBaseQS = DEF_MAX_BASE_QS;
	private short normalPloidy = DEF_PLOIDY;
	private boolean printSamplePloidy = false;
	
	private String knownSTRsFile = null;
	private String knownVariantsFile=null;
	private GenomicRegionSortedCollection<GenomicVariant> inputVariants = new GenomicRegionSortedCollection<>();
	
	// Model attributes
	private IndelRealignerPileupListener indelRealigner = new IndelRealignerPileupListener();
	
	private List<Sample> samples;
	private double coveredGenomeSize = 0;
	private long referenceGenomeSize = 0;
	
	//DEBUG
	private int posPrint = -1;
	
	
	
	// Get and set methods
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
		generator.setLog(log);
	}

	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
	
	public List<String> getInputFiles() {
		return inputFiles;
	}
	public void setInputFiles(List<String> inputFiles) {
		this.inputFiles = inputFiles;
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

	/**
	 * @return the outFilename
	 */
	public String getOutFilename() {
		return outFilename;
	}

	/**
	 * @param outFilename the outFilename to set
	 */
	public void setOutFilename(String outFilename) {
		this.outFilename = outFilename;
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
		generator.setMinMQ(minMQ);
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
		return ignoreLowerCaseRef;
	}
	public void setIgnoreLowerCaseRef(boolean ignoreLowerCaseRef) {
		this.ignoreLowerCaseRef = ignoreLowerCaseRef;
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
	public void setBasesToIgnore5P(String basesToIgnore5P) {
		setBasesToIgnore5P((byte)OptionValuesDecoder.decode(basesToIgnore5P, Byte.class));
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
	public void setBasesToIgnore3P(String basesToIgnore3P) {
		setBasesToIgnore3P((byte)OptionValuesDecoder.decode(basesToIgnore3P, Byte.class));
	}

	public double getHeterozygosityRate() {
		return heterozygosityRate;
	}
	public void setHeterozygosityRate(double heterozygosityRate) {
		this.heterozygosityRate = heterozygosityRate;
	}	
	public void setHeterozygosityRate(String value) {
		setHeterozygosityRate((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	public byte getMaxBaseQS() {
		return maxBaseQS;
	}
	public void setMaxBaseQS(byte maxBaseQS) {
		this.maxBaseQS = maxBaseQS;
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
		return minQuality;
	}
	public void setMinQuality(short minQuality) {
		this.minQuality = minQuality;
	}	
	public void setMinQuality(String value) {
		setMinQuality((short)OptionValuesDecoder.decode(value, Short.class));
	}

	public boolean isCallEmbeddedSNVs() {
		return callEmbeddedSNVs;
	}
	public void setCallEmbeddedSNVs(boolean callEmbeddedSNVs) {
		this.callEmbeddedSNVs = callEmbeddedSNVs;
	}	
	public void setCallEmbeddedSNVs(Boolean callEmbeddedSNVs) {
		setCallEmbeddedSNVs(callEmbeddedSNVs.booleanValue());
	}
	
	
	
	public double getMinAlleleDepthFrequency() {
		return minAlleleDepthFrequency;
	}
	public void setMinAlleleDepthFrequency(double minAlleleDepthFrequency) {
		this.minAlleleDepthFrequency = minAlleleDepthFrequency;
	}
	public void setMinAlleleDepthFrequency(String value) {
		setMinAlleleDepthFrequency((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	public List<Sample> getSamples() {
		return samples;
	}
	public void setSamples(List<Sample> samples) {
		this.samples = samples;
	}
	
	public static void main(String[] args) throws Exception {
		MultisampleVariantsDetector instance = new MultisampleVariantsDetector();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		for(;i<args.length;i++) {
			instance.inputFiles.add(args[i]);
		}
		instance.run();
	}
	
	public void run() throws IOException {
		logParameters();
		if (genome==null) throw new IOException("The reference genome file is a required parameter");
		referenceGenomeSize = genome.getTotalLength();
		QualifiedSequenceList sequences = genome.getSequencesMetadata();
		indelRealigner.setGenome(genome);
		generator.setSequencesMetadata(sequences);
		//TODO: assign sample ids if not in aln files
		if(samples == null) loadSamplesFromAlignmentHeaders();
		
		if(knownVariantsFile!=null) {
			log.info("Loading input variants");
			List<GenomicVariant> knownVariants = VCFFileReader.loadVariants(knownVariantsFile,true);
			log.info("Loaded "+knownVariants.size()+" input variants");
			inputVariants = new GenomicRegionSortedCollection<GenomicVariant>(sequences);
			inputVariants.addAll(knownVariants);
			indelRealigner.setInputVariants(inputVariants);
		} else if(knownSTRsFile!=null) {
			log.info("Loading input short tandem repeats from: "+knownSTRsFile);
			//TODO: STRs loader
			SimpleGenomicRegionFileHandler rfh = new SimpleGenomicRegionFileHandler();
			List<GenomicRegion> strs = rfh.loadRegions(knownSTRsFile);
			indelRealigner.setInputVariants(SingleSampleVariantsDetector.makeNonRedundantSTRs(genome,strs));
			log.info("Loaded "+strs.size()+" input short tandem repeats");
		}
		log.info("Finding variants");
		
		generator.addListener(indelRealigner);
		generator.addListener(this);
		try {
			outFile = new PrintStream(outFilename);
			vcfFileHeader = VCFFileHeader.makeDefaultEmptyHeader();
			for(Sample s:samples) vcfFileHeader.addSample(s, printSamplePloidy);
			writer.printHeader(vcfFileHeader, outFile);
			generator.processFiles(inputFiles);
		} finally {
			if(outFile!=null) outFile.close();
			dispose();
		}
		log.info("Multisample Variants Detector Completed");
	}
	
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		log.info("Input files: "+inputFiles);
		if (genome!=null) out.println("Loaded reference genome from: "+genome.getFilename());
		out.println("Output file: "+outFilename);
		if(knownVariantsFile!=null) out.println("File with known variants to genotype: " + knownVariantsFile);
		if(generator.getQuerySeq()!=null) {
			out.println("Analyze only region at "+getQuerySeq()+":"+getQueryFirst()+"-"+getQueryLast());
		}
		out.println("Ignore variants in lower case reference positions: " + isIgnoreLowerCaseRef());
		out.println("Maximum number of alignments starting at the same position: " + getMaxAlnsPerStartPos());
		out.println("Minimum mapping quality to consider an alignment unique: "+getMinMQ());
		out.println("Process non unique primary alignments: " + isProcessNonUniquePrimaryAlignments());
		out.println("Process secondary alignments: " + isProcessSecondaryAlignments());
		out.println("Base pairs to ignore from the 5' end of each read: " + getBasesToIgnore5P());
		out.println("Base pairs to ignore from the 3' end of each read: " + getBasesToIgnore3P());
		out.println("Prior heterozygosity rate: "+getHeterozygosityRate());
		out.println("Maximum base quality score (PHRED): " + getMaxBaseQS());
		if(knownSTRsFile!=null) out.println("File with known short tandem repeats: " + knownSTRsFile);
		out.println("Minimum variant quality score (PHRED): " + getMinQuality());
		out.println("Call SNVs within STRs: " + isCallEmbeddedSNVs());
		out.println("Normal ploidy: "+normalPloidy);
		out.println("Print header with sample ploidy in the vcf file: "+printSamplePloidy);
		log.info(os.toString());
	}

	
	private void loadSamplesFromAlignmentHeaders() throws IOException {
		Map<String, Sample> samplesMap = new TreeMap<>();
		for(String filename:inputFiles) {
			try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(filename)) {
				Map<String, String> samplesHeader = reader.getSampleIdsByReadGroup();
				for(String readGroup : samplesHeader.keySet()) {
					String sampleId = samplesHeader.get(readGroup);
					Sample sample = samplesMap.get(sampleId);
					if(sample==null) {
						sample = new Sample(sampleId);
						sample.setNormalPloidy(normalPloidy);
						samplesMap.put(sampleId, sample);
						log.info("Found sample: "+sampleId+" in file: "+filename);
					}
					sample.addReadGroup(readGroup);
				}
			}
		}
		samples = new ArrayList<>(samplesMap.values()); 
	}
	//Control attribute to avoid calling overlapping indels and to give an embedded status to SNVs within indels or STRs
	private int lastIndelEnd = 0;
	private int nextSIVIndex = 0;
	private List<GenomicVariant> seqInputVariants;
	@Override
	public void onPileup(PileupRecord pileup) {
		GenomicVariant variant = null;
		GenomicVariant inputVariant = null;
		if(inputVariants.size()==0) {
			if(pileup.isInputSTR()) lastIndelEnd = pileup.getPosition()+pileup.getReferenceSpan()-1;
			else if(pileup.getPosition()<=lastIndelEnd) pileup.setEmbedded(true);
			String referenceAllele = SingleSampleVariantPileupListener.calculateReferenceAlleleDiscovery(pileup,genome,callEmbeddedSNVs,ignoreLowerCaseRef);
			if(referenceAllele==null) return;
			
			variant = discoverPopulationVariant(pileup,referenceAllele);
			if(pileup.getPosition()==posPrint) System.out.println("Variant: "+variant);
			if(variant == null) return;
			List<CalledGenomicVariant> calls = genotypeVariant(variant, pileup);
			if (variant.getVariantQS()==0 || variant.getVariantQS() < minQuality) return;
			VCFRecord record = VCFRecord.createDefaultPopulationVCFRecord(variant, calls, vcfFileHeader);
			writer.printVCFRecord(record, outFile);
			if(!variant.isSNV()) {
				lastIndelEnd = variant.getLast();
			}
		} else if(nextSIVIndex<seqInputVariants.size()) {
			inputVariant = seqInputVariants.get(nextSIVIndex);
			while(inputVariant.getFirst() <= pileup.getPosition() ) {
				if(inputVariant.getFirst()==pileup.getPosition()) {
					List<CalledGenomicVariant> calls = genotypeVariant(variant, pileup);
					VCFRecord record = VCFRecord.createDefaultPopulationVCFRecord(variant, calls, vcfFileHeader);
					writer.printVCFRecord(record, outFile);
				}
				nextSIVIndex++;
				if(nextSIVIndex>=seqInputVariants.size()) break;
				inputVariant = seqInputVariants.get(nextSIVIndex);
			}
		}
		
		coveredGenomeSize++;
		if(progressNotifier!=null && coveredGenomeSize%10000==0) {
			int progress = (int)Math.round(100.0*coveredGenomeSize/referenceGenomeSize);
			generator.setKeepRunning(progressNotifier.keepRunning(progress));
		}
	}
	
	
	@Override
	public void onSequenceStart(QualifiedSequence sequence) {
		if(inputVariants.size()>0) seqInputVariants = inputVariants.getSequenceRegions(sequence.getName()).asList();
		nextSIVIndex = 0;
		lastIndelEnd = 0;
	}

	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {
		
	}
	
	public GenomicVariant discoverPopulationVariant(PileupRecord pileup, String referenceAllele) {
		GenomicVariant variant;
		if(referenceAllele.length()>1) {
			variant = discoverPopulationVariantWithSpan(pileup, referenceAllele);
		} else {
			variant = discoverPopulationSNV(pileup, referenceAllele.charAt(0));
		}
		if(variant != null) {
			if(variant.isSNV() && pileup.isEmbedded()) variant.setType(GenomicVariant.TYPE_EMBEDDED_SNV);
		}
		return variant;
	}
	public GenomicVariant discoverPopulationSNV(PileupRecord pileup, char reference) {
		List<PileupAlleleCall> alleleCalls = pileup.getAlleleCalls(1,(String)null);
		CountsHelper helperSNV = CountsHelper.calculateCountsSNV(alleleCalls, maxBaseQS);
		GenomicVariant variant = SingleSampleVariantPileupListener.createSNVVariantPool(pileup, helperSNV, reference, minAlleleDepthFrequency);
		if(variant == null) return null;
		while(variant.getAlleles().length > 2) {	
			List<CalledGenomicVariant> calls = genotypeVariant(variant, pileup, heterozygosityRate);
			GenomicVariant newVariant = makeNewVariant(variant, calls);
			if(newVariant!=variant) variant = newVariant;
			else break;
		}
		return variant;
	}
	
	private GenomicVariant discoverPopulationVariantWithSpan(PileupRecord pileup, String referenceAllele) {
		List<PileupAlleleCall> calls = pileup.getAlleleCalls(referenceAllele.length(),(String)null);
		AlleleCallClustersBuilder acBuilder = new AlleleCallClustersBuilder(pileup.getSequenceName(),pileup.getPosition());
		String [] alleles =  acBuilder.clusterAlleleCalls(pileup, calls, referenceAllele, maxBaseQS);
		CountsHelper helper = CountsHelper.calculateCounts(alleles, calls,maxBaseQS); 
		GenomicVariant variant = discoverPopulationIndel(pileup, helper);
		if(!pileup.isInputSTR() && variant==null) {
			if (pileup.isNewSTR()) {
				pileup.setSTR(false);
				pileup.setNewSTR(false);
			}
			//Try SNV if the indel alleles were not good to make a call
			variant = discoverPopulationSNV(pileup, referenceAllele.charAt(0));
		}
		return variant;
	}

	private GenomicVariant discoverPopulationIndel(PileupRecord pileup, CountsHelper helperIndel) {
		GenomicVariant variant = SingleSampleVariantPileupListener.createIndelVariantPool(pileup, helperIndel);
		if(variant==null) return null;
		// Redefine indel alleles based on genotype calls
		while(variant.getAlleles().length>2) {
			if(!pileup.isInputSTR() && SingleSampleVariantPileupListener.allelesSameLength(variant.getAlleles())) return null;
			List<CalledGenomicVariant> calls = genotypeVariant(variant, pileup, heterozygosityRate);
			if(variant.getVariantQS() < minQuality) return null;
			GenomicVariant newVariant = makeNewVariant(variant, calls);
			if(newVariant!=variant) variant = newVariant;
			else break;
		}
		if (pileup.isSTR()) {
			variant.setType(GenomicVariant.TYPE_STR);
		} else {
			variant.setType(GenomicVariant.TYPE_INDEL);
		}
		return variant;
	}
	
	/**
	 * 
	 * @param variant to evaluate
	 * @param calls Genotype calls for the given variant
	 * @return GenomicVariant with less alleles or the same variant if it can not be changed
	 */
	private GenomicVariant makeNewVariant (GenomicVariant variant, List<CalledGenomicVariant> calls) {
		if(variant.getAlleles().length<=2) return variant;
		Set<String> calledAllelesSet = new TreeSet<>();
		calledAllelesSet.add(variant.getReference());
		int n = samples.size();
		for(int i=0;i<n;i++) {
			CalledGenomicVariant call = calls.get(i);
			calledAllelesSet.addAll(Arrays.asList(call.getCalledAlleles()));
		}
		if(variant.getAlleles().length !=calledAllelesSet.size()) {
			variant = SingleSampleVariantPileupListener.makeNewVariant(variant,calledAllelesSet);
		}
		return variant;
		
	}

	/**
	 * Calls genotypes and updates the variant quality
	 * @param variant
	 * @param pileup
	 * @return List<CalledGenomicVariant> variants called for the samples included in this variants detector
	 */
	public List<CalledGenomicVariant> genotypeVariant(GenomicVariant variant, PileupRecord pileup) {
		return genotypeVariant(variant, pileup, heterozygosityRate);
	}
	/**
	 * Calls genotypes and updates the variant quality
	 * @param variant
	 * @param pileup
	 * @param h
	 * @return
	 */
	private List<CalledGenomicVariant> genotypeVariant(GenomicVariant variant, PileupRecord pileup, double h) {
		if(pileup.getPosition()==posPrint) System.out.println("Genotyping variant type: "+variant.getType()+" is SNV: "+variant.isSNV()+" alleles: "+Arrays.asList(variant.getAlleles()));
		List<CalledGenomicVariant> calls = new ArrayList<>();
		SingleSampleVariantPileupListener sampleDetector = new SingleSampleVariantPileupListener();
		sampleDetector.setMaxBaseQS(maxBaseQS);
		int n = samples.size();
		short variantQS = 0;
		for(int i=0;i<n;i++) {
			Sample sample = samples.get(i);
			CalledGenomicVariant call = sampleDetector.genotypeVariantSample(variant, pileup, sample, h);
			if(pileup.getPosition()==posPrint) System.out.println("Sample: "+call.getSampleId()+" Genotype: "+Arrays.asList(call.getCalledAlleles())+" GQ: "+call.getGenotypeQuality());
			if(!call.isUndecided() && !call.isHomozygousReference() && call.getGenotypeQuality()>variantQS) {
				variantQS = call.getGenotypeQuality();
			}
			calls.add(call);
		}
		//TODO: The variant could be genotyped again with a different heterozygosity rate
		variant.setVariantQS(variantQS);
		return calls;
	}
	
	private void dispose() {
		inputVariants =null;
		seqInputVariants = null;
		indelRealigner.setInputVariants(null);
	}


}
