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
import ngsep.math.NumberArrays;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.DiversityStatistics;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class MultisampleVariantsDetector implements PileupListener {

	private Logger log = Logger.getLogger(MultisampleVariantsDetector.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	public static final int DEF_MAX_ALNS_PER_START_POS = AlignmentsPileupGenerator.DEF_MAX_ALNS_PER_START_POS;
	public static final double DEF_MIN_ALLELE_FREQUENCY = 0;
	public static final double DEF_MIN_HETEROZYGOSITY_RATE_DIPLOID = VariantPileupListener.DEF_HETEROZYGOSITY_RATE_DIPLOID;
	public static final short DEF_MIN_QUALITY = 40;
	public static final short DEF_MIN_MQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	public static final byte DEF_MAX_BASE_QS = VariantPileupListener.DEF_MAX_BASE_QS;
	public static final byte DEF_PLOIDY = GenomicVariant.DEFAULT_PLOIDY;
	public static final String DEF_OUTPUT_FILE = "variants.vcf";
	
	
	
	private double coveredGenomeSize = 0;
	private long referenceGenomeSize = 0;
	
	private AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
	//Listeners
	private IndelRealignerPileupListener indelRealigner = new IndelRealignerPileupListener();

	
	private ReferenceGenome genome;
	private List<String> alignmentFiles = new ArrayList<>();
	private List<Sample> samples;
	
	//Output file variables
	private String outFilename = DEF_OUTPUT_FILE;
	private PrintStream outFile;
	private VCFFileHeader vcfFileHeader;
	private VCFFileWriter writer = new VCFFileWriter();
	
	private double heterozygosityRate = DEF_MIN_HETEROZYGOSITY_RATE_DIPLOID;
	private String knownSTRsFile = null; 
	private boolean ignoreLowerCaseRef = false;
	private boolean callEmbeddedSNVs = false;
	private double minAlleleFrequency = DEF_MIN_ALLELE_FREQUENCY;
	private short minQuality = DEF_MIN_QUALITY;
	private byte maxBaseQS = DEF_MAX_BASE_QS;
	private byte normalPloidy = DEF_PLOIDY;
	private boolean printSamplePloidy = false;
	private String knownVariantsFile=null;
	
	private GenomicRegionSortedCollection<GenomicVariant> inputVariants = new GenomicRegionSortedCollection<>();
	
	
	//Control attribute to avoid calling overlapping indels and to give an embedded status to SNVs within indels or STRs
	private int lastIndelEnd = 0;
	
	//DEBUG
	private int posPrint = -1;
	
	public static void main(String[] args) throws Exception {
		MultisampleVariantsDetector instance = new MultisampleVariantsDetector();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		for(;i<args.length;i++) {
			instance.alignmentFiles.add(args[i]);
		}
		instance.findVariants();
	}
	
	/**
	 * @return the log
	 */
	public Logger getLog() {
		return log;
	}

	/**
	 * @param log the log to set
	 */
	public void setLog(Logger log) {
		this.log = log;
		generator.setLog(log);
	}

	/**
	 * @return the progressNotifier
	 */
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}

	/**
	 * @param progressNotifier the progressNotifier to set
	 */
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}
	
	/**
	 * @return the alignmentFiles
	 */
	public List<String> getAlignmentFiles() {
		return alignmentFiles;
	}

	/**
	 * @param alignmentFiles the alignmentFiles to set
	 */
	public void setAlignmentFiles(List<String> alignmentFiles) {
		this.alignmentFiles = alignmentFiles;
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

	/**
	 * @return the heterozygosityRate
	 */
	public double getHeterozygosityRate() {
		return heterozygosityRate;
	}
	
	/**
	 * @param heterozygosityRate the heterozygosityRate to set
	 */
	public void setHeterozygosityRate(double heterozygosityRate) {
		this.heterozygosityRate = heterozygosityRate;
	}
	
	public void setHeterozygosityRate(String value) {
		setHeterozygosityRate((double)OptionValuesDecoder.decode(value, Double.class));
	}


	/**
	 * @return the minQuality
	 */
	public short getMinQuality() {
		return minQuality;
	}


	/**
	 * @param minQuality the minQuality to set
	 */
	public void setMinQuality(short minQuality) {
		this.minQuality = minQuality;
	}
	
	public void setMinQuality(String value) {
		setMinQuality((short)OptionValuesDecoder.decode(value, Short.class));
	}


	/**
	 * @return the maxBaseQS
	 */
	public byte getMaxBaseQS() {
		return maxBaseQS;
	}

	/**
	 * @param maxBaseQS the maxBaseQS to set
	 */
	public void setMaxBaseQS(byte maxBaseQS) {
		this.maxBaseQS = maxBaseQS;
	}
	
	public void setMaxBaseQS(String value) {
		setMaxBaseQS((byte)OptionValuesDecoder.decode(value, Byte.class));
	}

	/**
	 * @return the ignoreLowerCaseRef
	 */
	public boolean isIgnoreLowerCaseRef() {
		return ignoreLowerCaseRef;
	}

	/**
	 * @param ignoreLowerCaseRef the ignoreLowerCaseRef to set
	 */
	public void setIgnoreLowerCaseRef(boolean ignoreLowerCaseRef) {
		this.ignoreLowerCaseRef = ignoreLowerCaseRef;
	}
	
	public void setIgnoreLowerCaseRef(Boolean ignoreLowerCaseRef) {
		setIgnoreLowerCaseRef(ignoreLowerCaseRef.booleanValue());
	}

	/**
	 * @return the callEmbeddedSNVs
	 */
	public boolean isCallEmbeddedSNVs() {
		return callEmbeddedSNVs;
	}

	/**
	 * @param callEmbeddedSNVs the callEmbeddedSNVs to set
	 */
	public void setCallEmbeddedSNVs(boolean callEmbeddedSNVs) {
		this.callEmbeddedSNVs = callEmbeddedSNVs;
	}
	
	public void setCallEmbeddedSNVs(Boolean callEmbeddedSNVs) {
		setCallEmbeddedSNVs(callEmbeddedSNVs.booleanValue());
	}
	
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;	
	}
	
	/**
	 * @return the genome
	 */
	public ReferenceGenome getGenome() {
		return genome;
	}

	/**
	 * @return the minAlleleFrequency
	 */
	public double getMinAlleleFrequency() {
		return minAlleleFrequency;
	}

	/**
	 * @param minAlleleFrequency the minAlleleFrequency to set
	 */
	public void setMinAlleleFrequency(double minAlleleFrequency) {
		this.minAlleleFrequency = minAlleleFrequency;
	}
	
	public void setMinAlleleFrequency(String value) {
		setMinAlleleFrequency((double)OptionValuesDecoder.decode(value, Double.class));
	}
	
	/**
	 * @return the normalPloidy
	 */
	public byte getNormalPloidy() {
		return normalPloidy;
	}

	/**
	 * @param normalPloidy the normalPloidy to set
	 */
	public void setNormalPloidy(byte normalPloidy) {
		this.normalPloidy = normalPloidy;
	}
	
	public void setNormalPloidy(String value) {
		setNormalPloidy((byte)OptionValuesDecoder.decode(value, Byte.class));
	}
	
	/**
	 * @return the printSamplePloidy
	 */
	public boolean isPrintSamplePloidy() {
		return printSamplePloidy;
	}

	/**
	 * @param printSamplePloidy the printSamplePloidy to set
	 */
	public void setPrintSamplePloidy(boolean printSamplePloidy) {
		this.printSamplePloidy = printSamplePloidy;
	}

	public void setPrintSamplePloidy(Boolean printSamplePloidy) {
		this.setPrintSamplePloidy(printSamplePloidy.booleanValue());
	}

	/**
	 * @return
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
	 * @return
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
	 * @return
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
	
	/**
	 * @return
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
	 * @return
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
	 * @return
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

	/**
	 * @return
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
	
	

	/**
	 * @return the knownSTRsFile
	 */
	public String getKnownSTRsFile() {
		return knownSTRsFile;
	}

	/**
	 * @param knownSTRsFile the knownSTRsFile to set
	 */
	public void setKnownSTRsFile(String knownSTRsFile) {
		this.knownSTRsFile = knownSTRsFile;
	}
	
	public String getKnownVariantsFile() {
		return knownVariantsFile;
	}

	public void setKnownVariantsFile(String knownVariantsFile) {
		this.knownVariantsFile = knownVariantsFile;
	}
	
	public void printParameters() {
		log.info("Output file: "+outFilename);
		log.info("Heterozygosity rate: "+getHeterozygosityRate());
		if(generator.getQuerySeq()!=null) {
			log.info("Analyze only region at "+generator.getQuerySeq()+":"+generator.getQueryFirst()+"-"+generator.getQueryLast());
		}
		log.info("Maximum number of alignments starting at the same position: " + generator.getMaxAlnsPerStartPos());
		log.info("Ignore variants in lower case reference positions: " + isIgnoreLowerCaseRef());
		log.info("Process non unique primary alignments: " + generator.isProcessNonUniquePrimaryAlignments());
		log.info("Process secondary alignments: " + generator.isProcessSecondaryAlignments());
		log.info("Minimum variant quality score (PHRED): " + getMinQuality());
		log.info("Maximum base quality score (PHRED): " + getMaxBaseQS());
		log.info("Base pairs to ignore from the 5' end of each read: " + generator.getBasesToIgnore5P());
		log.info("Base pairs to ignore from the 3' end of each read: " + generator.getBasesToIgnore3P());
		log.info("File with known short tandem repeats: " + knownSTRsFile);
		log.info("File with known variants: " + knownVariantsFile);
		log.info("Call SNVs within STRs: " + isCallEmbeddedSNVs());
		
		log.info("Minimum mapping quality to consider an alignment unique: "+getMinMQ());
		log.info("Normal ploidy: "+normalPloidy);
		log.info("Print header with sample ploidy in the vcf file: "+printSamplePloidy);
	}

	public void findVariants() throws IOException {
		printParameters();
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
			indelRealigner.setInputVariants(VariantsDetector.makeNonRedundantSTRs(genome,strs));
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
			generator.processFiles(alignmentFiles);
		} finally {
			if(outFile!=null) outFile.close();
			dispose();
		}
		log.info("Multisample Variants Detector Completed");
	}

	
	private void loadSamplesFromAlignmentHeaders() throws IOException {
		Map<String, Sample> samplesMap = new TreeMap<>();
		log.info("Loading sample ids from: "+alignmentFiles);
		for(String filename:alignmentFiles) {
			
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

	private int nextSIVIndex = 0;
	private List<GenomicVariant> seqInputVariants;
	@Override
	public void onPileup(PileupRecord pileup) {
		GenomicVariant variant = null;
		GenomicVariant inputVariant = null;
		if(inputVariants.size()==0) {
			variant = findMultiallelicVariant(pileup);
		} else if(nextSIVIndex<seqInputVariants.size()) {
			inputVariant = seqInputVariants.get(nextSIVIndex);
			while(inputVariant.getFirst() <= pileup.getPosition() ) {
				if(inputVariant.getFirst()==pileup.getPosition()) {
					variant = inputVariant;
				}
				nextSIVIndex++;
				if(nextSIVIndex>=seqInputVariants.size()) break;
				inputVariant = seqInputVariants.get(nextSIVIndex);
			}
		}
		if(pileup.getPosition()==posPrint) System.out.println("Variant: "+variant);
		if(variant == null) return;
		
		List<CalledGenomicVariant> calls = genotypeVariant(variant, pileup, heterozygosityRate);
		if(inputVariant==null && (variant.getVariantQS()==0 || variant.getVariantQS() < minQuality)) return;
		//TODO: The variant could be genotyped again with a different heterozygosity rate
		
		DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(calls, false);
		int [] format = variant.isSNV()?VCFRecord.DEF_FORMAT_ARRAY_NGSEP_SNV:VCFRecord.DEF_FORMAT_ARRAY_NGSEP_NOSNV;
		VCFRecord record = new VCFRecord(variant, format, calls, vcfFileHeader);
		record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_SAMPLES_GENOTYPED, divStats.getNumSamplesGenotyped()));
		record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_NUMBER_ALLELES, divStats.getNumCalledAlleles()));
		record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_ALLELE_FREQUENCY_SPECTRUM, format(divStats.getAlleleCounts())));
		if(variant.isBiallelic()) record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_MAF, divStats.getMaf()));
		
		writer.printVCFRecord(record, outFile);
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
	

	private String format(int[] alleleCounts) {
		StringBuilder answer = new StringBuilder(""+alleleCounts[0]);
		for(int i=1;i<alleleCounts.length;i++) answer.append(","+alleleCounts[i]);
		return answer.toString();
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
			variant = makeNewVariant(variant,calledAllelesSet);
		}
		return variant;
		
	}

	private GenomicVariant makeNewVariant(GenomicVariant variant, Set<String> newAlleles) {
		if(variant.getFirst()==posPrint) System.out.println("Recoding alleles for "+variant.getFirst()+" alleles: "+Arrays.asList(variant.getAlleles()));
		List<String> alleles = new ArrayList<>(newAlleles.size());
		String reference = variant.getReference(); 
		alleles.add(reference);
		for(String allele:newAlleles) {
			if(!allele.equals(reference)) alleles.add(allele);
		}
		if(variant.getFirst()==posPrint) System.out.println("New alleles: "+alleles);
		if(variant.isSNV() && alleles.size()==2) {
			return new SNV(variant.getSequenceName(), variant.getFirst(), reference.charAt(0), alleles.get(1).charAt(0));
		}
		return new GenomicVariantImpl(variant.getSequenceName(), variant.getFirst(), alleles);
	}

	private GenomicVariant findMultiallelicVariant(PileupRecord pileup) {
		
		if(!callEmbeddedSNVs && pileup.isEmbedded()) return null;
		
		//Infer reference allele
		int last = pileup.getPosition()+pileup.getReferenceSpan()-1;
		CharSequence seq = genome.getReference(pileup.getSequenceName(), pileup.getPosition(), last);
		if(pileup.getPosition()==posPrint) System.out.println("Position: "+pileup.getPosition()+" Last: "+last);
		if(seq == null) return null;
		String referenceAllele = seq.toString();
		if(pileup.getPosition()==posPrint) System.out.println("Reference: "+referenceAllele);
		if(ignoreLowerCaseRef && Character.isLowerCase(referenceAllele.charAt(0))) return null;
		referenceAllele = referenceAllele.toUpperCase();
		//Avoid trying to call nested indels or SNVs within indels unless explicitly requested
		if(lastIndelEnd>=pileup.getPosition()) {
			if(!callEmbeddedSNVs) return null;
			referenceAllele = referenceAllele.substring(0,1);
			if (pileup.isSTR()) {
				pileup.setSTR(false);
			}
		}
		CountsHelper helperSNV = VariantDiscoverySNVQAlgorithm.calculateCountsSNV(pileup,maxBaseQS, null);
		if(pileup.getPosition()==posPrint) System.out.println("A count: "+helperSNV.getCount("A")+" total: "+helperSNV.getTotalCount() );
		if(pileup.getPosition()==posPrint) System.out.println("C count: "+helperSNV.getCount("C")+" total: "+helperSNV.getTotalCount() );
		if(pileup.getPosition()==posPrint) System.out.println("G count: "+helperSNV.getCount("G")+" total: "+helperSNV.getTotalCount() );
		if(pileup.getPosition()==posPrint) System.out.println("T count: "+helperSNV.getCount("T")+" total: "+helperSNV.getTotalCount() );
		GenomicVariant variant;
		if(referenceAllele.length()>1) {
			CountsHelper helperIndel = VariantDiscoverySNVQAlgorithm.calculateCountsIndel(pileup,null,referenceAllele, null); 
			variant = findMultiallelicIndel(pileup, helperIndel);
			if(variant!=null) {
				//System.out.println("Called indel at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+" variant type: "+calledVar.getType());
				lastIndelEnd = variant.getLast();
			} else {
				if (pileup.isNewSTR()) {
					pileup.setSTR(false);
					pileup.setNewSTR(false);
				}
				//Try SNV if the indel alleles were not good to make a call
				variant = findMultiallelicSNV(pileup, helperSNV, referenceAllele.charAt(0));
			}
		} else {
			variant = findMultiallelicSNV(pileup, helperSNV, referenceAllele.charAt(0));
		}
		if(variant != null) {
			if(variant.isSNV() && (pileup.isEmbedded() || variant.getFirst()<=lastIndelEnd)) variant.setType(GenomicVariant.TYPE_EMBEDDED_SNV);
		}
		return variant;
	}
	
	
	private GenomicVariant findMultiallelicSNV(PileupRecord pileup, CountsHelper helper, char reference) {
		if(helper.getTotalCount()==0) {
			return null;
		}
		int refIdx = DNASequence.BASES_STRING.indexOf(reference);
		if(refIdx<0) {
			//N reference can in principle be handled but it generates  many non variant sites
			return null;
		}
		//Simple method based on relative counts. To improve later
		int [] counts = helper.getCounts();
		int sum = NumberArrays.getSum(counts); 
		if(pileup.getPosition()==posPrint) System.out.println("Refidx: "+refIdx+" sum: "+sum);
		boolean [] allelesSupported = new boolean [ counts.length];
		List<String> alleles = new ArrayList<>();
		alleles.add(DNASequence.BASES_ARRAY[refIdx]);
		for(int i=0;i<counts.length;i++) {
			allelesSupported[i]=counts[i]>0 && (double)counts[i]/(double)sum >=minAlleleFrequency;
			if(allelesSupported[i] && i!=refIdx) {
				alleles.add(DNASequence.BASES_ARRAY[i]);
			}
		}
		if(pileup.getPosition()==posPrint) System.out.println("Alleles: "+alleles);
		GenomicVariant variant = null;
		if(alleles.size()==2) {
			variant = new SNV(pileup.getSequenceName(), pileup.getPosition(), reference, alleles.get(1).charAt(0));
			variant.setType(GenomicVariant.TYPE_BIALLELIC_SNV);
		} else if (alleles.size()>2){
			variant = new GenomicVariantImpl(pileup.getSequenceName(), pileup.getPosition(), alleles);
			variant.setType(GenomicVariant.TYPE_MULTIALLELIC_SNV);
			while(true) {	
				List<CalledGenomicVariant> calls = genotypeVariant(variant, pileup, heterozygosityRate);
				GenomicVariant newVariant = makeNewVariant(variant, calls);
				if(newVariant!=variant) variant = newVariant;
				else break;
			}
		}
		
		return variant;
	}

	private GenomicVariant findMultiallelicIndel(PileupRecord pileup, CountsHelper helperIndel) {
		List<String> alleles = helperIndel.getAllelesList();
		if(alleles.size() == 1) return null;
		if(helperIndel.getTotalCount()==0) return null;
		GenomicVariant variant = new GenomicVariantImpl(pileup.getSequenceName(),pileup.getPosition(),alleles);
		
		// Redefine indel alleles based on genotype calls
		while(true) {
			if(!pileup.isInputSTR() && allelesSameLength(variant.getAlleles())) return null;
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
	private boolean allelesSameLength(String[] alleles) {
		int l = alleles[0].length();
		for(String allele:alleles) {
			if(allele.length()!=l) return false;
		}
		return true;
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
		int n = samples.size();
		short variantQS = 0;
		for(int i=0;i<n;i++) {
			Sample sample = samples.get(i);
			CalledGenomicVariant call = genotypeVariantSample(variant, pileup, sample, h);
			if(pileup.getPosition()==posPrint) System.out.println("Sample: "+call.getSampleId()+" Genotype: "+Arrays.asList(call.getCalledAlleles())+" GQ: "+call.getGenotypeQuality());
			if(!call.isUndecided() && !call.isHomozygousReference() && call.getGenotypeQuality()>variantQS) {
				variantQS = call.getGenotypeQuality();
			}
			calls.add(call);
		}
		variant.setVariantQS(variantQS);
		return calls;
	}
	
	private CalledGenomicVariant genotypeVariantSample(GenomicVariant variant, PileupRecord pileup,  Sample sample, double h) {
		String referenceAllele = variant.getReference();
		
		CalledGenomicVariant calledVar = null;
		if(variant.isSNV()) {
			CountsHelper helperSNV = VariantDiscoverySNVQAlgorithm.calculateCountsSNV(pileup, maxBaseQS, sample.getReadGroups());
			calledVar = VariantDiscoverySNVQAlgorithm.callSNV(pileup, helperSNV, variant, referenceAllele.charAt(0), h, false);
		} else {
			CountsHelper helperIndel = VariantDiscoverySNVQAlgorithm.calculateCountsIndel(pileup,variant,referenceAllele, sample.getReadGroups()); 
			calledVar = VariantDiscoverySNVQAlgorithm.callIndel(pileup, helperIndel, variant, h, false);
		}
		if(calledVar==null) {
			calledVar = new CalledGenomicVariantImpl(variant, new byte[0]);
		}
		calledVar.setSampleId(sample.getId());
		calledVar.updateAllelesCopyNumberFromCounts(sample.getNormalPloidy());
		return calledVar;
	}

	private void dispose() {
		inputVariants =null;
		seqInputVariants = null;
		indelRealigner.setInputVariants(null);
	}


}
