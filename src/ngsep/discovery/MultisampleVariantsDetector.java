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
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.DiversityStatistics;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantAnnotation;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class MultisampleVariantsDetector implements PileupListener {

	private Logger log = Logger.getLogger(MultisampleVariantsDetector.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	public static final int DEF_MAX_ALNS_PER_START_POS = AlignmentsPileupGenerator.DEF_MAX_ALNS_PER_START_POS;
	public static final double DEF_MIN_ALLELE_FREQUENCY = 0.05;
	public static final double DEF_MIN_HETEROZYGOSITY_RATE_DIPLOID = VariantPileupListener.DEF_HETEROZYGOSITY_RATE_DIPLOID;
	public static final short DEF_MIN_QUALITY = VariantPileupListener.DEF_MIN_QUALITY;
	public static final short DEF_MIN_MQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	public static final short DEF_MAX_BASE_QS = VariantPileupListener.DEF_MAX_BASE_QS;
	public static final byte DEF_PLOIDY = GenomicVariant.DEFAULT_PLOIDY;
	
	
	
	private double coveredGenomeSize = 0;
	private long referenceGenomeSize = 0;
	
	private AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
	//Listeners
	private IndelRealignerPileupListener indelRealigner = new IndelRealignerPileupListener();

	
	private ReferenceGenome genome;
	private List<String> alignmentFiles = new ArrayList<>();
	private List<Sample> samples;
	private List<VariantPileupListener> genotypingListeners = new ArrayList<>();
	
	//Output file variables
	private String outFilename;
	private PrintStream outFile;
	private VCFFileHeader vcfFileHeader;
	private VCFFileWriter writer = new VCFFileWriter();
	
	private double heterozygosityRate = DEF_MIN_HETEROZYGOSITY_RATE_DIPLOID;
	private String knownSTRsFile = null; 
	private boolean ignoreLowerCaseRef = false;
	private boolean callEmbeddedSNVs = false;
	private double minAlleleFrequency = DEF_MIN_ALLELE_FREQUENCY;
	private short minQuality = DEF_MIN_QUALITY;
	private short maxBaseQS = DEF_MAX_BASE_QS;
	
	
	//Control attribute to avoid calling overlapping indels and to give an embedded status to SNVs within indels or STRs
	private int lastIndelEnd = 0;
	
	//DEBUG
	private int posPrint = -1;
	
	public static void main(String[] args) throws Exception {
		MultisampleVariantsDetector instance = new MultisampleVariantsDetector();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.genome = new ReferenceGenome(args[i++]);
		instance.outFilename = args[i++];
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
	public short getMaxBaseQS() {
		return maxBaseQS;
	}

	/**
	 * @param maxBaseQS the maxBaseQS to set
	 */
	public void setMaxBaseQS(short maxBaseQS) {
		this.maxBaseQS = maxBaseQS;
	}
	
	public void setMaxBaseQS(String value) {
		setMaxBaseQS((short)OptionValuesDecoder.decode(value, Short.class));
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

	public void findVariants() throws IOException {
		referenceGenomeSize = genome.getTotalLength();
		//TODO: assign sample ids if not in aln files
		if(samples == null) loadSamplesFromAlignmentHeaders();
		genotypingListeners.clear();
		for(Sample sample:samples) {
			VariantPileupListener genotypingListener = new VariantPileupListener();
			genotypingListener.setCallEmbeddedSNVs(callEmbeddedSNVs);
			genotypingListener.setGenome(genome);
			genotypingListener.setHeterozygosityRate(heterozygosityRate);
			genotypingListener.setIgnoreLowerCaseRef(ignoreLowerCaseRef);
			genotypingListener.setMaxBaseQS(maxBaseQS);
			genotypingListener.setMinQuality(minQuality);
			genotypingListener.setNormalPloidy(sample.getNormalPloidy());
			genotypingListener.setReadGroups(sample.getReadGroups());
			genotypingListeners.add(genotypingListener);
		}
		if(knownSTRsFile!=null) {
			log.info("Loading input short tandem repeats from: "+knownSTRsFile);
			//TODO: Choose the best format
			SimpleGenomicRegionFileHandler rfh = new SimpleGenomicRegionFileHandler();
			List<GenomicRegion> strs = rfh.loadRegions(knownSTRsFile);
			//TODO: STRs loader
			indelRealigner.setInputVariants(VariantsDetector.makeNonRedundantSTRs(genome,strs));
			log.info("Loaded "+strs.size()+" input short tandem repeats");
		}
		log.info("Finding variants");
		indelRealigner.setGenome(genome);
		generator.setSequencesMetadata(genome.getSequencesMetadata());
		generator.addListener(indelRealigner);
		generator.addListener(this);
		try {
			outFile = new PrintStream(outFilename);
			vcfFileHeader = VCFFileHeader.makeDefaultEmptyHeader();
			for(Sample s:samples) vcfFileHeader.addSample(s, true);
			writer.printHeader(vcfFileHeader, outFile);
			generator.processFiles(alignmentFiles);
		} finally {
			if(outFile!=null) outFile.close();
		}
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
						samplesMap.put(sampleId, sample);
						log.info("Found sample: "+sampleId+" in file: "+filename);
					}
					sample.addReadGroup(readGroup);
				}
			}
		}
		samples = new ArrayList<>(samplesMap.values()); 
	}

	@Override
	public void onPileup(PileupRecord pileup) {
		GenomicVariant variant = findVariant(pileup);
		if(pileup.getPosition()==posPrint) System.out.println("Variant: "+variant);
		if(variant == null) return;
		List<CalledGenomicVariant> calls;
		while(true) {
			calls = genotypeVariant(variant, pileup);
			int n = samples.size();
			if(pileup.getPosition()==posPrint) System.out.println("Num samples: "+n);
			short variantQS = 0;
			Set<String> calledAllelesSet = new TreeSet<>();
			calledAllelesSet.add(variant.getReference());
			for(int i=0;i<n;i++) {
				CalledGenomicVariant call = calls.get(i);
				if(!call.isUndecided() && !call.isHomozygousReference() && call.getGenotypeQuality()>variantQS) {
					variantQS = call.getGenotypeQuality();
				}
				calledAllelesSet.addAll(Arrays.asList(call.getCalledAlleles()));
			}
			if(variantQS == 0) return;
			if(variant.getAlleles().length>2 && variant.getAlleles().length !=calledAllelesSet.size()) {
				variant = makeNewVariant(variant,calledAllelesSet);
			} else {
				variant.setVariantQS(variantQS);
				break;
			}
		}
		
		DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(calls, false);
		int [] format = variant.isSNV()?VCFRecord.DEF_FORMAT_ARRAY_NGSEP_SNV:VCFRecord.DEF_FORMAT_ARRAY_NGSEP_NOSNV;
		VCFRecord record = new VCFRecord(variant, format, calls, vcfFileHeader);
		record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_SAMPLES_GENOTYPED, divStats.getNumSamplesGenotyped()));
		record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_NUMBER_ALLELES, divStats.getNumCalledAlleles()));
		record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_ALLELE_FREQUENCY_SPECTRUM, format(divStats.getAlleleCounts())));
		if(divStats.getNumCalledAlleles()==2) record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_MAF, divStats.getMaf()));
		
		writer.printVCFRecord(record, outFile);
		coveredGenomeSize++;
		if(progressNotifier!=null && coveredGenomeSize%10000==0) {
			int progress = (int)Math.round(100.0*coveredGenomeSize/referenceGenomeSize);
			generator.setKeepRunning(progressNotifier.keepRunning(progress));
		}
	}
	

	private String format(int[] alleleCounts) {
		StringBuilder answer = new StringBuilder(""+alleleCounts[0]);
		for(int i=1;i<alleleCounts.length;i++) answer.append(","+alleleCounts[i]);
		return answer.toString();
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

	private GenomicVariant findVariant(PileupRecord pileup) {
		
		if(!callEmbeddedSNVs && pileup.isEmbedded()) return null;
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
		CountsHelper helperSNV = VariantPileupListener.calculateCountsSNV(pileup,maxBaseQS, null);
		if(pileup.getPosition()==posPrint) System.out.println("A count: "+helperSNV.getCount("A")+" total: "+helperSNV.getTotalCount() );
		if(pileup.getPosition()==posPrint) System.out.println("C count: "+helperSNV.getCount("C")+" total: "+helperSNV.getTotalCount() );
		if(pileup.getPosition()==posPrint) System.out.println("G count: "+helperSNV.getCount("G")+" total: "+helperSNV.getTotalCount() );
		if(pileup.getPosition()==posPrint) System.out.println("T count: "+helperSNV.getCount("T")+" total: "+helperSNV.getTotalCount() );
		GenomicVariant variant;
		if(referenceAllele.length()>1) {
			CountsHelper helperIndel = VariantPileupListener.calculateCountsIndel(pileup,null,referenceAllele, null); 
			variant = callMultisampleIndel(pileup, helperIndel);
			if(variant!=null) {
				//System.out.println("Called indel at "+calledVar.getSequenceName()+":"+calledVar.getFirst()+" variant type: "+calledVar.getType());
				lastIndelEnd = variant.getLast();
			} else {
				if (pileup.isNewSTR()) {
					pileup.setSTR(false);
					pileup.setNewSTR(false);
				}
				//Try SNV if the indel alleles were not good to make a call
				variant = callMultisampleSNV(pileup, helperSNV, referenceAllele.charAt(0));
			}
		} else {
			variant = callMultisampleSNV(pileup, helperSNV, referenceAllele.charAt(0));
		}
		if(variant != null) {
			//System.out.println("Called SNV");
			if(variant.isSNV() && (pileup.isEmbedded() || variant.getFirst()<=lastIndelEnd)) variant.setType(GenomicVariant.TYPE_EMBEDDED_SNV);
		}
		return variant;
	}
	
	
	private GenomicVariant callMultisampleSNV(PileupRecord pileup, CountsHelper helper, char reference) {
		if(helper.getTotalCount()==0) {
			return null;
		}
		if(DNASequence.BASES_STRING.indexOf(reference)<0) {
			//N reference can in principle be handled but it generates  many non variant sites
			return null;
		}
		//Simple method based on relative counts. To improve later
		int [] counts = helper.getCounts();
		int sum = NumberArrays.getSum(counts); 
		int refIdx = DNASequence.BASES_STRING.indexOf(reference);
		if(pileup.getPosition()==posPrint) System.out.println("Refidx: "+refIdx+" sum: "+sum);
		boolean [] allelesSupported = new boolean [ counts.length];
		List<String> alleles = new ArrayList<>();
		if(refIdx>=0) alleles.add(DNASequence.BASES_ARRAY[refIdx]);
		else alleles.add(""+reference);
		for(int i=0;i<counts.length;i++) {
			allelesSupported[i]=(double)counts[i]/(double)sum >=minAlleleFrequency;
			if(allelesSupported[i] && i!=refIdx) {
				alleles.add(DNASequence.BASES_ARRAY[i]);
			}
		}
		if(pileup.getPosition()==posPrint) System.out.println("Alleles: "+alleles);
		if(alleles.size()==1) return null;
		GenomicVariant variant;
		if(alleles.size()==2) {
			variant = new SNV(pileup.getSequenceName(), pileup.getPosition(), reference, alleles.get(1).charAt(0));
			variant.setType(GenomicVariant.TYPE_BIALLELIC_SNV);
			//variant.setVariantQS(variantQS);
		} else {
			variant = new GenomicVariantImpl(pileup.getSequenceName(), pileup.getPosition(), alleles);
			variant.setType(GenomicVariant.TYPE_MULTIALLELIC_SNV);
			//variant.setVariantQS(variantQS);
		}
		
		return variant;
	}

	private GenomicVariant callMultisampleIndel(PileupRecord pileup, CountsHelper helperIndel) {
		List<String> alleles = helperIndel.getAllelesList();
		if(alleles.size() == 1) return null;
		GenomicVariantImpl variant = new GenomicVariantImpl(pileup.getSequenceName(),pileup.getPosition(),alleles);
		if (pileup.isSTR()) {
			variant.setType(GenomicVariant.TYPE_STR);
		} else {
			variant.setType(GenomicVariant.TYPE_INDEL);
		}
		return variant;
	}
	
	private List<CalledGenomicVariant> genotypeVariant(GenomicVariant variant, PileupRecord pileup) {
		List<CalledGenomicVariant> calls = new ArrayList<>();
		int n = samples.size();
		for(int i=0;i<n;i++) {
			Sample sample = samples.get(i);
			VariantPileupListener genotyper = genotypingListeners.get(i);
			CalledGenomicVariant call = genotyper.processPileup(pileup, variant);
			if(call == null) call = new CalledGenomicVariantImpl(variant, new byte[0]);
			//if(pileup.getPosition()==posPrint) System.out.println("Genotype sample "+sample.getId()+": "+((CalledSNV)call).getGenotype());
			call.setSampleId(sample.getId());
			calls.add(call);
		}
		return calls;
	}

	@Override
	public void onSequenceStart(QualifiedSequence sequence) {
		lastIndelEnd = 0;
	}

	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {
		
	}

}
