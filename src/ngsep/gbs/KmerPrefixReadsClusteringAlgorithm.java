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
package ngsep.gbs;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import ngsep.alignments.ReadAlignment;
import ngsep.discovery.CountsHelper;
import ngsep.discovery.PileupRecord;
import ngsep.discovery.VariantDiscoverySNVQAlgorithm;
import ngsep.discovery.VariantPileupListener;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.math.NumberArrays;
import ngsep.math.PhredScoreHelper;
import ngsep.sequences.DNASequence;
import ngsep.sequences.DNAShortKmer;
import ngsep.sequences.DNAShortKmerClusterMap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;
import ngsep.sequencing.ReadsDemultiplex;
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

/**
 * @author Jorge Gomez
 * @author Jorge Duitama
 * @author Andrea Parra
 */
public class KmerPrefixReadsClusteringAlgorithm {

	private Logger log = Logger.getLogger(ReadsDemultiplex.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	public static final double DEF_MIN_HETEROZYGOSITY_RATE_DIPLOID = VariantPileupListener.DEF_HETEROZYGOSITY_RATE_DIPLOID;
	public static final short DEF_MIN_QUALITY = 40;
	public static final double DEF_MIN_ALLELE_FREQUENCY = 0;
	public static final byte DEF_MAX_BASE_QS = VariantPileupListener.DEF_MAX_BASE_QS;
	public static final int DEF_KMER_LENGTH = 31;
	public static final int DEF_START = 8;
	public static final int DEF_MAX_READS_IN_MEMORY = 1000000;
	public static final String DEF_REGEXP_SINGLE="<S>.fastq.gz";
	public static final String DEF_REGEXP_PAIRED="<S>_<N>.fastq.gz";
	public static final byte DEF_PLOIDY = GenomicVariant.DEFAULT_PLOIDY;
	
	private String inputDirectory=".";
	private String outPrefix="./output";
	private int kmerLength = DEF_KMER_LENGTH;
	private Pattern regexp=Pattern.compile(DEF_REGEXP_SINGLE);
	private Map<String, String> filenamesBySampleId1=new HashMap<>();
	private Map<String, String> filenamesBySampleId2=new HashMap<>();
	private DNAShortKmerClusterMap kmersMap;
	
	private ProcessInfo processInfo = new ProcessInfo();
	private short minQuality = DEF_MIN_QUALITY;
	// FIXME What is this value?
	private int posPrint = -1;
	
	private double minAlleleFrequency = DEF_MIN_ALLELE_FREQUENCY;
	private double heterozygosityRate = DEF_MIN_HETEROZYGOSITY_RATE_DIPLOID;
	private byte maxBaseQS = DEF_MAX_BASE_QS;
	private byte normalPloidy = DEF_PLOIDY;
	private int numClustersWithVariants = 0;
	private int numClusteredFiles;
	private int numUnclusteredReads = 0;
	private int numTotalReads = 0;
	
	public static void main(String[] args) throws Exception {
		KmerPrefixReadsClusteringAlgorithm instance = new KmerPrefixReadsClusteringAlgorithm();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.inputDirectory = args[i++];
		instance.outPrefix = args[i++];
		instance.run();
	}
	
	public void setProgressNotifier(ProgressNotifier progressNotifier) { 
		this.progressNotifier = progressNotifier;
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	/**
	 * @return the kmerLength
	 */
	public int getKmerLength() {
		return kmerLength;
	}
	/**
	 * @param kmerLength the kmerLength to set
	 */
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
	}
	
	public void setKmerLength(String value) {
		setKmerLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	public void run() throws IOException {
		
//		processInfo.addTime(System.nanoTime(), "Load files start");
//		loadFilenamesAndSamples();
//		processInfo.addTime(System.nanoTime(), "Load files end");
//		processInfo.addTime(System.nanoTime(), "BuildKmersMap start");
//		log.info("Loaded "+filenamesBySampleId1.size()+" samples");
//		buildKmersMap();
//		processInfo.addTime(System.nanoTime(), "BuildKmersMap end");
//		processInfo.addTime(System.nanoTime(), "Cluster reads start");
//		log.info("Built kmers map with "+kmersMap.size()+" clusters");
//		List<String> clusteredReadsFilenames = clusterReads();
//		printStatistics("initial");
//		processInfo.addTime(System.nanoTime(), "Cluster reads end");
//		processInfo.addTime(System.nanoTime(), "Variant calling start");
		List<String> clusteredReadsFilenames = debug();
		this.numClusteredFiles = clusteredReadsFilenames.size();
		log.info("Clustered reads");
		callVariants(clusteredReadsFilenames);
		processInfo.addTime(System.nanoTime(), "Variant calling end");
		log.info("Called variants");
		printStatistics("final");
		log.info("Process finished");
	}
	
	private List<String> debug() {		
		log.info("Skipping to call variants");
		int stop = 26;
		String prefix = "run_9_clusteredReads_";
		String suffix = ".fastq.gz";
		List<String> clusteredReadsFilenames = new ArrayList<>();
		for(int i=0;i<=stop;i++) {
			clusteredReadsFilenames.add(prefix + i + suffix);
		}
		
		return clusteredReadsFilenames;
	}

	private void loadFilenamesAndSamples() {
		File[] files = (new File(inputDirectory)).listFiles();
		for(File f : files) {
			String filename = f.getName();
			//TODO: Use pattern
			int i = filename.indexOf(".fastq");
			if(i>=0) {
				String sampleId = filename.substring(0, i);
				filenamesBySampleId1.put(sampleId, f.getAbsolutePath());
			}
		}
	}
	
	
	public void buildKmersMap() throws IOException {
		kmersMap = new DNAShortKmerClusterMap();
		for(String filename:filenamesBySampleId1.values()) {
			addKmersFromFile(filename);
		}
	}

	private void addKmersFromFile(String filename) throws IOException {
		int readCount = 0;
		try (FastqFileReader openFile = new FastqFileReader(filename);) {
			Iterator<RawRead> reader = openFile.iterator();
			while(reader.hasNext()) {
				RawRead read = reader.next();
				String s = read.getSequenceString();
				if(DEF_START + kmerLength>s.length()) continue;
				String prefix = s.substring(DEF_START,DEF_START + kmerLength);
				if(DNASequence.isDNA(prefix)) {
					kmersMap.addOcurrance(new DNAShortKmer(prefix));
					readCount++;
				}
			}
		}
		log.info("Processed a total of " + readCount + " reads for file: "+filename);
		
	}
	public List<String> clusterReads() throws IOException {
		ClusteredReadsCache clusteredReadsCache = new ClusteredReadsCache(); 
		for(String sampleId:filenamesBySampleId1.keySet()) {
			String filename1 = filenamesBySampleId1.get(sampleId);
			String filename2 = filenamesBySampleId2.get(sampleId);
			if(filename2 == null) {
				log.info("Clustering reads from " + filename1);
				clusterReadsSingleFile (sampleId, filename1, clusteredReadsCache);
			} else {
				clusterReadsPairedEndFiles (sampleId, filename1, filename2, clusteredReadsCache);
			}
		}
		clusteredReadsCache.dump(outPrefix);
		return clusteredReadsCache.getClusteredReadFiles();
	}

	private void clusterReadsSingleFile(String sampleId, String filename, ClusteredReadsCache clusteredReadsCache) throws IOException {
		int unmatchedReads = 0;
		int count = 1;
		try (FastqFileReader openFile = new FastqFileReader(filename);) {
			Iterator<RawRead> reader = openFile.iterator();
			while(reader.hasNext()) {
				this.numTotalReads++;
				RawRead read = reader.next();
				String s = read.getSequenceString();
				if(DEF_START + kmerLength>s.length()) continue;
				String prefix = s.substring(DEF_START,DEF_START + kmerLength);
				if(!DNASequence.isDNA(prefix)) continue;
				Integer clusterId = kmersMap.getCluster(new DNAShortKmer(prefix));
				if(clusterId==null) {
					this.numUnclusteredReads++;
					unmatchedReads++;
					continue;
				}
				clusteredReadsCache.addSingleRead(clusterId, new RawRead(sampleId+"$"+clusterId+"$"+read.getName(), s, read.getQualityScores()));
				if(clusteredReadsCache.getTotalReads()>=DEF_MAX_READS_IN_MEMORY) {
					log.info("dumping reads");
					clusteredReadsCache.dump(outPrefix);
				}
				count++;
			}
			log.info(Integer.toString(unmatchedReads) + " reads remained unmatched for file: " + filename);
			log.info(Integer.toString(count) + " reads were succesfully matched for file: " + filename);
		}
	}

	private void clusterReadsPairedEndFiles(String sampleId, String filename1, String filename2, ClusteredReadsCache clusteredReadsCache) throws IOException {
		// TODO Auto-generated method stub
		return;
	}
	
	public void callVariants(List<String> clusteredReadsFilenames) throws IOException {
		int numberOfFiles = clusteredReadsFilenames.size();
		
		//process files in parallel
		FastqFileReader [] readers = new FastqFileReader[numberOfFiles];
		RawRead [] currentReads = new RawRead[numberOfFiles];
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		VCFFileWriter writer = new VCFFileWriter ();
		
		//TODO: Add sample ids to header
		Arrays.fill(currentReads, null);
		List<Iterator<RawRead>> iterators = new ArrayList<>();
		try (PrintStream outVariants = new PrintStream(outPrefix+"_variants.vcf");
				PrintStream clusterInfo = new PrintStream(outPrefix+"_clustersInfo.txt");
				PrintStream stream = new PrintStream(outPrefix+"_clustersVariantsInfoDEBUG_1.txt");) {
			int numNotNull = 0;
			int numCluster = 0;
			for(int i=0; i<numberOfFiles; i++) {
				readers[i] = new FastqFileReader(clusteredReadsFilenames.get(i));
				Iterator<RawRead> it = readers[i].iterator();
				iterators.add(it);
				if(it.hasNext()) {
					currentReads[i] = it.next();
					numNotNull++;
				}
			}
			stream.println("clusterNum\tvariatns1\tvariants2\trecords\tundecidedFlag\thomoxygousFlag\tgenotypeQSFlag");
			log.info("Processing a total of " + Integer.toString(numberOfFiles) + " clustered files.");
			clusterInfo.println("clusterId \t hammingDist \t totNumReads \t numVariants");
			while(numNotNull>0) {
				//gather reads next cluster
				ReadCluster nextCluster = new ReadCluster(numCluster);
				for(int i=0; i<numberOfFiles; i++) {
					addReadsToCluster(nextCluster, iterators.get(i), currentReads, i);
					if(currentReads[i]==null) numNotNull--;
				}
				List<VCFRecord> records = processCluster(nextCluster, header, stream);
				writer.printVCFRecords(records, outVariants);
				
				if(numCluster % 1000 == 0) {
					log.info("Done with cluster: " + Integer.toString(numCluster));
				}
				numCluster++;
			}
			
		} finally {
			for(FastqFileReader reader:readers) {
				if(reader!=null) reader.close();
			}
			// XXX Is it necessary to close the VCF writer?
		}
	}
		
	
	private void addReadsToCluster(ReadCluster nextCluster, Iterator<RawRead> iterator, 
			RawRead[] currentReads, int i) throws IOException {
		
		RawRead currentRead = currentReads[i];
		int numCluster = nextCluster.getClusterNumber();
		while(currentRead!=null) {
			String readIdWithCluster = currentRead.getName();
			String [] items = readIdWithCluster.split("\\$");
			if(items.length < 2) {
				continue;
			}
			String sampleId = items[0];
			int currentReadCluster = Integer.parseInt(items[1]);
			// TODO rename reads without cluster info
			String readId = items[2];
			if (currentReadCluster>numCluster) break;
			else if (currentReadCluster<numCluster) throw new RuntimeException("Disorgainzed file. Current cluster: "+numCluster+" found: "+currentReadCluster+" in read. "+readIdWithCluster);
			nextCluster.addRead(currentRead, sampleId);
			if(iterator.hasNext()) {
				currentReads[i] = iterator.next();
				currentRead = currentReads[i];
				
			} else {
				log.info("Done with file " + Integer.toString(i) + ".");
				currentReads[i] = null;
				currentRead = null;				
			}
		}
	}

	private List<VCFRecord> processCluster(ReadCluster readCluster, VCFFileHeader vcfFileHeader, PrintStream stream) throws IOException {
		
		List<VCFRecord> records = new ArrayList<>();
		List<ReadAlignment> readAlignments = new ArrayList<>();
		List<Sample> samples = new ArrayList<>();
		int clusterId = readCluster.getClusterNumber();
		String refSeq = readCluster.getRefSeq();
		String referenceId = Integer.toString(clusterId);
		QualifiedSequence refQS = new QualifiedSequence(referenceId, refSeq);
		ReferenceGenome singleSequenceGenome = new ReferenceGenome (refQS);
		VariantPileupListener variantsDetector = new VariantPileupListener();
		variantsDetector.setGenome(singleSequenceGenome);
		
		// For each read within the cluster create a ReadAlignment. Set characters and quality scores
		List<RawRead> reads = readCluster.getReads();
		List<String> sampleIds = readCluster.getSampleIds();
		
		// based on cluster depth, look for sub-clusters
		if(reads.size()>5000) {
//			List<ReadCluster> subclusters = findSubClusters(readCluster);
//			for(ReadCluster subcluster:subclusters) {
//				processCluster(subcluster, vcfFileHeader, stream);
//			}
			return records;

		}
		for(int i=0;i<reads.size();i++) {
			RawRead read = reads.get(i);
			String sampleId = sampleIds.get(i);
			Sample sample = new Sample(sampleId);
			sample.setNormalPloidy(normalPloidy);
			sample.addReadGroup(sampleId);
			int readLength = read.getLength();
			String CIGARString = Integer.toString(readLength) + "M"; 
			ReadAlignment readAlignment = new ReadAlignment(referenceId, 1, readLength, readLength, 0);
			readAlignment.setQualityScores(read.getQualityScores());
			readAlignment.setReadCharacters(read.getCharacters());
			readAlignment.setReadName(read.getName());
			readAlignment.setCigarString(CIGARString);
			readAlignment.setReadGroup(sampleId);
			readAlignments.add(readAlignment);	
			vcfFileHeader.addSample(sample, false);
			samples.add(sample);
			}

		// For each position in the representative sequence create a pileup record with cluster id as sequence name and position =i
		// FIXME Start at 1? if not, no readAlignments are added to the pileup CHECK.
		System.out.println("Cluster " + referenceId);
		for(int i=1; i<=refSeq.length(); i++) {
			
			PileupRecord clusterPileUp = new PileupRecord(referenceId, i);
			for(ReadAlignment readAlgn:readAlignments) {
				clusterPileUp.addAlignment(readAlgn);
			}
			
			GenomicVariant variant = findMultiallelicVariant(samples, clusterPileUp, refSeq.charAt(i-1), referenceId, heterozygosityRate);
			if(variant!=null) {
				List<CalledGenomicVariant> calls = genotypeVariant(samples, variant, clusterPileUp, heterozygosityRate);
				if(variant.getVariantQS()==0 || variant.getVariantQS() < minQuality) continue;
				DiversityStatistics divStats = DiversityStatistics.calculateDiversityStatistics(calls, false);
				int [] format = variant.isSNV()?VCFRecord.DEF_FORMAT_ARRAY_NGSEP_SNV:VCFRecord.DEF_FORMAT_ARRAY_NGSEP_NOSNV;
				VCFRecord record = new VCFRecord(variant, format, calls, vcfFileHeader);
				record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_SAMPLES_GENOTYPED, divStats.getNumSamplesGenotyped()));
				record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_NUMBER_ALLELES, divStats.getNumCalledAlleles()));
				//record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_ALLELE_FREQUENCY_SPECTRUM, format(divStats.getAlleleCounts())));
				if(variant.isBiallelic()) record.addAnnotation(new GenomicVariantAnnotation(variant, GenomicVariantAnnotation.ATTRIBUTE_MAF, divStats.getMaf()));
				records.add(record);
			}
		}
		
		if(records.size()>0) {
			System.out.println(records.size() +" records generated for cluster " + referenceId);	
		}		
		return records;
	}
	
	private List<ReadCluster> findSubClusters(ReadCluster readCluster) {
		List<ReadCluster> subclusters = new ArrayList<>();
		kmersMap = new DNAShortKmerClusterMap();
		List<RawRead> reads = readCluster.getReads();
		
		return subclusters;
	}

	// From MultisampleVariantsDetector.java
	private GenomicVariant findMultiallelicVariant(List<Sample> samples, PileupRecord clusterPileUp, char reference, String clusterNum, double h) {
		
		String referenceAllele = "" + reference;
		referenceAllele = referenceAllele.toUpperCase();
		
		CountsHelper helperSNV = VariantDiscoverySNVQAlgorithm.calculateCountsSNV(clusterPileUp,maxBaseQS, null);
		if(clusterPileUp.getPosition()==posPrint) System.out.println("A count: "+helperSNV.getCount("A")+" total: "+helperSNV.getTotalCount() );
		if(clusterPileUp.getPosition()==posPrint) System.out.println("C count: "+helperSNV.getCount("C")+" total: "+helperSNV.getTotalCount() );
		if(clusterPileUp.getPosition()==posPrint) System.out.println("G count: "+helperSNV.getCount("G")+" total: "+helperSNV.getTotalCount() );
		if(clusterPileUp.getPosition()==posPrint) System.out.println("T count: "+helperSNV.getCount("T")+" total: "+helperSNV.getTotalCount() );
		GenomicVariant variant;
		variant = findMultiallelicSNV(samples, clusterPileUp, helperSNV, referenceAllele.charAt(0), h);
		return variant;
	}
	
	private GenomicVariant findMultiallelicSNV(List<Sample> samples, PileupRecord pileup, CountsHelper helper, char reference, double h) {
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
				List<CalledGenomicVariant> calls = genotypeVariant(samples, variant, pileup, h);
				GenomicVariant newVariant = makeNewVariant(variant, calls);
				if(newVariant!=variant) variant = newVariant;
				else break;
			}
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
		int n = calls.size();
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
	
	
	private List<CalledGenomicVariant> genotypeVariant(List<Sample> samples, GenomicVariant variant, PileupRecord pileup, double h) {

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
			CountsHelper helperIndel = VariantDiscoverySNVQAlgorithm.calculateCountsIndel(pileup,variant,referenceAllele, maxBaseQS, sample.getReadGroups()); 
			calledVar = VariantDiscoverySNVQAlgorithm.callIndel(pileup, helperIndel, variant, h, false);
		}
		if(calledVar==null) {
			calledVar = new CalledGenomicVariantImpl(variant, new byte[0]);
		}
		calledVar.setSampleId(sample.getId());
		calledVar.updateAllelesCopyNumberFromCounts(sample.getNormalPloidy());
		return calledVar;
	}

	private void printStatistics(String ref) throws IOException {
		try(PrintStream processStats = new PrintStream(outPrefix+"_"+ref+"_processInfo.txt");){
			processStats.println("Process times:");
			processStats.print(this.processInfo.getTimes());
			processStats.println("Number of Files: " + Integer.toString(this.filenamesBySampleId1.size()));
			processStats.println("Number of Cluster Files: " + Integer.toString(this.numClusteredFiles));
			processStats.println("Number of Clusters: " + Integer.toString(kmersMap.size()));
			processStats.println("Number of Clusters with variants: " + Integer.toString(this.numClustersWithVariants));
			processStats.println("Number of Reads: " + Integer.toString(this.numTotalReads));
			processStats.println("Number of Unclustered Reads: " + Integer.toString(this.numUnclusteredReads));
		}
	}
	
	/**
	 * 
	 */
	public void readClustersFromFile(String filename) throws IOException {
		int readCount = 0;
		try (FastqFileReader openFile = new FastqFileReader(filename);) {
			Iterator<RawRead> reader = openFile.iterator();
			while(reader.hasNext()) {
				RawRead read = reader.next();
				String s = read.getSequenceString();
				if(DEF_START + kmerLength>s.length()) continue;
				String prefix = s.substring(DEF_START,DEF_START + kmerLength);
				if(DNASequence.isDNA(prefix)) {
					kmersMap.addOcurrance(new DNAShortKmer(prefix));
					readCount++;
				}
			}
		}
		log.info("Processed a total of " + readCount + " reads for file: "+filename);	
	}
}

class ProcessInfo {
	private List<Long> processTimes = new ArrayList<>(); 
	private List<String> processTimesLabels = new ArrayList<>();
	
	public void addTime(long time, String label) {
		processTimes.add(time);
		processTimesLabels.add(label);
	}
	
	public String getTimes() {
		String timeInfo = "";
		for(int i=0; i<processTimes.size(); i++) {
			timeInfo += processTimesLabels.get(i) + ": " + Long.toString(processTimes.get(i)) + "\n";
		}
		return timeInfo;
	}
}

class ClusteredReadsCache {
	private Map<Integer,List<RawRead>> clusteredReadsCache = new TreeMap<>();
	private Map<Integer, Map<String, String>> clusterInfo = new TreeMap<>();
	private int totalReads = 0;
	private List<String> outFiles = new ArrayList<>();
	
	public void addSingleRead(int k, RawRead read) {
		List<RawRead> readsClusterK = clusteredReadsCache.get(k);
		if(readsClusterK==null) {
			readsClusterK = new ArrayList<>();
			clusteredReadsCache.put(k, readsClusterK);
		}
		readsClusterK.add(read);
		totalReads++;
	}
	
	public void addClusterInfo(int k, String key, String value) {
		Map<String, String> info = clusterInfo.get(k);
		if(clusterInfo == null) {
			info = new HashMap<>();		
		}
		info.put(key, value);
	}
	
	public Map<String, String> getClusterInfo(int k) {
		return clusterInfo.get(k);
	}
	
	public List<String> getClusteredReadFiles() {
		return outFiles;
	}
	
	public List<RawRead> getClusterReads(int k) {
		return clusteredReadsCache.get(k);
	}
	
	/**
	 * @return the totalReads
	 */
	public int getTotalReads() {
		return totalReads;
	}
	
	/**
	 * Dumps the cache to the file with the given name and clears this cache
	 * @param outPrefix prefix of the file to dump the cache
	 */
	public void dump(String outPrefix) throws IOException {
		int number = outFiles.size();
		String singleFilename = outPrefix+"_clusteredReads_"+number+".fastq.gz";
		outFiles.add(singleFilename);
		try (OutputStream os = new FileOutputStream(singleFilename);
			 GZIPOutputStream gos = new GZIPOutputStream(os);
			 PrintStream out = new PrintStream(gos);) {
			for(List<RawRead> readsCluster:clusteredReadsCache.values()) {
				for(RawRead read:readsCluster) {
					read.save(out);
				}
			}
		}
		clusteredReadsCache.clear();
		totalReads = 0;
	}
}
