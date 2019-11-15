package ngsep.gbs;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import ngsep.alignments.ReadAlignment;
import ngsep.discovery.CountsHelper;
import ngsep.discovery.PileupRecord;
import ngsep.discovery.VariantDiscoverySNVQAlgorithm;
import ngsep.discovery.VariantPileupListener;
import ngsep.genome.ReferenceGenome;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;
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

public class ProcessClusterVCFTask extends Thread {
	//Results
	private boolean hasFinished = false;
	
	//Data
	private ReadCluster readCluster;
	private VCFFileHeader vcfFileHeader;
	private VCFFileWriter vcfWriter;
	private PrintStream outVariants;
	private PrintStream outConsensus;
	
	private KmerPrefixReadsClusteringAlgorithm parent;
	
	
	public ProcessClusterVCFTask(ReadCluster readCluster, VCFFileHeader vcfFileHeader, VCFFileWriter writer, KmerPrefixReadsClusteringAlgorithm parent, PrintStream outVariants) {
		this.readCluster = readCluster;
		this.vcfFileHeader = vcfFileHeader;
		this.vcfWriter = writer;
		this.outVariants = outVariants;
		this.parent = parent;
	}
	
	/**
	 * 
	 * @return true if the task is finished, false otherwise
	 */
	public boolean hasFinished() {
		return hasFinished;
	}
	
	@Override
	public void run() {
		List<VCFRecord> generatedRecords = generateRecordsForCluster();
		
		synchronized (outConsensus) {
			generateConsensusFasta();
		}
		
		//Writing synchronously to the centralized vcf writter
		synchronized (vcfWriter) {
			vcfWriter.printVCFRecords(generatedRecords, outVariants);
		}
		
		//Writing synchronously to statistics
		synchronized (parent) {
			parent.countVariants(generatedRecords);
		}
		
		hasFinished = true;
	}
	
	private List<VCFRecord> generateRecordsForCluster() {
		List<VCFRecord> records = new ArrayList<>();
		List<ReadAlignment> readAlignments = new ArrayList<>();
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
		
		for(int i=0;i<reads.size();i++) {
			RawRead read = reads.get(i);
			String sampleId = sampleIds.get(i);
			int readLength = read.getLength();
			String CIGARString = Integer.toString(readLength) + "M"; 
			ReadAlignment readAlignment = new ReadAlignment(referenceId, 1, readLength, readLength, 0);
			readAlignment.setQualityScores(read.getQualityScores());
			readAlignment.setReadCharacters(read.getCharacters());
			readAlignment.setReadName(read.getName());
			readAlignment.setCigarString(CIGARString);
			readAlignment.setReadGroup(sampleId);
			readAlignments.add(readAlignment);	
		}

		// For each position in the representative sequence create a pileup record with cluster id as sequence name and position =i
		for(int i=1; i<=refSeq.length(); i++) {
			
			PileupRecord clusterPileUp = new PileupRecord(referenceId, i);
			for(ReadAlignment readAlgn:readAlignments) {
				clusterPileUp.addAlignment(readAlgn);
			}
			double h = parent.getHeterozygosityRate();
			List<Sample> samples = parent.getSamples();
			GenomicVariant variant = findMultiallelicVariant(samples, clusterPileUp, refSeq.charAt(i-1), referenceId, h);
			if(variant!=null) {
				List<CalledGenomicVariant> calls = genotypeVariant(samples, variant, clusterPileUp, h);
				if(variant.getVariantQS()==0 || variant.getVariantQS() < parent.getMinQuality()) continue;
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
		
		return records;
	}
	
	private void generateConsensusFasta() {
		outConsensus.println(">Cluster_" + readCluster.getClusterNumber());
		outConsensus.println(readCluster.getRefSeq());
	}
	// From MultisampleVariantsDetector.java
	private GenomicVariant findMultiallelicVariant(List<Sample> samples, PileupRecord clusterPileUp, char reference, String clusterNum, double h) {
		
		String referenceAllele = "" + reference;
		referenceAllele = referenceAllele.toUpperCase();
		
		CountsHelper helperSNV = VariantDiscoverySNVQAlgorithm.calculateCountsSNV(clusterPileUp,parent.getMaxBaseQS(), null);
		//posPrint used to be here
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
		//posPrint used to be here
		boolean [] allelesSupported = new boolean [ counts.length];
		List<String> alleles = new ArrayList<>();
		alleles.add(DNASequence.BASES_ARRAY[refIdx]);
		for(int i=0;i<counts.length;i++) {
			allelesSupported[i]=counts[i]>0 && (double)counts[i]/(double)sum >=parent.getMinAlleleFrequency();
			if(allelesSupported[i] && i!=refIdx) {
				alleles.add(DNASequence.BASES_ARRAY[i]);
			}
		}
		
		//posPrint used to be here
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
		//posPrint used to be here
		List<String> alleles = new ArrayList<>(newAlleles.size());
		String reference = variant.getReference(); 
		alleles.add(reference);
		for(String allele:newAlleles) {
			if(!allele.equals(reference)) alleles.add(allele);
		}

		//posPrint used to be here
		if(variant.isSNV() && alleles.size()==2) {
			return new SNV(variant.getSequenceName(), variant.getFirst(), reference.charAt(0), alleles.get(1).charAt(0));
		}
		return new GenomicVariantImpl(variant.getSequenceName(), variant.getFirst(), alleles);
	}
	
	private List<CalledGenomicVariant> genotypeVariant(List<Sample> samples, GenomicVariant variant, PileupRecord pileup, double h) {
		//posPrint used to be here
		List<CalledGenomicVariant> calls = new ArrayList<>();
		int n = samples.size();
		short variantQS = 0;
		for(int i=0;i<n;i++) {
			Sample sample = samples.get(i);
			CalledGenomicVariant call = genotypeVariantSample(variant, pileup, sample, h);
			
			//posPrint used to be here
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
			CountsHelper helperSNV = VariantDiscoverySNVQAlgorithm.calculateCountsSNV(pileup, parent.getMaxBaseQS(), sample.getReadGroups());
			calledVar = VariantDiscoverySNVQAlgorithm.callSNV(pileup, helperSNV, variant, referenceAllele.charAt(0), h, false);
		} else {
			CountsHelper helperIndel = VariantDiscoverySNVQAlgorithm.calculateCountsIndel(pileup,variant,referenceAllele, parent.getMaxBaseQS(), sample.getReadGroups()); 
			calledVar = VariantDiscoverySNVQAlgorithm.callIndel(pileup, helperIndel, variant, h, false);
		}
		if(calledVar==null) {
			calledVar = new CalledGenomicVariantImpl(variant, new byte[0]);
		}
		calledVar.setSampleId(sample.getId());
		calledVar.updateAllelesCopyNumberFromCounts(sample.getNormalPloidy());
		return calledVar;
	}
}
