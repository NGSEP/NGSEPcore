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

/**
 * 
 * @author Jorge Gomez
 * @author Andrea Parra
 *
 */
public class ProcessClusterVCFTask extends Thread {
	//Results
	private boolean hasFinished = false;
	
	//Data
	private ReadCluster readCluster;
	private VCFFileHeader vcfFileHeader;
	private VCFFileWriter vcfWriter;
	private PrintStream outVariants;
	private PrintStream outConsensus = null;
	private PrintStream clusterDetails = null;
	private boolean isPairedEnd;
	
	private KmerPrefixReadsClusteringAlgorithm parent;
	
	
	public ProcessClusterVCFTask(ReadCluster readCluster, VCFFileHeader vcfFileHeader, VCFFileWriter writer, KmerPrefixReadsClusteringAlgorithm parent, PrintStream outVariants, PrintStream outConsensus) {
		this.readCluster = readCluster;
		this.vcfFileHeader = vcfFileHeader;
		this.vcfWriter = writer;
		this.outVariants = outVariants;
		this.outConsensus = outConsensus;
		this.parent = parent;
	}
	
	public ProcessClusterVCFTask(ReadCluster readCluster, VCFFileHeader vcfFileHeader, VCFFileWriter writer, KmerPrefixReadsClusteringAlgorithm parent, PrintStream outVariants, PrintStream outConsensus, PrintStream clusterDetails) {
		this.readCluster = readCluster;
		this.vcfFileHeader = vcfFileHeader;
		this.vcfWriter = writer;
		this.outVariants = outVariants;
		this.outConsensus = outConsensus;
		this.clusterDetails = clusterDetails;
		this.parent = parent;
	}
	
	public boolean isPairedEnd() {
		return isPairedEnd;
	}

	public void setPairedEnd(boolean isPairedEnd) {
		this.isPairedEnd = isPairedEnd;
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
		
		if (outConsensus != null) {
			synchronized (outConsensus) {
				writeConsensusFasta();
				if(clusterDetails != null) {
					writeClusterDetails();
				}
			}
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
	
	private void writeClusterDetails() {
		clusterDetails.println(readCluster.getClusterNumber()+"\t"+readCluster.getNumberOfTotalReads());
	}

	private List<VCFRecord> generateRecordsForCluster() {
		List<VCFRecord> records = new ArrayList<>();
		List<ReadAlignment> readAlignments = new ArrayList<>();
		int clusterId = readCluster.getClusterNumber();
		readCluster.buildAlignment();
		String consensus = readCluster.getConsensusSequence().toUpperCase();
		int consensusLength = consensus.length();
		if(consensusLength<40) return records;
		String referenceId = Integer.toString(clusterId);
		QualifiedSequence refQS = new QualifiedSequence(referenceId, consensus);
		ReferenceGenome singleSequenceGenome = new ReferenceGenome (refQS);
		VariantPileupListener variantsDetector = new VariantPileupListener();
		variantsDetector.setGenome(singleSequenceGenome);
		
		// For each read within the cluster create a ReadAlignment. Set characters and quality scores
		List<RawRead> alignedReads = readCluster.getAlignedReads();
		List<String> sampleIds = readCluster.getSampleIds();
		if (alignedReads.size()!=sampleIds.size()) {
			System.err.println("Inconsistent number of reads and samples. reads: "+alignedReads.size()+" samples: "+sampleIds);
			return records;
		}
		
		for(int i=0;i<alignedReads.size();i++) {
			RawRead read = alignedReads.get(i);
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
		
		if(readCluster.getBreakPosition() != null) {
			consensusLength = readCluster.getBreakPosition();
		} 
		
		for(int i=1; i<=consensusLength; i++) {
			
			PileupRecord clusterPileUp = new PileupRecord(referenceId, i);
			for(ReadAlignment readAlgn:readAlignments) {
				clusterPileUp.addAlignment(readAlgn);
			}
			double h = parent.getHeterozygosityRate();
			List<Sample> samples = parent.getSamples();
			GenomicVariant variant = findMultiallelicVariant(samples, clusterPileUp, consensus.charAt(i-1), referenceId, h);
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
		records = filterRecords(records);
		return records;
	}

	private void writeConsensusFasta() {
		outConsensus.println(">Cluster_" + readCluster.getClusterNumber());
		if(readCluster.getBreakPosition() != null) {
			outConsensus.println(readCluster.getConsensusSequence().substring(0, readCluster.getBreakPosition()));
		} else {
			outConsensus.println(readCluster.getConsensusSequence());
		}
		
	}
	
	
	// From MultisampleVariantsDetector.java
	private GenomicVariant findMultiallelicVariant(List<Sample> samples, PileupRecord clusterPileUp, char reference, String clusterNum, double h) {
		
		CountsHelper helperSNV = VariantDiscoverySNVQAlgorithm.calculateCountsSNV(clusterPileUp,parent.getMaxBaseQS(), null);
		//posPrint used to be here
		GenomicVariant variant;
		variant = findMultiallelicSNV(samples, clusterPileUp, helperSNV, reference, h);
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
			allelesSupported[i]=counts[i]>0 && (double)counts[i]/(double)sum >=parent.getMinAlleleDepthFrequency();
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
	private List<VCFRecord> filterRecords(List<VCFRecord> records) {
		boolean allSNV = true;
		int n = records.size();
		for(VCFRecord record: records) {
			if (!(record.getVariant() instanceof SNV)) {
				allSNV = false;
				break;
			}
		}
		if (n<3 && allSNV) return records;
		List<VCFRecord> filtered = new ArrayList<VCFRecord>();
		for(int i=0;i<n;i++) {
			VCFRecord record = records.get(i);
			int pos = record.getFirst();
			int numClose = 0;
			for(int j=i+1;j<n && numClose<=1;j++) {
				if(records.get(j).getFirst()-pos < 10) numClose++;
			}
			if(numClose>1) break;
			filtered.add(record);
		}
		return filtered;
	}
}
