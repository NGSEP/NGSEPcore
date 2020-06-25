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
import java.util.List;

import ngsep.alignments.ReadAlignment;
import ngsep.discovery.MultisampleVariantsDetector;
import ngsep.discovery.PileupRecord;
import ngsep.sequences.RawRead;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;
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
		MultisampleVariantsDetector mvd = new MultisampleVariantsDetector();
		double h = parent.getHeterozygosityRate();
		mvd.setHeterozygosityRate(h);
		mvd.setMaxBaseQS(parent.getMaxBaseQS());
		mvd.setMinAlleleDepthFrequency(parent.getMinAlleleDepthFrequency());
		List<Sample> samples = parent.getSamples();
		mvd.setSamples(samples);
		for(int i=1; i<=consensusLength; i++) {
			PileupRecord clusterPileUp = new PileupRecord(referenceId, i);
			for(ReadAlignment readAlgn:readAlignments) {
				clusterPileUp.addAlignment(readAlgn);
			}
			
			GenomicVariant variant = mvd.discoverPopulationVariant(clusterPileUp, consensus.substring(i-1,i));
			if(variant!=null) {
				List<CalledGenomicVariant> calls = mvd.genotypeVariant(variant, clusterPileUp);
				if(variant.getVariantQS()==0 || variant.getVariantQS() < parent.getMinQuality()) continue;
				VCFRecord record = VCFRecord.createDefaultPopulationVCFRecord(variant, calls, vcfFileHeader);
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
