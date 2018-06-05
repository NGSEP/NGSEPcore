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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegion;
import ngsep.genome.ReferenceGenome;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledGenomicVariantImpl;
import ngsep.variants.GenomicVariant;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class MultisampleVariantsDetector implements PileupListener {

	private Logger log = Logger.getLogger(MultisampleVariantsDetector.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	private double coveredGenomeSize = 0;
	private long referenceGenomeSize = 0;
	
	private AlignmentsPileupGenerator generator = new AlignmentsPileupGenerator();
	//Listeners
	private IndelRealignerPileupListener indelRealigner = new IndelRealignerPileupListener();

	
	private ReferenceGenome genome;
	private List<String> alignmentFiles = new ArrayList<>();
	private List<Sample> samples = new ArrayList<>();
	private List<VariantPileupListener> genotypingListeners = new ArrayList<>();
	
	//Output file variables
	private String outFilename;
	private PrintStream outFile;
	private VCFFileHeader vcfFileHeader;
	private VCFFileWriter writer = new VCFFileWriter();
	
	private String knownSTRsFile = null;
	private short maxBaseQS=0; 
	private boolean ignoreLowerCaseRef = false;
	private boolean callEmbeddedSNVs = false;
	
	//Control attribute to avoid calling overlapping indels and to give an embedded status to SNVs within indels or STRs
	private int lastIndelEnd = 0;
	
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
	
	public void findVariants() throws IOException {
		//TODO: assign sample ids if not in aln files
		if(samples == null) loadSamplesFromAlignmentHeaders();
		referenceGenomeSize = genome.getTotalLength();
		if(knownSTRsFile!=null) {
			log.info("Loading input short tandem repeats");
			//TODO: Choose the best format
			SimpleGenomicRegionFileHandler rfh = new SimpleGenomicRegionFileHandler();
			List<GenomicRegion> strs = rfh.loadRegions(knownSTRsFile);
			//TODO: STRs loader
			//indelRealigner.setInputVariants(makeNonRedundantSTRs(strs));
			log.info("Loaded "+strs.size()+" input short tandem repeats");
		}
		log.info("Finding variants");
		indelRealigner.setGenome(genome);
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
		for(String filename:alignmentFiles) {
			try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(filename)) {
				Map<String, String> samplesHeader = reader.getSampleIdsByReadGroup();
				for(String readGroup : samplesHeader.keySet()) {
					String sampleId = samplesHeader.get(readGroup);
					Sample sample = samplesMap.get(sampleId);
					if(sample==null) {
						sample = new Sample(sampleId);
						samplesMap.put(sampleId, sample);
					}
					sample.addReadGroup(readGroup);
				}
			}
		}
	}

	@Override
	public void onPileup(PileupRecord pileup) {
		GenomicVariant variant = findVariant(pileup);
		if(variant == null) return;
		List<CalledGenomicVariant> calls = new ArrayList<>();
		int n = samples.size();
		for(int i=0;i<n;i++) {
			Sample sample = samples.get(i);
			VariantPileupListener genotyper = genotypingListeners.get(i);
			CalledGenomicVariant call = genotyper.processPileup(pileup, variant);
			if(call == null) call = new CalledGenomicVariantImpl(variant, -1);
			call.setSampleId(sample.getId());
			calls.add(call);
		}
		VCFRecord record = new VCFRecord(variant, VCFRecord.DEF_FORMAT_ARRAY_NGSEP_SNV, calls, vcfFileHeader);
		writer.printVCFRecord(record, outFile);
		coveredGenomeSize++;
		if(progressNotifier!=null && coveredGenomeSize%10000==0) {
			int progress = (int)Math.round(100.0*coveredGenomeSize/referenceGenomeSize);
			generator.setKeepRunning(progressNotifier.keepRunning(progress));
		}
	}

	private GenomicVariant findVariant(PileupRecord pileup) {
		if(!callEmbeddedSNVs && pileup.isEmbedded()) return null;
		int last = pileup.getPosition()+pileup.getReferenceSpan()-1;
		CharSequence seq = genome.getReference(pileup.getSequenceName(), pileup.getPosition(), last);
		if(seq == null) return null;
		String referenceAllele = seq.toString();
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
	
	private GenomicVariant callMultisampleSNV(PileupRecord pileup, CountsHelper helperSNV, char reference) {
		// TODO Auto-generated method stub
		return null;
	}

	private GenomicVariant callMultisampleIndel(PileupRecord pileup, CountsHelper helperIndel) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void onSequenceStart(QualifiedSequence sequence) {
		lastIndelEnd = 0;
	}

	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {
		
	}

}
