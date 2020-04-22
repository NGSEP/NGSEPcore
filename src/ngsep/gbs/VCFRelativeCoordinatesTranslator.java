package ngsep.gbs;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.DNASequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class VCFRelativeCoordinatesTranslator {

	private String outFile="output.vcf";
	private String filenameAlignmentBAM;
	private String filenameRelativeVCF;
	private static ReferenceGenome refGenome;
	private Map<String, ReadAlignment> alignmentsHash;
	int numTranslatedRecords = 0;
	int totalRecords = 0;
	int skippedRecords = 0;
	int nullReadName = 0;
	int notDNA = 0;
	int TruePosNotAlign = 0;
	int unmappedRead = 0;
	int notSNV = 0;
	int notRefSeq = 0;
	int refNotInAlleles = 0;
	int biallelic = 0;
	int triallelic = 0;
	int nonVariant = 0;
	
//	public static void main(String[] args) throws Exception {
//		VCFRelativeCoordinatesTranslator instance = new VCFRelativeCoordinatesTranslator();
//		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
//		instance.filenameAlignmentBAM = args[i++];
//		instance.filenameRelativeVCF = args[i++];
//		String refGenomeFile = args[i++];
//		instance.outFile = args[i++];
//		instance.run(refGenomeFile);
//	}
	public static void main(String[] args) throws Exception {
		VCFRelativeCoordinatesTranslator instance = new VCFRelativeCoordinatesTranslator();
		instance.filenameAlignmentBAM = args[0];
		instance.filenameRelativeVCF = args[1];
		String refGenomeFile = args[2];
		instance.outFile = args[3];
		instance.run(refGenomeFile);
	}

	
	/**
	 * Load reference genome based on provided file
	 * @param genomeFileName
	 * @return
	 * @throws IOException
	 */
	private static ReferenceGenome loadReferenceGenome(String genomeFileName) throws IOException {
		ReferenceGenome genome = new ReferenceGenome(genomeFileName);
		System.out.println(genome.getTotalLength());
		System.out.println(genome.getNumSequences());
		return genome; 
	}

	public void run(String refGenomeFile) throws Exception {
		refGenome = loadReferenceGenome(refGenomeFile);
		getAlignmentHash();
		System.out.println("Loaded a total of " + alignmentsHash.size() + " alignments.");
		
		translate();
//		explore();
		}
	
	/**
	 *  For each record on the RelativeVCF, this method translates and rewrites. 
	 * @throws Exception
	 */
	private void translate() throws Exception {
		List<VCFRecord> translatedRecords = new ArrayList<>();
		List<Sample> samples;
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		
		try (VCFFileReader vcfOpenFile = new VCFFileReader(filenameRelativeVCF)) {
			samples = vcfOpenFile.getHeader().getSamples();
			header.addSamples(samples);
			Iterator<VCFRecord> vcfReader = vcfOpenFile.iterator();
			// Iterate over vcfRecords
			while(vcfReader.hasNext()) {
				VCFRecord translatedRecord = null;
				VCFRecord record = vcfReader.next();
				
				
					
				
				int clusterID = Integer.parseInt(record.getSequenceName());
				ReadAlignment alignment =  alignmentsHash.get(clusterID);
				if(alignment != null ) {
					// If variant is SNP
					if(record.getVariant().isSNV()) {
						translatedRecord = translateRecord(alignment, record, header);
					} else {
						notSNV++;
					}
					if(translatedRecord != null) {
						translatedRecords.add(translatedRecord);
						numTranslatedRecords++;
					}
				} else {
					skippedRecords++;
				}
				totalRecords++;
			}
			
		}
		Collections.sort(translatedRecords,new GenomicRegionComparator(refGenome.getSequencesMetadata()));
		VCFFileWriter writer = new VCFFileWriter ();
		try (PrintStream mappedVCF = new PrintStream(outFile)) {
			writer.printHeader(header, mappedVCF);
			writer.printVCFRecords(translatedRecords, mappedVCF);
		}
		
		
		System.out.println("Total number of records in relative VCF: " + totalRecords);
		System.out.println("Number of translated records: " + numTranslatedRecords);
		System.out.println("Number of skipped records due to lack of a primary alignment: " + skippedRecords);
		System.out.println("Not DNA: " + notDNA); 
		System.out.println("True pos not aligned to ref: " + TruePosNotAlign);
		System.out.println("Unmapped read: " + unmappedRead);
		System.out.println("Not SNV: " + notSNV);
		System.out.println("RefAllele not found in alleles: " + refNotInAlleles);
		System.out.println("Triallelic: " + triallelic);
		System.out.println("Biallelic: " + biallelic);
		System.out.println("Non variant: " + nonVariant);
		System.out.println("Reference seq is null: " + notRefSeq);
	}
	
	/**
	 * This method takes in an alignment to the reference and a record 
	 * and translates each record to match the reference.
	 * @param algn
	 * @param record
	 * @return translated record
	 */
	private VCFRecord translateRecord(ReadAlignment algn, VCFRecord record, VCFFileHeader header) {
		
		VCFRecord translatedRecord = null;
		List<CalledGenomicVariant> trueCalls = new ArrayList<>();

		// True position of variant with respect to reference forward strand (1-based)
		int truePos;
		
		// Reference nucleotide as interpreted by the alignment to the reference
		char trueRef;
		
		// Position of variant with respect to the de-novo alignment 1-based
		// Relative reference base (the base at position of variant in the consensus seq)
		int relativePos = record.getFirst();
		
		// Variant as found in de-novo vcf
		GenomicVariant relativeVar = record.getVariant();
		
		if(!(relativeVar instanceof SNV)) notSNV++;
		
		// Alternative Allele based on de-novo VCF
		String[] relativeAlleles = relativeVar.getAlleles();
		
		int[] filedsFormat = record.getFieldsFormat();
		
		// IF the reference is not in the relativeAlleles, the variant is likely triallelic
		boolean refInRelativeAllels = false;
		
		// Will be something like Cluster-####
		String seqName = algn.getSequenceName();
		
		if(algn.isReadUnmapped() || seqName==null) {
			unmappedRead++;
			return null;
		}
		
		
		//-1 if the read position does not align with the reference
//		XXX: THIS DOESNT WORK BECAUSE THERE IS NO WAY TO GET REFBASE WITH ABSOLUTE POSITION. 
//		truePos = algn.getReferencePosition(relativePos);
		
		//This allows is to get the true position, not absolute, but relative to the genomic region.
		// This plus the sequence index should get us the reference Allel.
		if(algn.isNegativeStrand()) {
			// truePos = algn.getFirst() + (seqLength - relativePos);
			truePos = algn.getLast() - (relativePos - 1);
//			System.out.println(truePos + "\t" + (algn.getFirst() + (algn.getReadLength() - relativePos)));
		} else {
			truePos = algn.getFirst() + relativePos - 1;
		}
		
		//Get index of alignment -> make it 0-based. Basically the index of the genomic region
		int algnIndex = algn.getSequenceIndex();
		if(algnIndex == -1) {
			System.out.println("algnIndex = -1");
		} 
		
		CharSequence trueRefSeq = refGenome.getReference(algnIndex, truePos, truePos);
		if(trueRefSeq == null) {
			notRefSeq++;
			return null;
		}
		
		if(!DNASequence.isDNA(trueRefSeq)) {
			notDNA++;
			return null;
		} else {
			trueRef = trueRefSeq.charAt(0);
		}
		//TODO: Translate indels
		List<String> refBasedAlleles = new ArrayList<>();
		refBasedAlleles.add(""+trueRef);
		for(String allele:relativeAlleles) {
			if(!DNASequence.isDNA(allele)) continue;
			if(algn.isNegativeStrand()) {
				allele = DNASequence.getReverseComplement(allele).toString();
			}
			if(allele.charAt(0)==trueRef) refInRelativeAllels = true;
			if(!refBasedAlleles.contains(allele)) {
				refBasedAlleles.add(allele);
			}
		}
		GenomicVariant variant;
		if(refBasedAlleles.size()== 2) {
			variant = new SNV(seqName, truePos, trueRef, refBasedAlleles.get(1).charAt(0));
			biallelic++;
		} else if (refBasedAlleles.size()>= 3) {
			variant = new GenomicVariantImpl(seqName, truePos, refBasedAlleles );
			triallelic++;
		} else {
			nonVariant++;
			return null;
		}
		
		if(!refInRelativeAllels) {
			refNotInAlleles++;
			
		}
		// get the calls for this record.
		List<CalledGenomicVariant> calls = record.getCalls();
		// Exclude reference from allels
		
			
		// Translate genotypes 
		//-1 for undecided, 0 for homozygous reference, 1 for heterozygous, 2 for homozygous variant
		for(CalledGenomicVariant relativeCall: calls) {
			String [] calledAlleles = relativeCall.getCalledAlleles();
			if (algn.isNegativeStrand()) {
				for(int i=0;i<calledAlleles.length;i++) {
					calledAlleles[i] = DNASequence.getReverseComplement(calledAlleles[i]).toString();
				}
			}
			if(variant instanceof SNV) {
				byte genotype = CalledGenomicVariant.GENOTYPE_UNDECIDED;
				if(calledAlleles.length==2) genotype = CalledGenomicVariant.GENOTYPE_HETERO;
				else if (calledAlleles.length==1) {
					genotype = CalledGenomicVariant.GENOTYPE_HOMOREF;
					if(calledAlleles[0].charAt(0)!=trueRef) {
						genotype = CalledGenomicVariant.GENOTYPE_HOMOALT;
					}
				}
				CalledGenomicVariant trueCall = new CalledSNV((SNV)variant, genotype);
				trueCalls.add(trueCall);
			}
		}
		if(trueCalls.size()>0) {
			translatedRecord = new VCFRecord(variant, filedsFormat, trueCalls, header);
		}
		return translatedRecord;
	}
	
	/**
	 * This method builds a hash with the provided alignments to be consulted at translation time
	 * @throws IOException
	 */
	private void getAlignmentHash() throws IOException {
		this.alignmentsHash = new HashMap<>();
		try (ReadAlignmentFileReader bamOpenFile = new ReadAlignmentFileReader(filenameAlignmentBAM);) {
			Iterator<ReadAlignment> bamReader = bamOpenFile.iterator();
			while(bamReader.hasNext()) {
				ReadAlignment algn = bamReader.next();
				String algnName = "";
				if(algn.getReadName().contains("_")) {
					algnName = algn.getReadName().split("_")[1];
				} else {
					algnName = algn.getReadName();
				}
				//int clusterId = Integer.parseInt(algnName);
				if(!algn.isSecondary()) {
					this.alignmentsHash.put(algnName, algn);
				}
			}
		}
		
	}
	
	
	private void explore() throws IOException {
		try (ReadAlignmentFileReader bamOpenFile = new ReadAlignmentFileReader(filenameAlignmentBAM);) {
			Iterator<ReadAlignment> bamReader = bamOpenFile.iterator();
			int maxIndex = 0;
			int minIndex = 100;
			while(bamReader.hasNext()) {
				ReadAlignment algn = bamReader.next();
//				System.out.println("Read Name: " + algn.getReadName());
				if(algn.getSequenceIndex() >= maxIndex) {
					maxIndex = algn.getSequenceIndex();
				}
				if(algn.getSequenceIndex() <= minIndex) {
					minIndex = algn.getSequenceIndex();
				}
//				System.out.println("Sequence index: " + algn.getSequenceIndex());
			}
			System.out.println("Max index: " + maxIndex);
			System.out.println("Min index: " + minIndex);
		}
	}
	
}
