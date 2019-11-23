package ngsep.gbs;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.genome.GenomicRegion;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.sequences.DNASequence;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.SNV;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class VCFTranslator {

	private String outPrefix="./output";
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
	int refNotInAlleles = 0;
	int biallelic = 0;
	int triallelic = 0;
	private static Map<Character, Character> complement = new HashMap<>();

	public static void main(String[] args) throws Exception {
		VCFTranslator instance = new VCFTranslator();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.filenameAlignmentBAM = args[i++];
		instance.filenameRelativeVCF = args[i++];
		String refGenomeFile = args[i++];
		instance.outPrefix = args[i++];
		instance.run(refGenomeFile);
	}

	private static void createComplementMap() {
		complement.put('A', 'T');
		complement.put('T', 'A');
		complement.put('C', 'G');
		complement.put('G', 'C');

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
		createComplementMap();
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
		try (VCFFileReader vcfOpenFile = new VCFFileReader(filenameRelativeVCF);
				PrintStream mappedVCF = new PrintStream(outPrefix+"_mapped.vcf");
				PrintStream log = new PrintStream(outPrefix + "_translation.log");) {
			
			boolean writeHeader = true;
			VCFFileWriter writer = new VCFFileWriter ();
			Iterator<VCFRecord> vcfReader = vcfOpenFile.iterator();
			
			// Iterate over vcfRecords
			while(vcfReader.hasNext()) {
				VCFRecord translatedRecord = null;
				VCFRecord record = vcfReader.next();
				
				// Write header on translatedVCF
				if(writeHeader) {
					writeHeader = false;
					VCFFileHeader header = record.getHeader();
					writer.printHeader(header, mappedVCF);
				}
				
				String clusterID = record.getSequenceName().split("_")[1];
				ReadAlignment alignment =  alignmentsHash.get(clusterID);
				if(alignment != null ) {
					// If variant is SNP
					if(record.getVariant().isSNV()) {
						translatedRecord = translateRecord(alignment, record);
					} else {
						notSNV++;
					}
					if(translatedRecord != null) {
						writer.printVCFRecord(translatedRecord, mappedVCF);
						numTranslatedRecords++;
					}
				} else {
					skippedRecords++;
				}
				totalRecords++;
			}
			log.print("Total number of records in relative VCF: " + totalRecords + "\n" +
					"Number of translated records: " + numTranslatedRecords + "\n" +
					"Number of skipped records due to lack of a primary alignment: " + skippedRecords + "\n" +
					"Not DNA: " + notDNA + "\n" + 
					"True pos not aligned to ref: " + TruePosNotAlign +"\n"+
					"Unmapped read: " + unmappedRead + "\n" +
					"Not SNV: " + notSNV + "\n" +
					"RefAllele not found in alleles: " + refNotInAlleles + "\n" +
					"Tirallelic: " + triallelic + "\n" +
					"Biallelic: " + biallelic);
		}
	}
	
	/**
	 * This method takes in an alignment to the reference and a record 
	 * and translates each record to match the reference.
	 * @param algn
	 * @param record
	 * @return translated record
	 */
	private VCFRecord translateRecord(ReadAlignment algn, VCFRecord record) {
		
		VCFRecord translatedRecord = null;
		List<CalledGenomicVariant> trueCalls = new ArrayList<>();

		// True position of variant with respect to reference forward strand (1-based)
		int truePos;
		
		// Reference nucleotide as interpreted by the alignment to the reference
		char trueRef;
		
		// Position of variant with respect to the de-novo alignment 1-based
		// Relative reference base (the base at position of variant in the consensus seq)
		int relativePos = record.getFirst();
		char relativeRef = algn.getReadCharacters().charAt(relativePos - 1);
		int seqLength = record.length();
		
		
		//-1 for undecided, 0 for homozygous reference, 1 for heterozygous, 2 for homozygous variant
		byte genotype = -1;
		
		// describes the agreement between the reported reference in relative VCF and 
		// true reference in the reference genome
		boolean agreement = false;
		
		// Variant as found in de-novo vcf
		GenomicVariant relativeVar = record.getVariant();
		
		// Alternative Allele based on de-novo VCF
		String[] relativeAlleles = relativeVar.getAlleles();
		
		VCFFileHeader header = record.getHeader();
		int[] filedsFormat = record.getFieldsFormat();
		char trueAlt = 'N';
		
		//TODO: set up logic to deal with triallelic
		boolean isTriallelic = false;
		boolean isBiallelic = false;
		
		// IF the reference is not in the relativeAlleles, the variant is likely triallelic
		boolean refInRelativeAllels = false;
		
		// Will be something like Cluster-####
		String seqName = algn.getSequenceName();
		
		if(algn.isReadUnmapped()) {
			unmappedRead++;
			return null;
		}
		
		//-1 if the read position does not align with the reference
//		XXX: THIS DOESNT WORK BECAUSE THERE IS NO WAY TO GET REFBASE WITH ABSOLUTE POSITION. 
//		truePos = algn.getReferencePosition(relativePos);
		
		//This allows is to get the true position, not absolute, but relative to the genomic region.
		// This plus the sequence index should get us the reference Allel.
		if(algn.isNegativeStrand()) {
			truePos = algn.getFirst() + (seqLength - relativePos);
		} else {
			truePos = algn.getFirst() + relativePos - 1;
		}
		
		//Get index of alignment -> make it 0-based. Basically the index of the genomic region
		int algnIndex = algn.getSequenceIndex();
		if(algnIndex == -1) {
			System.out.println("algnIndex = -1");
		}
		
		CharSequence trueRefSeq = refGenome.getReference(algnIndex, truePos, truePos);
		
		// get the calls for this record.
		List<CalledGenomicVariant> calls = record.getCalls();

		if(!DNASequence.isDNA(trueRefSeq)) {
			notDNA++;
			return null;
		} else {
			trueRef = trueRefSeq.charAt(0);
		}
		if(trueRef == relativeRef) {
			agreement = true;
		} else {
			agreement = false;
		}
		
		if(algn.isNegativeStrand()) {
			trueRef = complement.get(trueRef);
		} 
		
		for(String allel: relativeAlleles) {
			if(allel.charAt(0) == trueRef) {
				refInRelativeAllels = true;
				continue;
			} else {
				trueAlt = allel.charAt(0);
			}
		}
		
		if(relativeAlleles.length == 2) {
			if(refInRelativeAllels) {
				biallelic++;
				isBiallelic = true;
			} else {
				triallelic++;
				isTriallelic = true;
			}
		} else if((relativeAlleles.length == 3) && refInRelativeAllels) {
			triallelic++;
			isTriallelic = true;
		} else {
			return null;
		}		
		
		if(!refInRelativeAllels) {
			refNotInAlleles++;
			
		}
		
		// Exclude reference from allels
		if(seqName != null && (!isTriallelic)) {
			SNV variant = new SNV(seqName, truePos, trueRef, trueAlt);
			
			// Translate genotypes 
			//-1 for undecided, 0 for homozygous reference, 1 for heterozygous, 2 for homozygous variant
			for(CalledGenomicVariant relativeCall: calls) {
				if(relativeCall.isHomozygousReference() && agreement) {
					genotype = 0;
				} else if(relativeCall.isHomozygousReference() && !agreement) {
					genotype = 2;
				} else if(relativeCall.isHomozygous() && agreement) {
					genotype = 2;
				} else if(relativeCall.isHomozygous() && !agreement) {
					genotype = 0;
				} else if(relativeCall.isHeterozygous()) {
					genotype = 1;
				} else {
					genotype = -1;
				}
				
				CalledGenomicVariant trueCall = new CalledSNV(variant, genotype);
				trueCalls.add(trueCall);
			}
			translatedRecord = new VCFRecord(variant, filedsFormat, trueCalls, header);
		
		} else {
			nullReadName++;
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
				String algnName = algn.getReadName().split("_")[1];
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
