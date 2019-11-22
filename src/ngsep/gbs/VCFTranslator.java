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
		return genome; 
	}

	public void run(String refGenomeFile) throws Exception {
		createComplementMap();
		refGenome = loadReferenceGenome(refGenomeFile);
		getAlignmentHash();
		System.out.println("Loaded a total of " + alignmentsHash.size() + " alignments.");
//		System.out.println(refGenome.getReference(0, 16320346, 16320346));
		
		translate();
		//explore();
		}
	
	/**
	 *  For each record on the RelativeVCF, this method translates and rewrites. 
	 * @throws Exception
	 */
	private void translate() throws Exception {
		try (ReadAlignmentFileReader bamOpenFile = new ReadAlignmentFileReader(filenameAlignmentBAM);
				VCFFileReader vcfOpenFile = new VCFFileReader(filenameRelativeVCF);
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
				ReadAlignment alignments =  alignmentsHash.get(clusterID);
				if(alignments != null ) {
					// If variant is SNP
					if(record.getVariant().isSNV()) {
						translatedRecord = translateRecord(alignments, record);
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
					"Not SNV: " + notSNV);
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
		// Reference base (the base at position of variant in the consensus seq)
		int relativePos = record.getFirst();
		char relativeRef = algn.getReadCharacters().charAt(relativePos - 1);
		
		//-1 for undecided, 0 for homozygous reference, 1 for heterozygous, 2 for homozygous variant
		byte genotype = -1;
		
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
		
		// IF the reference is not in the relativeAlleles, the variant is likely triallelic
		boolean refInRelativeAllels = false;
		
		// Will be something like Cluster-####
		String seqName = algn.getSequenceName();
		
		if(algn.isReadUnmapped()) {
			unmappedRead++;
			return null;
		}
		
		//-1 if the read position does not align with the reference
		truePos = algn.getReferencePosition(relativePos);
		if(truePos == -1) {
			TruePosNotAlign++;
			return null;
		}
		
		// get the calls for this record.
		List<CalledGenomicVariant> calls = record.getCalls();
//		System.out.println(truePos);
		// TODO: Check index. Not clear if Im using it correctly
		CharSequence trueRefSeq = refGenome.getReference(0, truePos, truePos);
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
			if(allel == Character.toString(trueRef)) {
				refInRelativeAllels = true;
				continue;
			} else {
				trueAlt = allel.charAt(0);
			}
		}
		if(!refInRelativeAllels) {
			isTriallelic = true;
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
		try (VCFFileReader vcfOpenFile = new VCFFileReader(filenameRelativeVCF);) {
			Iterator<VCFRecord> VCFReader = vcfOpenFile.iterator();
			while(VCFReader.hasNext()) {
				VCFRecord algn = VCFReader.next();
				//System.out.print("Read Name: " + algn.getReadName());
				System.out.print("Sequence Name: " + algn.getSequenceName() + "\n");
			}
		}
	}
	
}
