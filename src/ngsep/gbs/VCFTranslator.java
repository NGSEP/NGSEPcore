package ngsep.gbs;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.main.CommandsDescriptor;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;
import ngsep.variants.SNV;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileReader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class VCFTranslator {
	private String filenameConsensusFasta;
	private String outPrefix="./output";
	private String filenameAlignmentBAM;
	private String filenameRelativeVCF;
	private Map<String, ArrayList<ReadAlignment>> alignmentsHash;
	int numTranslatedRecords = 0;
	int totalRecords = 0;
	int skippedRecords = 0;
	int nullReadName = 0;
	private static Map<Character, Character> complement = new HashMap<>();

	public static void main(String[] args) throws Exception {
		VCFTranslator instance = new VCFTranslator();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.filenameAlignmentBAM = args[i++];
		instance.filenameRelativeVCF = args[i++];
		instance.outPrefix = args[i++];
		createComplementMap();
		instance.run();
	}

	private static void createComplementMap() {
		complement.put('A', 'T');
		complement.put('T', 'A');
		complement.put('C', 'G');
		complement.put('G', 'C');

	}

	public void run() throws Exception {
		getAlignmentHash();
		System.out.println("Loaded a total of " + alignmentsHash.size() + " records to be translated.");
		translate();
		System.out.println("Translated a total of " + numTranslatedRecords + " alignments.");
		System.out.println("Skipped a total of " + skippedRecords + ".");
		System.out.println("Total null ReadNames " + nullReadName);
		//explore();
		}
	private void translate() throws Exception {
		try (ReadAlignmentFileReader bamOpenFile = new ReadAlignmentFileReader(filenameAlignmentBAM);
				VCFFileReader vcfOpenFile = new VCFFileReader(filenameRelativeVCF);
				PrintStream mappedVCF = new PrintStream(outPrefix+"_mapped.vcf");
				PrintStream log = new PrintStream(outPrefix + "_translation.log");) {
			boolean writeHeader = true;
			VCFFileWriter writer = new VCFFileWriter ();
			Iterator<VCFRecord> vcfReader = vcfOpenFile.iterator();
			while(vcfReader.hasNext()) {
				ArrayList<VCFRecord> translatedRecords = new ArrayList<>();
				VCFRecord record = vcfReader.next();
				if(writeHeader) {
					writeHeader = false;
					VCFFileHeader header = record.getHeader();
					writer.printHeader(header, mappedVCF);
				}
				String sequenceName = record.getSequenceName();
				ArrayList<ReadAlignment> alignments =  alignmentsHash.get(sequenceName);
				if(alignments != null ) {
					// If variant is SNP
					if(record.getVariant().isSNV()) {
						translatedRecords = translateSNP(alignments, record);
					}
					writer.printVCFRecords(translatedRecords, mappedVCF);
					numTranslatedRecords++;
				} else {
					skippedRecords++;
				}
				totalRecords++;
			}
		}
	}
	
	private ArrayList<VCFRecord> translateSNP(ArrayList<ReadAlignment> alignments, VCFRecord record) {
		ArrayList<VCFRecord> mappedRecords = new ArrayList<>();

		// Position of alignment with respect to reference forward strand (1-based)
		int refPos;
		// True position of variant with respect to reference forward strand (1-based)
		int truePos;
		// Reference nucleotide as interpreted by the alignment to the reference
		char trueRef;
		// Position of variant with respect to the de-novo alignment
		int relativePos = record.getFirst();
		// Reference nucleotide as interpreted by the de-novo alignment
		char relativeRef;
		byte genotype = 0;
		// Variant as found in de-novo vcf
		GenomicVariant relativeVar = record.getVariant();
		// Alternative Allele based on de-novo VCF
		String[] relativeAllels = relativeVar.getAlleles();
		char relativeAlt = relativeAllels[1].charAt(0);
		VCFFileHeader header = record.getHeader();
		int[] filedsFormat = record.getFieldsFormat();

		boolean refInRelativeAllels = false;
		// True alternative allel, set to relativeAlt to init
		char trueAlt = relativeAlt;
		for(ReadAlignment algn: alignments) {
			int read_length = algn.getReadCharacters().length();
			String seqName = algn.getSequenceName();
			refPos = algn.getFirst();
			if(algn.isNegativeStrand()) {
				// If read aligns to reverse strand
				truePos = refPos + read_length - relativePos;
				relativeRef = algn.getReadCharacters().charAt(read_length - relativePos);
				trueRef = complement.get(relativeRef);
			} else {
				// If read aligns to forward strand
				truePos = relativePos + refPos - 1;
				relativeRef = algn.getReadCharacters().charAt(relativePos - 1);
				trueRef = relativeRef;
			}
			// Exclude reference from allels
			if(relativeVar.isBiallelic()) {
				for(String allel: relativeAllels) {
					if(allel == Character.toString(trueRef)) {
						refInRelativeAllels = true;
						continue;
					} else {
						trueAlt = allel.charAt(0);
					}
				}
//				if(!refInRelativeAllels) {
//					throw new RuntimeException("Reported reference allel not found in reported alleles");
//				}
			}
			if(seqName != null ) {
				SNV variant = new SNV(seqName, truePos, trueRef, trueAlt);
				CalledSNV call = new CalledSNV(variant, genotype);
				VCFRecord translatedRecord = new VCFRecord(variant, filedsFormat, call, header);
				mappedRecords.add(translatedRecord);
			} else {
				nullReadName++;
			}
		}
		return mappedRecords;
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

	private void getAlignmentHash() throws IOException {
		this.alignmentsHash = new HashMap<>();
		try (ReadAlignmentFileReader bamOpenFile = new ReadAlignmentFileReader(filenameAlignmentBAM);) {
			Iterator<ReadAlignment> bamReader = bamOpenFile.iterator();
			while(bamReader.hasNext()) {
				ReadAlignment algn = bamReader.next();
				if(algn.isUnique()) {
					ArrayList<ReadAlignment> algns = this.alignmentsHash.get(algn.getReadName());
					if(algns == null) {
						algns = new ArrayList<>();
					}
					algns.add(algn);
					this.alignmentsHash.put(algn.getReadName(), algns);
				}
			}
		}
		
	}
	
}
