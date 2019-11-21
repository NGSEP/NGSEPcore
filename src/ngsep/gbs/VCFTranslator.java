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
	private Map<String, ReadAlignment> alignmentsHash;
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
				VCFRecord translatedRecord = null;
				VCFRecord record = vcfReader.next();
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
		}
	}
	
	private VCFRecord translateRecord(ReadAlignment algn, VCFRecord record) {
		
		VCFRecord translatedRecord = null;
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
		String[] relativeAlleles = relativeVar.getAlleles();
		char relativeAlt = relativeAlleles[1].charAt(0);
		VCFFileHeader header = record.getHeader();
		int[] filedsFormat = record.getFieldsFormat();
		char trueAlt = 'N';
		boolean isTriallelic = false;
		boolean refInRelativeAllels = false;
	
		
		int read_length = algn.getReadCharacters().length();
		String seqName = algn.getSequenceName();
		truePos = algn.getReferencePosition(relativePos);
		CharSequence alleleCall = algn.getAlleleCall(truePos);
		
		// if alleleCall is null, then refChar == readChar at truePos
		// if there is variant, we must get char from the ref-genome. 
		if(alleleCall == null) {
			relativeRef = algn.getReadCharacters().charAt(relativePos);
		} else {
			//TODO: get actual relativeRef
			relativeRef = algn.getReadCharacters().charAt(relativePos);
		}
		
		if(algn.isNegativeStrand()) {
			trueRef = complement.get(relativeRef);
		} else {
			trueRef = relativeRef;
		}
		// is biallelic in relativeVCF
		if(relativeVar.isBiallelic()) {
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
		}
		// Exclude reference from allels
		if(seqName != null && (!isTriallelic)) {
			SNV variant = new SNV(seqName, truePos, trueRef, trueAlt);
			CalledSNV call = new CalledSNV(variant, genotype);
			translatedRecord = new VCFRecord(variant, filedsFormat, call, header);
		} else {
			nullReadName++;
		}
		return translatedRecord;
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
				String algnName = algn.getReadName().split("_")[1];
				if(!algn.isSecondary()) {
					this.alignmentsHash.put(algnName, algn);
				}
			}
		}
		
	}
	
}
