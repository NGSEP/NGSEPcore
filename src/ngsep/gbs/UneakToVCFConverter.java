package ngsep.gbs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.SNV;
import ngsep.variants.Sample;
import ngsep.vcf.VCFFileHeader;
import ngsep.vcf.VCFFileWriter;
import ngsep.vcf.VCFRecord;

public class UneakToVCFConverter {

	public static void main(String[] args) throws Exception {
		UneakToVCFConverter converter = new UneakToVCFConverter();
		String hapmapFile = args[0];
		String consensusFile = args[1];
		String outPrefix = args[2];
		converter.process(hapmapFile,consensusFile,outPrefix);
	}

	public void process(String hapmapFile, String consensusFile, String outPrefix) throws IOException {
		FastaSequencesHandler fastaHandler = new FastaSequencesHandler();
		QualifiedSequenceList seqs = fastaHandler.loadSequences(consensusFile);
		List<VCFRecord> records = new ArrayList<VCFRecord>();
		VCFFileHeader header = VCFFileHeader.makeDefaultEmptyHeader();
		List<Sample> samples = new ArrayList<Sample>();
		List<QualifiedSequence> consensus = new ArrayList<QualifiedSequence>(seqs.size());
		try (FileReader reader = new FileReader(hapmapFile);
			 BufferedReader in = new BufferedReader(reader)) {
			String line = in.readLine();
			String [] items = line.split("\t");
			for(int i=11;i<items.length;i++) {
				Sample s = new Sample(items[i]);
				samples.add(s);
				header.addSample(s, false);
			}
			line = in.readLine();
			for (int i=0;line!=null;i+=2) {
				items = line.split("\t");
				//System.out.println("Processing line: "+i);
				char a1 = items[1].charAt(0);
				char a2 = items[1].charAt(2);
				QualifiedSequence s1 = seqs.get(i);
				QualifiedSequence s2 = seqs.get(i+1);
				String seqName = items[0];
				int pos = calculatePosition (a1,a2,s1.getCharacters().toString(),s2.getCharacters().toString());
				consensus.add(new QualifiedSequence(seqName,s1.getCharacters()));
				SNV snv = new SNV(seqName, pos, a1, a2);
				List<CalledGenomicVariant> calls = new ArrayList<CalledGenomicVariant>(samples.size());
				for(int j =11;j<items.length;j++) {
					CalledSNV call = makeCalledSNV(snv,items[j],samples.get(j-11));
					calls.add(call);
				}
				VCFRecord record = new VCFRecord(snv, VCFRecord.DEF_FORMAT_ARRAY_MINIMAL, calls, header);
				records.add(record);
				
				line = in.readLine();
			}
		}
		try (PrintStream out = new PrintStream(outPrefix+"_consensus.fa")) {
			fastaHandler.saveSequences(consensus, out, 100);
		}
		VCFFileWriter writer = new VCFFileWriter();
		try (PrintStream out = new PrintStream(outPrefix+".vcf")) {
			writer.printHeader(header, out);
			writer.printVCFRecords(records, out);
		}
		
		
	}

	private int calculatePosition(char a1, char a2, String s1, String s2) {
		for(int i=0;i<s1.length();i++) {
			char s1c = s1.charAt(i);
			char s2c = s2.charAt(i);
			if(s1c == a1 && s2c == a2) return i+1;
		}
		return 0;
	}

	private CalledSNV makeCalledSNV(SNV snv, String genotype, Sample sample) {
		char g = genotype.charAt(0);
		CalledSNV call;
		if(g=='N') call =  new CalledSNV(snv, CalledGenomicVariant.GENOTYPE_UNDECIDED);
		else if (g==snv.getReferenceBase()) call = new CalledSNV(snv, CalledGenomicVariant.GENOTYPE_HOMOREF);
		else if (g==snv.getAlternativeBase()) call = new CalledSNV(snv, CalledGenomicVariant.GENOTYPE_HOMOALT);
		else call= new CalledSNV(snv, CalledGenomicVariant.GENOTYPE_HETERO);
		call.setSampleId(sample.getId());
		return call;
	}

}
