package ngsep.vcf;

import java.io.PrintStream;
import java.util.Iterator;

public class VCFGenerateADField {

	public static void main(String[] args) throws Exception {
		VCFFileWriter writer = new VCFFileWriter();
		int [] format = {VCFRecord.FORMAT_IDX_GT,VCFRecord.FORMAT_IDX_PL,VCFRecord.FORMAT_IDX_GQ,VCFRecord.FORMAT_IDX_DP,VCFRecord.FORMAT_IDX_ADP_GATK,VCFRecord.FORMAT_IDX_ACN};
		PrintStream out = System.out;
		try (VCFFileReader reader = new VCFFileReader(args[0])) {
			VCFFileHeader header = reader.getHeader();
			VCFHeaderLine hl = new VCFHeaderLine("FORMAT", VCFRecord.FORMAT_ALLELE_DEPTH_GATK,"\"Counts for observed alleles, including the reference allele\"","R","Integer");
			header.addHeaderLine(hl);
			writer.printHeader(header, out);
			Iterator<VCFRecord> it = reader.iterator();
			while (it.hasNext()) {
				VCFRecord record = it.next();
				record.setFieldsFormat(format);
				writer.printVCFRecord(record, out);
			}
		}

	}

}
