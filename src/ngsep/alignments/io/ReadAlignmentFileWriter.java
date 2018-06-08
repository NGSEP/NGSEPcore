package ngsep.alignments.io;

import java.io.File;
import java.io.PrintStream;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import ngsep.alignments.ReadAlignment;

public class ReadAlignmentFileWriter {
	
	private SAMFileWriter writer;
	private SAMFileHeader samFileHeader;

	public ReadAlignmentFileWriter (PrintStream out)
	{
		samFileHeader = new SAMFileHeader();
//		writer= new SAMFileWriterFactory().makeBAMWriter(samFileHeader, true, out);
		writer= new SAMFileWriterFactory().makeBAMWriter(samFileHeader, false, out);


	}
	
	public void write(ReadAlignment readAlignment)
	{
		SAMRecord samRecord= new SAMRecord(samFileHeader);
		
		//QNAME
		samRecord.setReadName(readAlignment.getReadName());
		
		//FLAG
		samRecord.setFlags(readAlignment.getFlags());
		
		//RNAME
		samRecord.setReferenceName(readAlignment.getSequenceName());
		
		//POS
		samRecord.setAlignmentStart(readAlignment.getFirst());
		
		//MAPQ
		samRecord.setMappingQuality(readAlignment.getAlignmentQuality());
		
		//CIGAR
		samRecord.setCigarString(readAlignment.getCigarString());

		//RNEXT
		samRecord.setMateReferenceName(readAlignment.getMateSequenceName());
		
		//PNEXT
		samRecord.setMateAlignmentStart(readAlignment.getMateFirst());
		
		//TLEN
		samRecord.setInferredInsertSize(readAlignment.getInferredInsertSize());
		
		//SEQ
		String basesString = readAlignment.getReadCharacters().toString();
		samRecord.setReadBases(basesString.getBytes());
		
		//QUAL
		samRecord.setBaseQualityString(readAlignment.getQualityScores());
		
		writer.addAlignment(samRecord);
	}
}
