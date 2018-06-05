package ngsep.alignments.io;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import ngsep.alignments.ReadAlignment;

public class ReadAlignmentFileWriter {
	
	public void write(ReadAlignment readAlignment)
	{
		SAMFileHeader samFileHeader = new SAMFileHeader();
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
		//samRecord.setReadBases(readAlignment.getReadCharacters());
		
		//QUAL
		samRecord.setBaseQualityString(readAlignment.getQualityScores());

	}
}
