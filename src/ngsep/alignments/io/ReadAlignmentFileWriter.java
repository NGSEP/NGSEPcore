package ngsep.alignments.io;

import java.io.File;
import java.io.PrintStream;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import ngsep.alignments.ReadAlignment;

public class ReadAlignmentFileWriter {
	
	private SAMFileWriter writer;
	private SAMFileHeader samFileHeader;

	public ReadAlignmentFileWriter (PrintStream out)
	{
		samFileHeader = new SAMFileHeader();
		SAMSequenceDictionary sequenceDictionary = new SAMSequenceDictionary();
		
		SAMSequenceRecord sequenceRecord = new SAMSequenceRecord("chrI", 100000);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrII", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrIII", 50);
		sequenceDictionary.addSequence(sequenceRecord);

		
		sequenceRecord = new SAMSequenceRecord("chrVI", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrIV", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrV", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		
		sequenceRecord = new SAMSequenceRecord("chrIX", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrXII", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrVIII", 50);
		sequenceDictionary.addSequence(sequenceRecord);
//		
		sequenceRecord = new SAMSequenceRecord("chrXVI", 50);
		sequenceDictionary.addSequence(sequenceRecord);
//		
//		
		sequenceRecord = new SAMSequenceRecord("chrXIII", 50);
		sequenceDictionary.addSequence(sequenceRecord);
//		
		sequenceRecord = new SAMSequenceRecord("chrXV", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrX", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrXIV", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrVII", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrXI", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		sequenceRecord = new SAMSequenceRecord("chrM", 50);
		sequenceDictionary.addSequence(sequenceRecord);
		
		
		
		samFileHeader.setSequenceDictionary(sequenceDictionary);
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
		// error because alignment is null 
		//samRecord.setCigarString(readAlignment.getCigarString());

		//RNEXT
		String mateReferenceName = readAlignment.getMateSequenceName()!=null?readAlignment.getMateSequenceName():SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
		samRecord.setMateReferenceName(mateReferenceName);

		
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
