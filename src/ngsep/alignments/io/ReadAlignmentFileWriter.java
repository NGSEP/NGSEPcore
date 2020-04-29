/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.alignments.io;

import java.io.Closeable;
import java.io.PrintStream;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMValidationError;
import ngsep.alignments.ReadAlignment;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

/**
 * @author German Andrade
 * @author Jorge Duitama 
 *
 */
public class ReadAlignmentFileWriter implements Closeable {
	
	private SAMFileWriter writer=null;
	private SAMFileHeader samFileHeader;
	private String sampleId;
	private PrintStream out;

	public ReadAlignmentFileWriter (QualifiedSequenceList sequences, PrintStream out)
	{
		samFileHeader = new SAMFileHeader();
		SAMSequenceDictionary sequenceDictionary = new SAMSequenceDictionary();
		for(QualifiedSequence seq:sequences) {
			SAMSequenceRecord sequenceRecord = new SAMSequenceRecord(seq.getName(), seq.getLength());
			sequenceDictionary.addSequence(sequenceRecord);
		}
		samFileHeader.setSequenceDictionary(sequenceDictionary);
		this.out = out;
	}
	
	public void write(ReadAlignment readAlignment)
	{
		if(writer == null) writer= new SAMFileWriterFactory().makeBAMWriter(samFileHeader, false, out);
		SAMRecord samRecord= new SAMRecord(samFileHeader);
		samRecord.setReadName(readAlignment.getReadName());
		samRecord.setFlags(readAlignment.getFlags());
		
		if(readAlignment.getSequenceName()!=null) {
			samRecord.setReferenceName(readAlignment.getSequenceName());
			samRecord.setAlignmentStart(readAlignment.getFirst());
			samRecord.setCigarString(readAlignment.getCigarString());
		} else {
			samRecord.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
			samRecord.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
			samRecord.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
		}
		samRecord.setMappingQuality(readAlignment.getAlignmentQuality());
		
		if(readAlignment.getMateSequenceName()!=null) {
			samRecord.setMateReferenceName(readAlignment.getMateSequenceName());
			samRecord.setMateAlignmentStart(readAlignment.getMateFirst());
			samRecord.setInferredInsertSize(readAlignment.getInferredInsertSize());
		} else {
			samRecord.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
			samRecord.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
		}
		
		//SEQ
		samRecord.setReadString(readAlignment.getReadCharacters().toString());
		
		//QUAL
		String qs = readAlignment.getQualityScores();
		if(qs!=null) samRecord.setBaseQualityString(qs);
		else samRecord.setBaseQualityString(SAMRecord.NULL_QUALS_STRING);
		
		//Read group
		samRecord.setAttribute(SAMTag.RG.toString(), sampleId);
		samRecord.setAttribute(SAMTag.NM.toString(), Integer.valueOf(readAlignment.getNumMismatches()));
		
		//System.out.println("Bases: "+samRecord.getReadString()+" qual: "+samRecord.getBaseQualityString());
		List<SAMValidationError> errors= samRecord.isValid();
		if(errors!=null) System.out.println("errors: "+errors.size()+errors.get(0));
		writer.addAlignment(samRecord);
		//System.out.println("Saved record ");	
	}
	public void close() {
		writer.close();
		out.close();
	}
	public void setSampleInfo(String sampleId, ReadAlignment.Platform platform) {
		this.sampleId = sampleId;
		SAMReadGroupRecord sampleRecord = new SAMReadGroupRecord(sampleId);
		sampleRecord.setSample(sampleId);
		sampleRecord.setPlatform(platform.toString());
		samFileHeader.addReadGroup(sampleRecord);
	}
}
