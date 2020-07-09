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
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.logging.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import ngsep.alignments.ReadAlignment;
import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class ReadAlignmentFileReader implements Iterable<ReadAlignment>,Closeable {
	
	public static final int LOAD_MODE_MINIMAL = 1;
	public static final int LOAD_MODE_ALIGNMENT = 2;
	public static final int LOAD_MODE_SEQUENCE = 3;
	public static final int LOAD_MODE_FULL = 4;
	
	
	private static final String ATTRIBUTE_NUMALNS="NH";
	
	private Logger log = Logger.getLogger(ReadAlignmentFileReader.class.getName());
	
	private SAMRecordIterator it; 
	
	private ReadAlignmentFileIterator currentIterator = null;
	private QualifiedSequenceList sequences = new QualifiedSequenceList();
	//Memory saver to avoid loading the read group for each alignment
	private QualifiedSequenceList readGroupIds = new QualifiedSequenceList();
	private Map<String,String> sampleIdsByReadGroup = new HashMap<>();
	
	private int requiredFlags = 0;
	private int filterFlags = 0;
	//TODO: Change to false by default
	private int loadMode = LOAD_MODE_FULL;
	private boolean validateHeader = true;
	private int minMQ = ReadAlignment.DEF_MIN_MQ_UNIQUE_ALIGNMENT;
	
	public ReadAlignmentFileReader (String filename) throws IOException {
		init(null,new File(filename));
	}
	public ReadAlignmentFileReader (String filename, ReferenceGenome genome, boolean validateHeader) throws IOException {
		sequences = genome.getSequencesMetadata();
		this.validateHeader = validateHeader;
		init(null,new File(filename));
	}
	public ReadAlignmentFileReader (String filename, QualifiedSequenceList sequences, boolean validateHeader) throws IOException {
		this.sequences = sequences;
		this.validateHeader = validateHeader;
		init(null,new File(filename));
	}
	//TODO: Make more constructors
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		if (log == null) throw new NullPointerException("Log can not be null");
		this.log = log;
	}
	public int getRequiredFlags() {
		return requiredFlags;
	}
	public void setRequiredFlags(int requiredFlags) {
		this.requiredFlags = requiredFlags;
	}
	public int getFilterFlags() {
		return filterFlags;
	}
	public void setFilterFlags(int filterFlags) {
		this.filterFlags = filterFlags;
	}
	public QualifiedSequenceList getSequences() {
		return sequences;
	}
	
	public List<String> getReadGroups() {
		return readGroupIds.getNamesStringList();
	}
	
	/**
	 * @return the sampleIdsByReadGroup
	 */
	public Map<String, String> getSampleIdsByReadGroup() {
		return Collections.unmodifiableMap(sampleIdsByReadGroup);
	}
	public int getLoadMode() {
		return loadMode;
	}
	public void setLoadMode(int loadMode) {
		this.loadMode = loadMode;
	}
	public int getMinMQ() {
		return minMQ;
	}
	public void setMinMQ(int minMQ) {
		this.minMQ = minMQ;
	}
	@Override
	public void close() throws IOException {
		it.close();
	}

	@Override
	public Iterator<ReadAlignment> iterator() {
		if (it == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (currentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        currentIterator = new ReadAlignmentFileIterator(); 
		return currentIterator;
	}
	
	private void init (InputStream stream, File file) throws IOException {
		if (stream != null && file != null) throw new IllegalArgumentException("Stream and file are mutually exclusive");
		if(file!=null) {
			stream = new FileInputStream(file);
		}
		//TODO: Use stream
		SamReader reader = SamReaderFactory.makeDefault().open(file);
		SAMFileHeader header = reader.getFileHeader();
		if(header != null) loadHeader(header);
		it = reader.iterator();
	}
	private void loadHeader(SAMFileHeader header) throws IOException {
		for(SAMReadGroupRecord rgRecord:header.getReadGroups()) {
			String id = rgRecord.getId();
			String sampleId = rgRecord.getSample();
			if(sampleId ==null) sampleId = id;
			if(!sampleIdsByReadGroup.containsKey(id)) {
				readGroupIds.addOrLookupName(id);
				sampleIdsByReadGroup.put(id, sampleId);
			}
		}
		if(!validateHeader) return;
		SAMSequenceDictionary dict = header.getSequenceDictionary();
		if(dict == null) return;
		for(SAMSequenceRecord seqRecord:dict.getSequences()) {
			QualifiedSequence seq = new QualifiedSequence(seqRecord.getSequenceName());
			seq.setLength(seqRecord.getSequenceLength());
			boolean added;
			try {
				added = sequences.add(seq);
			} catch (RuntimeException e) {
				throw new IOException("Inconsistent file header. Sequence "+seqRecord.getSequenceName()+" not present in the reference sequences");
			}
			if(!added) {
				QualifiedSequence seq2 = sequences.get(seq.getName());
				if(seq2!=seq && seq2.getLength()>0 && seq2.getLength()!=seq.getLength()) throw new IOException("Inconsistent length in file header. Sequence "+seqRecord.getSequenceName()+" has length: "+seq2.getLength()+" in the refrence sequences and length: "+seq.getLength()+" in the file header");
			}
		}
	}
	private int currentSequenceIdx = -1;
	private String currentSequenceName = null;
	
	private ReadAlignment loadAlignment(SAMRecord alnRecord) {
		String refName = alnRecord.getReferenceName();
		String sequenceName;
		if(currentSequenceIdx == alnRecord.getReferenceIndex()) {
			sequenceName = currentSequenceName;
		} else {
			sequenceName = loadSequenceName(refName,sequences);
			currentSequenceIdx =  alnRecord.getReferenceIndex();
			currentSequenceName = sequenceName;
		}
		int flags = alnRecord.getFlags();
		if(isMultiple(alnRecord)) flags += ReadAlignment.FLAG_MULTIPLE_ALN;
		String mateSequenceName = sequenceName;
		boolean differentSequence = alnRecord.getReferenceIndex()!= alnRecord.getMateReferenceIndex();
		if(differentSequence) {
			mateSequenceName = loadSequenceName(alnRecord.getMateReferenceName(),mateSeqNamesList);
			flags += ReadAlignment.FLAG_MATE_DIFFERENT_SEQUENCE;
		}
		
		
		ReadAlignment answer = new ReadAlignment(sequenceName, alnRecord.getAlignmentStart(), alnRecord.getAlignmentEnd(), alnRecord.getReadLength(), flags);
		answer.setSequenceIndex(currentSequenceIdx);
		if(loadMode != LOAD_MODE_MINIMAL) {
			String cigar = null;
			if(!SAMRecord.NO_ALIGNMENT_CIGAR.equals(alnRecord.getCigarString())) cigar = alnRecord.getCigarString();
			answer.setCigarString(cigar);
		}
		int alnQuality = alnRecord.getMappingQuality();
		if(alnQuality>127) alnQuality = 127;
		answer.setAlignmentQuality((byte) alnQuality);
		answer.setMateSequenceName(mateSequenceName);
		answer.setMateFirst(alnRecord.getMateAlignmentStart());
		answer.setInferredInsertSize(alnRecord.getInferredInsertSize());
		SAMReadGroupRecord readGroup = alnRecord.getReadGroup();
		if(readGroup!=null) {
			String id = readGroupIds.addOrLookupName(readGroup.getId()).getName();
			answer.setReadGroup(id);
		}
		Integer mismatches = alnRecord.getIntegerAttribute(SAMTag.NM.toString());
		if (mismatches!=null) answer.setNumMismatches((short) Math.min(mismatches,32000));
		if (loadMode >= LOAD_MODE_SEQUENCE) {
			answer.setReadCharacters(alnRecord.getReadString());
			String qs = alnRecord.getBaseQualityString();
			if(qs!=null && !SAMRecord.NULL_QUALS_STRING.equals(qs)) {
				answer.setQualityScores(qs);
			}
		}
		if(loadMode == LOAD_MODE_FULL) answer.setReadName(alnRecord.getReadName());
		return answer;
	}
	private QualifiedSequenceList mateSeqNamesList = new QualifiedSequenceList();
	private static String loadSequenceName(String refName, QualifiedSequenceList namesList) {
		String sequenceName;
		if(refName == null || refName.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
			sequenceName = null;
		} else {
			sequenceName = namesList.addOrLookupName(refName).getName();
		}
		return sequenceName;
	}
	
	private boolean isMultiple(SAMRecord aln) {
		if (aln.getNotPrimaryAlignmentFlag()) return true;
		Integer alns = aln.getIntegerAttribute(ATTRIBUTE_NUMALNS);
		if(alns!=null && alns > 1) return true;
		if(alns!=null && alns == 1) return false;
		if (aln.getMappingQuality()<minMQ) return true;
		return false;
	}
	private static boolean isSameAlignment (SAMRecord aln1, SAMRecord aln2) {
		if(aln1 == null) return aln2==null;
		if(aln2 == null) return aln1==null;
		if(aln1.getAlignmentStart()!=aln2.getAlignmentStart()) return false;
		/*String name1 = aln1.getReadName();
		int l1 = name1!=null?name1.length():0;
		String name2 = aln2.getReadName();
		int l2 = name2!=null?name2.length():0;
		if(l1!=l2) return false;
		if(l1>0 && name1.charAt(l1-1)!=name2.charAt(l2-1)) return false;*/
		if(aln1.getReadPairedFlag() != aln2.getReadPairedFlag()) return false;
		if(aln1.getReadPairedFlag() && aln1.getFirstOfPairFlag()!=aln2.getFirstOfPairFlag()) return false;
		if(!aln1.getReadName().equals(aln2.getReadName()))return false;
		return true;
	}
	
	private class ReadAlignmentFileIterator implements Iterator<ReadAlignment> {
		private ReadAlignment nextAln;
		public ReadAlignmentFileIterator() {
			nextAln = loadAln();
		}
		@Override
		public boolean hasNext() {
			return nextAln!=null;
		}

		@Override
		public ReadAlignment next() {
			if(nextAln==null) throw new NoSuchElementException();
			ReadAlignment answer = nextAln;
			nextAln = loadAln();
			return answer;
		}
		private SAMRecord lastRecord = null;
		private ReadAlignment loadAln() {
			SAMRecord alnRecord;
			while(true) {
				if(!it.hasNext()) return null;
				try {
					alnRecord = it.next();
				} catch (Exception e) {
					//e.printStackTrace();
					log.warning(e.getMessage());
					continue;
				}
				if(isSameAlignment(lastRecord, alnRecord)) continue;
				lastRecord = alnRecord;
				
				ReadAlignment answer;
				try {
					answer = loadAlignment(alnRecord);
				} catch (RuntimeException e) {
					e.printStackTrace();
					log.warning(e.getMessage());
					continue;
				}
				int flags = answer.getFlags();
				//System.out.println("Flags: "+answer.getFlags()+" filter: "+filterFlags+". filtered: "+(flags & filterFlags)+" required: "+requiredFlags+". Reqpass: "+(flags & requiredFlags));
				if((flags & filterFlags)!=0) continue;
				if((flags & requiredFlags)!=requiredFlags) continue;
				return answer;
			} 
		}
		@Override
		public void remove() {
			throw new UnsupportedOperationException("Remove not supported by ReadAlignmentFileIterator");
		}
	}
	

}
