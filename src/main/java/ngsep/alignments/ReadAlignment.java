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
package ngsep.alignments;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.assembly.AlignmentConstantGap;
import ngsep.genome.GenomicRegion;
import ngsep.math.NumberArrays;
import ngsep.sequences.RawRead;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;
/**
 * 
 * @author Jorge Duitama
 *
 */
public class ReadAlignment implements GenomicRegion {

	public static final int FLAG_PAIRED = 0x1;
	public static final int FLAG_PROPER = 0x2;
	public static final int FLAG_READ_UNMAPPED = 0x4;
	public static final int FLAG_MATE_UNMAPPED = 0x8;
	public static final int FLAG_READ_REVERSE_STRAND = 0x10;
	public static final int FLAG_MATE_REVERSE_STRAND = 0x20;
	public static final int FLAG_FIRST_OF_PAIR = 0x40;
	public static final int FLAG_SECOND_OF_PAIR = 0x80;
	public static final int FLAG_SECONDARY = 0x100;
	public static final int FLAG_FAILS_QC = 0x200;
	public static final int FLAG_PCR_DUP = 0x400;
	public static final int FLAG_SUPPLEMENTARY = 0x800;
	//Internal flag to know if the alignment is unique
	public static final int FLAG_MULTIPLE_ALN = 0x1000;
	//Internal flag to know if the pair maps to a different chromosome
	public static final int FLAG_MATE_DIFFERENT_SEQUENCE = 0x2000;
	
	public static final byte ALIGNMENT_DELETION = 1; //D operation in SAM format
	public static final byte ALIGNMENT_MATCH = 3; //M and = operations in SAM format
	public static final byte ALIGNMENT_INSERTION = 2; //I  operation in SAM format
	public static final byte ALIGNMENT_MISMATCH = 7; //X operation in SAM format
	public static final byte ALIGNMENT_SKIPFROMREAD = 6; //S operation in SAM format
	public static final byte ALIGNMENT_SKIPFROMREF = 5; //N operation in SAM format
	public static final byte ALIGNMENT_HARDCLIP = 0; //H operation in SAM format
	public static final byte ALIGNMENT_PADDING = 4; //P operation in SAM format
	
	public static String ALIGNMENT_CHAR_CODES = "HDIMPNSX";
	
	public static final int DEF_MIN_MQ_UNIQUE_ALIGNMENT = 20; 
	
	public static final String DEF_READ_GROUP = "";
	
	public enum Platform {
		ILLUMINA,
		IONTORRENT,
		PACBIO,
		ONT;
		
		public boolean isLongReads () {
			return this == PACBIO || this == ONT;
		}
	}
	
	
	private int sequenceIndex;
	private String sequenceName;
	private int first;
	private int last;
	private int flags = 0;
	private byte alignmentQuality = 0;
    private String mateSequenceName = null;
    private int mateFirst = 0;
    private int inferredInsertSize = 0;
    private short numMismatches = 0;
    
    //Replaces the CIGAR. value / 8 is the length and value % 8 is the operation. See constants above
    private int [] alignment;
    private int readLength = 0; //Read length provided in the constructor or consistent with the alignment
    
    //Quick retrieval data
    private short [] alleleCallLength; //Length of the allele call at each read position. 0 for skipped bases
    private Map<Integer,GenomicVariant> indelCalls; //Indel calls indexed by the last reference position before the event. Null for alignments without indels 
    private boolean alleleCallsUpdated = false;
	
	//Read information loaded on demand 
	private char [] readName=null;
	private CharSequence readCharacters=null;
	private byte [] qualityScores=null;
	private int readNumber;
	
	//Attributes to ignore bases
	private byte basesToIgnoreCloseToIndel = 2;
	private short basesToIgnoreStart = 0;
	private short basesToIgnoreEnd = 0;
	
	//Optional information stored
	private String readGroup = DEF_READ_GROUP;
	
	/**
	 * Creates a read alignment with the given information
	 * @param sequenceName Name of the sequence (e.g. chromosome) to which the read aligned
	 * @param first position of the sequence to which the read aligned
	 * @param last position of the sequence to which the read aligned
	 * @param readLength length of the read
	 * @param flags indicating information of the alignment according to the SAM format
	 */
	public ReadAlignment(String sequenceName, int first, int last, int readLength, int flags) {
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		this.readLength = readLength;
		this.setFlags(flags);
		updateAlleleCallsInfo();
	}
	
	/**
	 * Creates a read alignment with the given information
	 * @param sequenceIndex Index of the sequence to which the read aligned
	 * @param first position of the sequence to which the read aligned
	 * @param last position of the sequence to which the read aligned
	 * @param readLength length of the read
	 * @param flags indicating information of the alignment according to the SAM format
	 */
	public ReadAlignment(int sequenceIndex, int first, int last, int readLength, int flags) {
		this.sequenceIndex = sequenceIndex;
		this.first = first;
		this.last = last;
		this.readLength = readLength;
		this.setFlags(flags);
		updateAlleleCallsInfo();
	}
	
	/**
	 * @return the sequenceIndex
	 */
	public int getSequenceIndex() {
		return sequenceIndex;
	}



	/**
	 * @param sequenceIndex the sequenceIndex to set
	 */
	public void setSequenceIndex(int sequenceIndex) {
		this.sequenceIndex = sequenceIndex;
	}



	/**
	 * @return the readNumber
	 */
	public int getReadNumber() {
		return readNumber;
	}



	/**
	 * @param readNumber the readNumber to set
	 */
	public void setReadNumber(int readNumber) {
		this.readNumber = readNumber;
	}



	@Override
	public String getSequenceName() {
		return sequenceName;
	}
	
	/**
	 * @param sequenceName the sequenceName to set
	 */
	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}

	@Override
	public int getFirst() {
		return first;
	}

	@Override
	public int getLast() {
		if(this.isReadUnmapped()) throw new RuntimeException("Operation not allowed on an unmapped read");
		return last;
	}

	@Override
	public int length() {
		if(this.isReadUnmapped()) throw new RuntimeException("Operation not allowed on an unmapped read");
		return last - first + 1;
	}

	@Override
	public boolean isPositiveStrand() {
		return (flags & FLAG_READ_REVERSE_STRAND)==0;
	}

	@Override
	public boolean isNegativeStrand() {
		return (flags & FLAG_READ_REVERSE_STRAND)!=0;
	}
	
	/**
	 * Changes the negative strand status
	 * @param negativeStrand true if the read aligned to the negative reference strand
	 */
	public void setNegativeStrand (boolean negativeStrand) {
		if(negativeStrand) flags = flags | FLAG_READ_REVERSE_STRAND;
		else flags = flags & ~FLAG_READ_REVERSE_STRAND;
		alleleCallsUpdated = false;
	}
	
	/**
	 * Information of the alignment according to the SAM format
	 * @return int information flags
	 */
	public int getFlags() {
		return flags;
	}

	/**
	 * Changes the flags with information of the alignment
	 * @param flags new flags
	 */
	public void setFlags(int flags) {
		//TODO: Check consistency
		this.flags = flags;
		if(isReadUnmapped()) {
			alignment = null;
			last = 0;
		}
		alleleCallsUpdated = false;
	}
	
	/**
	 * Provides the alignment quality as a phred score
	 * @return short Alignment quality
	 */
	public byte getAlignmentQuality() {
		return alignmentQuality;
	}

	/**
	 * Changes the alignment quality
	 * @param alignmentQuality New Quality as a phred score
	 */
	public void setAlignmentQuality(byte alignmentQuality) {
		this.alignmentQuality = alignmentQuality;
	}
	
	/**
	 * Provides the sequence to which the mate of this read aligned
	 * @return String Reference sequence name
	 */
	public String getMateSequenceName() {
		return mateSequenceName;
	}

	/**
	 * Changes the sequence to which the mate of this read aligned
	 * @param mateSequenceName New reference sequence name for the mate
	 */
	public void setMateSequenceName(String mateSequenceName) {
		this.mateSequenceName = mateSequenceName;
	}

	/**
	 * Provides the position to which the mate of this read aligned
	 * @return int Reference position
	 */
	public int getMateFirst() {
		return mateFirst;
	}

	/**
	 * Changes the position to which the mate of this read aligned
	 * @param mateFirst New reference postion
	 */
	public void setMateFirst(int mateFirst) {
		this.mateFirst = mateFirst;
	}

	/**
	 * Provides the insert size inferred from the alignment of this read and its mate
	 * @return int distance between the start of the left read and the end of the right read.
	 * 0 if one read did not aligned or if the reads aligned to different chromosomes
	 */
	public int getInferredInsertSize() {
		return inferredInsertSize;
	}

	/**
	 * Changes the inferred insert size
	 * @param inferredInsertSize New insert size
	 */
	public void setInferredInsertSize(int inferredInsertSize) {
		this.inferredInsertSize = inferredInsertSize;
	}

	/**
	 * Tells if this read is paired
	 * @return boolean true if the read is paired
	 */
	public boolean isPaired() {
		return (flags & FLAG_PAIRED)!=0;
	}
	/**
	 * Changes the paired status
	 * @param paired true if the read is paired
	 */
	public void setPaired (boolean paired) {
		if(paired) flags = flags | FLAG_PAIRED;
		else flags = flags & ~FLAG_PAIRED;
	}

	/**
	 * Tells if this read and its mate were properly aligned according to an expected distance and orientation
	 * @return boolean true if the read and its mate were properly aligned
	 */
	public boolean isProperPair() {
		return (flags & FLAG_PROPER)!=0;
	}
	/**
	 * Changes the read unmapped status
	 * @param properPair true if the read is not mapped
	 */
	public void setProperPair (boolean properPair) {
		if(properPair) flags = flags | FLAG_PROPER;
		else flags = flags & ~FLAG_PROPER;
	}
	/**
	 * Tells if this alignment is invalid because the read did not actually align to any position of the reference
	 * @return boolean true if the read did not align to the reference
	 */
	public boolean isReadUnmapped() {
		return (flags & FLAG_READ_UNMAPPED)!=0;
	}
	/**
	 * Changes the read unmapped status
	 * @param mateUnmapped true if the read is not mapped
	 */
	public void setReadUnmapped (boolean readUnmapped) {
		if(readUnmapped) flags = flags | FLAG_READ_UNMAPPED;
		else flags = flags & ~FLAG_READ_UNMAPPED;
		alleleCallsUpdated = false;
	}
	/**
	 * Tells if the mate of this read did not align to any position of the reference
	 * @return boolean true if the mate did not align to the reference
	 */
	public boolean isMateUnmapped() {
		return (flags & FLAG_MATE_UNMAPPED)!=0;
	}
	/**
	 * Changes the set mate unmapped status
	 * @param mateUnmapped true if the mate is unmapped
	 */
	public void setMateUnmapped (boolean mateUnmapped) {
		if(mateUnmapped) flags = flags | FLAG_MATE_UNMAPPED;
		else flags = flags & ~FLAG_MATE_UNMAPPED;
	}
	
	/**
	 * Tells if the mate aligned to the positive reference strand
	 * @return true if the mate aligned to the positive reference strand
	 */
	public boolean isMatePositiveStrand() {
		return (flags & FLAG_MATE_REVERSE_STRAND)==0;
	}
	
	/**
	 * Tells if the mate aligned to the negative reference strand
	 * @return true if the mate aligned to the negative reference strand
	 */
	public boolean isMateNegativeStrand() {
		return (flags & FLAG_MATE_REVERSE_STRAND)!=0;
	}
	
	/**
	 * Changes the mate negative strand status
	 * @param mateNegativeStrand true if the mate aligned to the negative reference strand
	 */
	public void setMateNegativeStrand (boolean mateNegativeStrand) {
		if(mateNegativeStrand) flags = flags | FLAG_MATE_REVERSE_STRAND;
		else flags = flags & ~FLAG_MATE_REVERSE_STRAND;
	}
	
	/**
	 * Tells if this read is the first of its pair
	 * @return boolean true if the read is the first of the pair
	 */
	public boolean isFirstOfPair () {
		return (flags & FLAG_FIRST_OF_PAIR)!=0;
	}
	/**
	 * Changes the first of pair status
	 * @param firstOfPair true if the read is the first of the pair
	 */
	public void setFirstOfPair (boolean firstOfPair) {
		if(firstOfPair) flags = flags | FLAG_FIRST_OF_PAIR;
		else flags = flags & ~FLAG_FIRST_OF_PAIR;
	}
	/**
	 * Tells if this read is the second of its pair
	 * @return boolean true if the read is the second of the pair
	 */
	public boolean isSecondOfPair () {
		return (flags & FLAG_SECOND_OF_PAIR)!=0;
	}
	/**
	 * Changes the second of pair status
	 * @param secondOfPair true if the read is the second of the pair
	 */
	public void setSecondOfPair (boolean secondOfPair) {
		if(secondOfPair) flags = flags | FLAG_SECOND_OF_PAIR;
		else flags = flags & ~FLAG_SECOND_OF_PAIR;
	}
	/**
	 * Tells if the alignment is secondary
	 * @return boolean true if the alignment is secondary
	 */
	public boolean isSecondary () {
		return (flags & FLAG_SECONDARY)!=0;
	}
	/**
	 * Changes the secondary status
	 * @param secondary true if the read is the second of the pair
	 */
	public void setSecondary (boolean secondary) {
		if(secondary) flags = flags | FLAG_SECONDARY;
		else flags = flags & ~FLAG_SECONDARY;
	}
	/**
	 * Tells if the alignment fails QC
	 * @return boolean true if the alignment fails QC
	 */
	public boolean failsQC () {
		return (flags & FLAG_FAILS_QC)!=0;
	}
	/**
	 * Tells if the read is a PCR duplicate
	 * @return boolean true if the read is marked as a PCR duplicate
	 */
	public boolean isPCRDuplicate () {
		return (flags & FLAG_PCR_DUP)!=0;
	}
	/**
	 * Tells if the alignment is supplementary
	 * @return boolean true if the alignment is supplementary
	 */
	public boolean isSupplementary () {
		return (flags & FLAG_SUPPLEMENTARY)!=0;
	}
	/**
	 * Tells if the alignment is unique
	 * @return boolean true if this alignment is unique
	 */
	public boolean isUnique () {
		//TODO: Make rules to be unique part of the constructor
		return (flags & FLAG_MULTIPLE_ALN)==0;
	}
	/**
	 * Tells if the mate of this read aligned to a different sequence
	 * @return boolean true if the mate of this read aligned to another sequence
	 */
	public boolean isMateDifferentSequence () {
		return (flags & FLAG_MATE_DIFFERENT_SEQUENCE)!=0;
	}
	/**
	 * Changes the mate in different sequence status
	 * @param mateDifferentSequence true if the mate was aligned in a different sequence
	 */
	public void setMateDifferentSequence (boolean mateDifferentSequence) {
		if(mateDifferentSequence) flags = flags | FLAG_MATE_DIFFERENT_SEQUENCE;
		else flags = flags & ~FLAG_MATE_DIFFERENT_SEQUENCE;
	}
	/**
	 * Tells if the mate of this read aligned to the same sequence
	 * @return boolean true if the mate of this read aligned to the same sequence
	 */
	public boolean isMateSameSequence () {
		return (flags & FLAG_MATE_DIFFERENT_SEQUENCE)==0;
	}
	
	/**
	 * Provides the read name
	 * @return String read name
	 */
	public String getReadName() {
		if(readName == null) return null;
		return new String (readName);
	}

	/**
	 * Changes the read name
	 * @param readName new name
	 */
	public void setReadName(String readName) {
		if(readName == null) this.readName = null;
		this.readName = readName.toCharArray();
	}

	/**
	 * Provides the characters of the read as a CharSequence object
	 * @return CharSequence read characters 
	 */
	public CharSequence getReadCharacters() {
		return readCharacters;
	}

	/**
	 * Changes the characters of this read
	 * @param readCharacters new characters
	 */
	public void setReadCharacters(CharSequence readCharacters) {
		if(readCharacters!=null && readLength>0 &&readCharacters.length()!=readLength) throw new IllegalArgumentException("Input Read length: "+readCharacters.length()+" inconsistent with the expected length "+readLength);		 
		this.readCharacters = readCharacters; 
		if(readCharacters == null) this.qualityScores = null;
		else if (readLength==0) readLength = readCharacters.length();
	}
	
	/**
	 * Provides the length of this read
	 * @return int read length
	 */
	public int getReadLength() {
		return readLength;
	}

	/**
	 * Provides the base quality scores as a String object
	 * @return String quality scores in phred+33 format according to the SAM format specification
	 */
	public String getQualityScores() {
		if(qualityScores == null) return null;
		char [] qs = new char[qualityScores.length];
		for(int i=0;i<qs.length;i++) {
			qs[i] = (char)qualityScores[i];
		}
		return new String(qs);
	}

	/**
	 * Changes the base quality scores
	 * @param qualityScores new quality scores in phred+33 format according to the SAM format specification
	 */
	public void setQualityScores(String qualityScores) {
		if(qualityScores==null) {
			this.qualityScores = null;
			return;
		} 
		int l = qualityScores.length();
		this.qualityScores = new byte [readLength];
		Arrays.fill(this.qualityScores, (byte)38);
		for(int i=0;i<l;i++) {
			int sig = (int) qualityScores.charAt(i);
			if(sig>127) sig = 127;
			this.qualityScores[i] = (byte) sig;  
		}
	}
	
	/**
	 * Changes the basepairs to ignore close to an indel event
	 * @param basesToIgnoreCloseToIndel new number of bases to ignore
	 */
	public void setBasesToIgnoreCloseToIndel(byte basesToIgnoreCloseToIndel) {
		if(basesToIgnoreCloseToIndel<1) throw new IllegalArgumentException("Bases to ignore close to indel must be at least 1");
		this.basesToIgnoreCloseToIndel = basesToIgnoreCloseToIndel;
		this.alleleCallsUpdated = false;
	}
	
	/**
	 * Changes the base pairs to ignore at the 5 prime end of this read
	 * @param basesToIgnore5P new number of base pairs to ignore
	 */
	public void setBasesToIgnore5P(short basesToIgnore5P) {
		if(isNegativeStrand()) this.basesToIgnoreEnd = basesToIgnore5P;
		else this.basesToIgnoreStart = basesToIgnore5P;
		this.alleleCallsUpdated = false;
	}

	/**
	 * Changes the base pairs to ignore at the 3 prime end of this read
	 * @param basesToIgnore3P new number of base pairs to ignore
	 */
	public void setBasesToIgnore3P(short basesToIgnore3P) {
		if(isNegativeStrand()) this.basesToIgnoreStart = basesToIgnore3P;
		else this.basesToIgnoreEnd = basesToIgnore3P;
		this.alleleCallsUpdated = false;
	}
	
	/**
	 * Provides the base pairs to ignore at the start of this alignment
	 * @return short basePairs to ignore
	 */
	public short getBasesToIgnoreStart() {
		return basesToIgnoreStart;
	}
	/**
	 * Changes the base pairs to ignore at the start of this alignment
	 * @param basesToIgnoreStart new number of base pairs to ignore
	 */
	public void setBasesToIgnoreStart(short basesToIgnoreStart) {
		this.basesToIgnoreStart = basesToIgnoreStart;
		this.alleleCallsUpdated = false;
	}

	/**
	 * Provides the base pairs to ignore at the end of the alignment
	 * @return short base pairs to ignore
	 */
	public short getBasesToIgnoreEnd() {
		return basesToIgnoreEnd;
	}
	/**
	 * Changes the base pairs to ignore at the end of this alignment
	 * @param basesToIgnoreEnd new number of base pairs to ignore
	 */
	public void setBasesToIgnoreEnd(short basesToIgnoreEnd) {
		this.basesToIgnoreEnd = basesToIgnoreEnd;
		this.alleleCallsUpdated = false;
	}
	/**
	 * Provides the read group of this read
	 * @return String read group
	 */
	public String getReadGroup() {
		return readGroup;
	}

	/**
	 * Changes the read group of this read
	 * @param readGroup New read group
	 */
	public void setReadGroup(String readGroup) {
		if(readGroup!=null) this.readGroup = readGroup;
		else this.readGroup = DEF_READ_GROUP;
	}
	
	/**
	 * @return Approximate number of mismatches
	 */
	public short getNumMismatches() {
		return numMismatches;
	}

	/**
	 * Changes the approximate number of mismatches
	 * @param numMismatches New number of mismatches
	 */
	public void setNumMismatches(short numMismatches) {
		this.numMismatches = numMismatches;
	}

	/**
	 * Tells if the alignment is partial (if it has soft clipped basepairs)
	 * @param minClipLength Minimum length of the clip to be considered 
	 * @return boolean true if the alignment starts or ends with softclipped basepairs 
	 */
	public boolean isPartialAlignment(int minClipLength) {
	
		if(alignment == null || alignment.length<2) return false;
		int s1 = 0;
		byte op1 = getOperator(alignment[0]); 
		if(op1==ALIGNMENT_SKIPFROMREAD || op1==ALIGNMENT_HARDCLIP) s1 = getOperationLength(alignment[0]);
		int s2 = 0;
		int l = alignment.length-1;
		byte op2 = getOperator(alignment[l]);
		if(op2==ALIGNMENT_SKIPFROMREAD || op2==ALIGNMENT_HARDCLIP) s2 = getOperationLength(alignment[l]);
		return s1>=minClipLength || s2>=minClipLength;
	}

	private void updateAlleleCallsInfo() {
		if(alleleCallsUpdated) return;
		if(alignment == null) {
			alleleCallLength = null;
			indelCalls = null;
			alleleCallsUpdated = true;
			return;
		}
		int currentRefPos = first;
		int currentReadPos = 0;
		//Reset quick retrieval info
		alleleCallLength = new short [readLength];
		Arrays.fill(alleleCallLength, (short)0);
		indelCalls = null;
		boolean previousIsIndel=false;
		for(int i=0;i<alignment.length;i++) {
			int length = getOperationLength(alignment[i]);
			byte operator = getOperator(alignment[i]);
			boolean cRef = consumesReferenceBases(alignment[i]);
			boolean cRead = consumesReadBases(alignment[i]);
			byte nextOperator = -1;
			int nextOpLen = 0;
			boolean nextIsIndel = false;
			int nextReadConsumption = 0;
			if(i<alignment.length-1) {
				nextOperator = getOperator(alignment[i+1]);
				nextOpLen = getOperationLength(alignment[i+1]);
				nextIsIndel = isIndel(nextOperator);
				nextReadConsumption = consumesReadBases(nextOperator)?nextOpLen:0;
			}
			if(cRef) {
				if(cRead) {
					//Fill pileups with matches or mismatches. Check in the last position if an indel is coming to call it 
					for(int j=0;j<length;j++) {
						//Skip bases too close to the start
						boolean skip = currentReadPos < basesToIgnoreStart;
						//Skip bases close to the end
						skip = skip || (readLength - currentReadPos) <= basesToIgnoreEnd;
						//Skip bases right after the previous indel event
						skip = skip || previousIsIndel && j<basesToIgnoreCloseToIndel;
						//Skip bases before the next indel event
						skip = skip || (nextIsIndel && j<length-1 && j>=length-basesToIgnoreCloseToIndel);
						//Skip the the indel itself if too close to the end of the read
						skip = skip || (nextIsIndel && j==length-1 && (currentReadPos < basesToIgnoreCloseToIndel || readLength - currentReadPos - nextReadConsumption < basesToIgnoreCloseToIndel));
						//Skip the indel if the end falls into the bases to ignore
						int readPosAfterIndel = currentReadPos + nextReadConsumption + 1;
						skip = skip || (nextIsIndel && j==length-1 && (readLength - readPosAfterIndel < basesToIgnoreEnd));
						//System.out.println("Read id: "+aln.getReadName()+". Read start: "+aln.getAlignmentStart()+". Read pos: "+readPos+ " ReferencePos: "+referencePos+" skip: "+skip+" IgnoreStart: "+basesToIgnoreStart+" IgnoreEnd: "+basesToIgnoreEnd+" Reversed: "+aln.getReadNegativeStrandFlag()+" call "+read.charAt(readPos)+" CIGAR: "+aln.getCigarString());
						if(!skip) {
							if (j==length-1 && nextIsIndel) {
								int refLast = currentRefPos+1;
								if(nextOperator==ALIGNMENT_INSERTION) {
									//Process insertion call. TODO: Verify that nextOpLen < 2E15
									alleleCallLength[currentReadPos] = (short) (nextOpLen+2);
								} else {
									//Process deletion call
									alleleCallLength[currentReadPos] = (short)2;
									refLast += nextOpLen;
								}
								if(indelCalls==null) {
									indelCalls = new TreeMap<Integer,GenomicVariant>();
								}
								GenomicVariantImpl indel = new GenomicVariantImpl(sequenceName, currentRefPos, refLast, GenomicVariant.TYPE_INDEL);
								indel.setLength(nextOpLen);
								indelCalls.put(currentRefPos, indel);
							} else {
								//Process match or mismatch call
								alleleCallLength[currentReadPos] = (short)1;
							} 
						}
						currentRefPos++;
						currentReadPos++;
					}
				} else {
					currentRefPos+=length;
				}
			} else if(cRead) {
				currentReadPos+=length;
			}
			previousIsIndel = isIndel(operator);
		}
		alleleCallsUpdated = true;
		//throw new RuntimeException("Called method to update allele calls");
	}

	/**
	 * Returns the read position that corresponds to the given reference coordinate
	 * @param referencePos 1-based reference coordinate
	 * @return 0-based position of the read aligning with the given reference, starting from the start of the alignment.
	 * -1 if the reference position does not align with the read
	 */
	public int getAlignedReadPosition (int referencePos) {
		failIfReadUnmappedOrInconsistentAlignment();
		
		int currentRefPos = first;
		int currentReadPos = 0;
		if(referencePos<first || referencePos>last) return -1;
		for(int i=0;i<alignment.length;i++) {
			int length = getOperationLength(alignment[i]);
			boolean cRef = consumesReferenceBases(alignment[i]);
			boolean cRead = consumesReadBases(alignment[i]);
			if(cRef && cRead) {
				if(referencePos<currentRefPos) return -1;
				else if(currentRefPos+length>referencePos) {
					int answer = currentReadPos + referencePos-currentRefPos;
					if(answer <0 || answer >= readLength) {
						System.err.println("WARN: Inconsistent CIGAR for read "+getReadName()+" at "+first+" "+last+" length: "+readLength+" CIGAR: "+getCigarString());
						return -1;
					}
					return answer; 
				}
			}
			if(cRef) {
				currentRefPos += length;
			}
			if(cRead) {
				currentReadPos += length;
			}
		}
		return -1;
	}

	private void failIfReadUnmappedOrInconsistentAlignment() {
		if(this.isReadUnmapped()) throw new RuntimeException("Operation not allowed on an unmapped read");
		if(alignment==null) throw new RuntimeException("Inconsistent null alignment for an unmapped read. Read name: "+new String(readName));
	}
	/**
	 * Returns the reference position that corresponds to the given read position
	 * @param readPos 0-based position of the original read before alignment
	 * @return int 1-based reference coordinate aligning with the given position. -1 if the read position does not align with the reference
	 */
	public int getReferencePosition (int readPos) {
		failIfReadUnmappedOrInconsistentAlignment();
		if(!isNegativeStrand()) {
			return getReferencePositionAlignedRead(readPos);
		} else {
			return getReferencePositionReverse(readPos);
		}
	}

	/**
	 * Returns the reference position that corresponds to the given position of the aligned read
	 * @param alignedReadPos 0-based position of the read, as aligned to the reference (left to right)
	 * @return int 1-based reference coordinate aligning with the given position. -1 if the read position does not align with the reference
	 */
	public int getReferencePositionAlignedRead(int alignedReadPos) {
		int currentRefPos = first;
		int currentAlnPos = 0;
		for(int i=0;i<alignment.length;i++) {
			int length = getOperationLength(alignment[i]);
			boolean cRef = consumesReferenceBases(alignment[i]);
			boolean cRead = consumesReadBases(alignment[i]);
			if(cRef && cRead) {
				if(alignedReadPos<currentAlnPos) return -1;
				else if(currentAlnPos+length>alignedReadPos) {
					return currentRefPos + alignedReadPos-currentAlnPos; 
				}
			}
			if(cRef) {
				currentRefPos += length;
			}
			if(cRead) {
				currentAlnPos += length;
			}
		}
		return -1;
	}
	
	private int getReferencePositionReverse (int readPos) {
		failIfReadUnmappedOrInconsistentAlignment();
		int currentRefPos = last;
		int currentReadPos = 0;
		for(int i=alignment.length-1;i>=0;i--) {
			int length = getOperationLength(alignment[i]);
			boolean cRef = consumesReferenceBases(alignment[i]);
			boolean cRead = consumesReadBases(alignment[i]);
			if(cRef && cRead) {
				if(readPos<currentReadPos) return -1;
				else if(currentReadPos-length<readPos) {
					return currentRefPos - readPos-currentReadPos; 
				}
			}
			if(cRef) {
				currentRefPos -= length;
			}
			if(cRead) {
				currentReadPos += length;
			}
		}
		return -1;
	}
	
	
	/**
	 * Provides the allele call (if any) at the given reference position taking into account possible
	 * base pairs to ignore
	 * @param referencePos Position to consider in the reference sequence to which this read aligned
	 * @return CharSequence Object with the base pair(s) starting from the reference position.
	 * It can be more than one character, especially in the case of insertions.  
	 */
	public CharSequence getAlleleCall (int referencePos) {
		if(readCharacters ==null) return null;
		int readPos = getAlignedReadPosition(referencePos);
		if(readPos<0) return null;
		updateAlleleCallsInfo();
		//if(readPos>=alleleCallLength.length) System.err.println("Invalid read position: "+readPos+" Read: "+getReadName()+" Located at "+sequenceName+":"+first+" last:"+last+". Length: "+readCharacters.length()+" CIGAR: "+getCigarString()+" expected length: "+expectedReadLength+" reference pos: "+referencePos);
		//if(readPos>=alleleCallLength.length) System.err.println("Alignment: "+alignment[0]+" "+alignment[1]+" "+alignment[2]);
		int length = alleleCallLength[readPos];
		if(length == 0) return null;
		//if(referencePos==-1) System.out.println("ReadAlignment. Read id: "+getReadName()+". Sequence: "+readCharacters.toString()+". readpos: "+readPos+". end: "+(readPos+length)+" subseq: "+readCharacters.subSequence(readPos, readPos+length)+". length subseq: "+readCharacters.subSequence(readPos, readPos+length).length());
		return readCharacters.subSequence(readPos, readPos+length);
	}
	/**
	 * Provides the allele call (if any) at the given reference coordinates taking into account possible
	 * base pairs to ignore
	 * @param referenceFirst First position to consider in the reference sequence to which this read aligned
	 * @param referenceFirst Last position to consider in the reference sequence to which this read aligned
	 * @return CharSequence Object with the base pair(s) aligning between the given reference positions.
	 */
	public CharSequence getAlleleCall (int referenceFirst, int referenceLast) {
		if(readCharacters == null) return null;
		updateAlleleCallsInfo();
		int readFirst = getAlignedReadPosition(referenceFirst);
		int readLast = getAlignedReadPosition(referenceLast);
		if(readFirst<0 || readLast<0 || readLast < readFirst) return null;
		if(withinIgnoreRegions(readFirst, readLast)) return null;
		return readCharacters.subSequence(readFirst, readLast+1);
	}
	/**
	 * Returns the base quality score in phred+33 format of the base pair aligning to the given position
	 * @param referencePos Position to consider in the reference sequence to which this read aligned
	 * @return char base quality score
	 */
	public char getBaseQualityScore (int referencePos) {
		int readPos = getAlignedReadPosition(referencePos);
		if(readPos<0) return 33;
		if(qualityScores == null) return '+';
		return (char) qualityScores[readPos];
	}
	/**
	 * Returns the base quality scores in phred+33 format of the base pairs aligning to the given coordinates
	 * @param referenceFirst First position to consider in the reference sequence to which this read aligned
	 * @param referenceFirst Last position to consider in the reference sequence to which this read aligned
	 * @return String base quality scores
	 */
	public String getBaseQualityScores (int referenceFirst, int referenceLast) {
		int readFirst = getAlignedReadPosition(referenceFirst);
		int readLast = getAlignedReadPosition(referenceLast);
		if(readFirst<0 || readLast<0  || readLast < readFirst) return null;
		if(withinIgnoreRegions(readFirst, readLast)) return null;
		if(qualityScores == null) return RawRead.generateFixedQSString('+', readLast-readFirst+1);
		return getQualityScores().substring(readFirst,readLast+1);
	}
	private boolean withinIgnoreRegions (int readFirst, int readLast) {
		return readFirst<basesToIgnoreStart || readLength - readLast <= basesToIgnoreEnd;
	}
	/**
	 * Returns start sites in this alignment for indel events
	 * @return Map<Integer,GenomicVariant> Map with reference positions as keys and indel events as values
	 */
	public Map<Integer,GenomicVariant> getIndelCalls () {
		failIfReadUnmappedOrInconsistentAlignment();
		updateAlleleCallsInfo();
		return indelCalls;
	}
	/**
	 * Provides the indel call (if any) at the given reference position
	 * @param referencePos Position to consider in the reference sequence to which this read aligned
	 * @return GenomicVariant object with the indel call
	 */
	public GenomicVariant getIndelCall(int referencePos) {
		failIfReadUnmappedOrInconsistentAlignment();
		updateAlleleCallsInfo();
		if(indelCalls == null) return null;
		return indelCalls.get(referencePos);
	}
	
	/**
	 * Modifies this alignment moving the start of a given indel event
	 * @param indelRefPos Start position of the indel event to move
	 * @param newIndelRefPos New start position of the indel event
	 * @return boolean true if the indel could be safely moved
	 */
	public boolean moveIndelStart (int indelRefPos, int newIndelRefPos) {
		int posPrint = -1;
		failIfReadUnmappedOrInconsistentAlignment();
		int displacement = newIndelRefPos - indelRefPos;
		if(first == posPrint) System.out.println("Read "+getReadName()+" CIGAR: "+getCigarString()+" Displacement: "+displacement+" new pos: "+newIndelRefPos);
		if(displacement == 0) return true;
		
		int [] newAlignment = new int [alignment.length];
		int indelNextIdx = -1;
		int currentRefPos = first;
		for(int i=0;i<alignment.length;i++) {
			int length = getOperationLength(alignment[i]);
			byte operator = getOperator(alignment[i]);
			if(isIndel(operator) && currentRefPos == indelRefPos+1) {
				if(i==0 || i==alignment.length -1 ) return false;
				byte opBefore = getOperator(alignment[i-1]);
				int lenBefore = getOperationLength(alignment[i-1]);
				if(!consumesReadBases(opBefore) || !consumesReferenceBases(opBefore) || lenBefore<=-displacement) return false;
				byte opAfter = getOperator(alignment[i+1]);
				int lenAfter = getOperationLength(alignment[i+1]);
				if(!consumesReadBases(opAfter) || !consumesReferenceBases(opAfter) || lenAfter<=displacement) return false;
				if(first == posPrint) System.out.println("Op before: "+lenBefore+" "+opBefore+" Op after: "+lenAfter+" "+opAfter);
				newAlignment[i-1] = getAlnValue(lenBefore+displacement, opBefore);
				newAlignment[i] = alignment[i];
				newAlignment[i+1] = getAlnValue(lenAfter-displacement, opAfter);
				indelNextIdx = i+1;
			} else if (i!=indelNextIdx) {
				//Avoids resetting to the original value the next element after the indel
				newAlignment[i] = alignment[i];
			}
			if(consumesReferenceBases(alignment[i])) {
				currentRefPos += length;
			}
		}
		if(indelNextIdx<0) return false;
		if(first == posPrint) System.out.println("New alignment for read: "+getReadName()+": "+newAlignment[0]+" "+newAlignment[1]+" "+newAlignment[2]);
		alignment = newAlignment;
		alleleCallsUpdated = false;
		return true;
	}
	/**
	 * Gives the number of items making the alignment encoded by the CIGAR
	 * For example, the number of items of the CIGAR 50M10I40M is 3
	 * @return int number of CIGAR items 
	 */
	public int getNumCigarItems () {
		return alignment.length;
	}
	/**
	 * Returns the length of the CIGAR item with the given index
	 * @param itemIndex (0-based) index of the item to query
	 * @return int length of the element at position index of the CIGAR
	 */
	public int getCigarItemLength (int itemIndex) {
		return getOperationLength(alignment[itemIndex]);
	}
	
	/**
	 * Returns the operator of the CIGAR item with the given index
	 * @param itemIndex (0-based) index of the item to query
	 * @return int operator of item at position index of the CIGAR. See constants ALIGNMENT*
	 */
	public int getCigarItemOperator (int itemIndex) {
		return getOperator(alignment[itemIndex]);
	}
	//Alignment encoding / decoding methods
	public static int getOperationLength(int alnValue) {
		return alnValue/8;
	}

	public static byte getOperator(int alnValue) {
		return (byte)(alnValue & 0x7);
	}
	
	public static boolean consumesReadBases (int alnValue) {
		return (alnValue & 0x2) != 0;
	}
	
	public static boolean consumesReferenceBases (int alnValue) {
		return (alnValue & 0x1) != 0;
	}
	
	public static int getAlnValue (int operationLength , byte operator) {
		return 8*operationLength + operator;
	}
	
	public static boolean isIndel(byte operator) {
		return operator == ALIGNMENT_DELETION || operator == ALIGNMENT_INSERTION;
	}
	
	/**
	 * Provides the CIGAR describing the alignment as a String object
	 * @return String CIGAR
	 */
	public String getCigarString() {
		StringBuilder cigar = new StringBuilder();
		for(int i=0;i<alignment.length;i++) {
			int length = getOperationLength(alignment[i]);
			byte operator = getOperator(alignment[i]);
			cigar.append(length);
			cigar.append(ALIGNMENT_CHAR_CODES.charAt(operator));
		}
		return cigar.toString();
	}
	/**
	 * Changes this alignment according to the given CIGAR String
	 * @param cigarString object describing the new alignment
	 */
	public void setCigarString(String cigarString) {
		int posPrint = -1;
		if(cigarString==null) {
			if(!isReadUnmapped()) throw new IllegalArgumentException("Invalid null CIGAR for an aligned read at: "+sequenceName+":"+first);
			return;
		}
		if(first == posPrint) System.out.println(" CIGAR: "+cigarString);
		List<Integer> alignmentList = new ArrayList<Integer>();
		int nextLengthStart=0;
		int expectedReadLength = 0;
		int expectedEnd = first;
		for(int i=0;i<cigarString.length();i++) {
			char c = cigarString.charAt(i);
			if(!Character.isDigit(c)) {
				if(i==nextLengthStart) new IllegalArgumentException("Malformed CIGAR. Operator: "+c+" is not preceded by a length");
				String lengthStr = cigarString.substring(nextLengthStart,i);
				int length;
				try {
					length = Integer.parseInt(lengthStr);
				} catch (NumberFormatException e) {
					throw new IllegalArgumentException("Malformed CIGAR:"+cigarString+". Invalid length: "+lengthStr+". CIGAR substring: "+nextLengthStart+"-"+i,e);
				}
				if(first == posPrint) System.out.println("Next CIGAR length: "+length+" operation: "+c);
				byte operator = (byte) ALIGNMENT_CHAR_CODES.indexOf(c);
				if(operator<0 && c == '=') operator = ALIGNMENT_MATCH;
				if(operator==-1) new IllegalArgumentException("Malformed CIGAR. Unrecognized operator: "+c);
				alignmentList.add(getAlnValue(length, operator));
				if(consumesReadBases(operator)) expectedReadLength+=length;
				if(consumesReferenceBases(operator)) expectedEnd+=length;
				nextLengthStart = i+1;
			}
		}
		if(readCharacters!=null && expectedReadLength!=readCharacters.length())  new IllegalArgumentException("Malformed CIGAR. Expected read length: "+expectedReadLength+" different than read length: "+readCharacters.length());
		alignment = NumberArrays.toIntArray(collapseEqualEvents (alignmentList));
		if(expectedReadLength != this.readLength) {
			System.out.println("WARN. New CIGAR changes read length for read at "+sequenceName+":"+first+" current length: "+this.readLength+" expected length: "+expectedReadLength);
			this.readLength = expectedReadLength;
		}
		if(expectedEnd -1 != this.last) {
			System.out.println("WARN. New CIGAR changes last alignment position for read at "+sequenceName+":"+first+" currentLast: "+this.last+" expectedLast: "+(expectedEnd-1));
			this.last = expectedEnd - 1;
		}
		if(first == posPrint && alignment.length >=3) System.out.println("Alignment codes: "+alignment[0]+" "+alignment[1]+" "+alignment[2]+" Expected read length: "+expectedReadLength);
		alleleCallsUpdated = false;
	}
	
	/**
	 * Changes the current alignment of this read
	 * @param alignmentCodesList List with encoded alignment operations. The last three bits of each number encode the alignment operation.
	 * The remaining bits encode the length
	 */
	public void setAlignment (List<Integer> alignmentCodesList) {
		alignment = NumberArrays.toIntArray(collapseEqualEvents (alignmentCodesList));
		//Check consistency of alignment last position and expected read length
		int expectedReadLength = 0;
		int expectedEnd = first;
		for (int i=0;i<alignment.length;i++) {
			int nextAlnCode = alignment[i];
			int nextLength = getOperationLength(nextAlnCode);
			if(consumesReferenceBases(nextAlnCode))expectedEnd+=nextLength;
			if(consumesReadBases(nextAlnCode))expectedReadLength+=nextLength;
			
		}
		if(expectedReadLength != this.readLength) {
			System.out.println("WARN. New alignment changes read length for read at "+sequenceName+":"+first+" current length: "+this.readLength+" expected length: "+expectedReadLength);
			this.readLength = expectedReadLength;
		}
		if(expectedEnd -1 != this.last) {
			System.out.println("WARN. New alignment changes last alignment position for read at "+sequenceName+":"+first+" currentLast: "+this.last+" expectedLast: "+(expectedEnd-1));
			this.last = expectedEnd - 1;
		}
		alleleCallsUpdated = false;
	}
	
	private List<Integer> collapseEqualEvents(List<Integer> alignmentList) {
		List<Integer> answer = new ArrayList<Integer>();
		byte lastOperator=-1;
		int totalLength = 0;
		for(int alnValue:alignmentList) {
			int length = getOperationLength(alnValue);
			byte op = getOperator(alnValue);
			if(op!=lastOperator) {
				if(totalLength>0) answer.add(getAlnValue(totalLength, lastOperator));
				totalLength = 0;
				lastOperator = op;
			}
			totalLength+=length;
		}
		if(totalLength>0) answer.add(getAlnValue(totalLength, lastOperator));
		return answer;
	}
	
	public int getSoftClipStart() {
		if(getCigarItemOperator(0) == ALIGNMENT_SKIPFROMREAD) {
			return getCigarItemLength(0);
		}
		return 0;
	}

	public int getSoftClipEnd() {
		if(getCigarItemOperator(alignment.length-1) == ALIGNMENT_SKIPFROMREAD ) {
			return getCigarItemLength(alignment.length-1);
		}
		return 0;
	}

	/**
	 * Realigns the start of this alignment
	 * @param newAlnFirst 1-based New reference coordinate matching the start of the read
	 * @param firstMatchLength Length of the match that will become the new first element of the alignment
	 * @param refPosAfter 1-based reference coordinate of the first position that should remain with the same alignment.
	 * @param alnPosAfter 0-based first alignment position that should remain with the same alignment. It must be larger than firstMatchLength 
	 */
	public void realignStart(int newAlnFirst, int firstMatchLength, int refPosAfter, int alnReadPosAfter) {
		int posPrint = -1;
		assert alnReadPosAfter>=firstMatchLength;
		if(first==posPrint) System.out.println("Read "+getReadName()+" at "+first+"-"+last+" CIGAR: "+getCigarString()+" New aln first: "+newAlnFirst+" first match: "+firstMatchLength+" aln pos after: "+alnReadPosAfter);
		List<Integer> alignmentList = new ArrayList<Integer>();
		alignmentList.add(getAlnValue(firstMatchLength, ALIGNMENT_MATCH));
		int nextRefPos = newAlnFirst+firstMatchLength;
		int unknownBpRead = alnReadPosAfter-firstMatchLength;
		int unknownBpRef = refPosAfter-nextRefPos;
		if(first==posPrint) System.out.println("Read "+getReadName()+" Next ref pos: "+nextRefPos+" unknown bp read: "+unknownBpRead+" ref pos after: "+refPosAfter+" unknown bp ref: "+unknownBpRef);
		if(unknownBpRef<0 || unknownBpRead<0) {
			System.err.println("WARN: Can not realing start of alignment for read: "+getReadName()+" at "+first+"-"+last+". Reference first: "+nextRefPos+" Reference after: "+refPosAfter);
			return;
		}
		int difference = unknownBpRead - unknownBpRef;
		if(difference==0) {
			if(unknownBpRead>0 ) alignmentList.add(getAlnValue(unknownBpRead, ALIGNMENT_MATCH));
		} else if (difference>0) {
			alignmentList.add(getAlnValue(difference, ALIGNMENT_INSERTION));
			if(unknownBpRef>0) alignmentList.add(getAlnValue(unknownBpRef, ALIGNMENT_MATCH));
		} else {
			alignmentList.add(getAlnValue(-difference, ALIGNMENT_DELETION));
			if(unknownBpRead>0) alignmentList.add(getAlnValue(unknownBpRead, ALIGNMENT_MATCH));
		}
		
		int currentReadPos = 0;
		boolean copy = false;
		for(int i=0;i<alignment.length;i++) {
			int length = getOperationLength(alignment[i]);
			byte operator = getOperator(alignment[i]);
			boolean cRead = consumesReadBases(alignment[i]);
			if(copy) {
				alignmentList.add(alignment[i]);
			}
			if(cRead) {
				if(!copy && alnReadPosAfter<currentReadPos+length) {
					copy = true;
					int diff = currentReadPos+length-alnReadPosAfter;
					alignmentList.add(getAlnValue(diff, operator));
				}
				currentReadPos += length;
			}
		}
		first = newAlnFirst;
		alignment = NumberArrays.toIntArray(alignmentList);
		alleleCallsUpdated = false;
	}

	/**
	 * Realigns the end of this alignment
	 * @param refPosBefore 1-based genomic coordinate of the last position of the reference that should not be realigned
	 * @param alnPosBefore 0-based last aligned read position that should not be realigned
	 * @param finalMatchRefStart 1-based position of the reference for the last match
	 * @param finalMatchLength Number of match bp that must appear at the end. it must be smaller than length - readPosBefore 
	 */
	public void realignEnd(int refPosBefore, int alnPosBefore, int finalMatchRefStart, int finalMatchLength) {
		int bpEndRead = readLength - alnPosBefore - 1;
		assert  bpEndRead>=finalMatchLength;
		List<Integer> alignmentList = new ArrayList<Integer>();
		int currentReadPos = 0;
		boolean copy = true;
		for(int i=0;i<alignment.length;i++) {
			int length = getOperationLength(alignment[i]);
			byte operator = getOperator(alignment[i]);
			boolean cRead = consumesReadBases(alignment[i]);
			if(cRead) {
				if(copy && alnPosBefore<currentReadPos+length) {
					int diff = currentReadPos+length-alnPosBefore-1;
					alignmentList.add(getAlnValue(length-diff, operator));
					copy = false;
				}
				currentReadPos += length;
			}
			if(copy) {
				alignmentList.add(alignment[i]);
			}
		}
		int unknownBpRef = finalMatchRefStart-refPosBefore-1;
		int unknownBpRead = bpEndRead - finalMatchLength;
		if(unknownBpRef<0 || unknownBpRead<0) {
			System.err.println("WARN: Can not realing end of alignment for read: "+getReadName()+" at "+first+"-"+last+". Reference before: "+refPosBefore+" final match start: "+finalMatchRefStart);
			return;
		}
		int difference = unknownBpRead - unknownBpRef;
		if(difference==0) {
			if(unknownBpRead>0 ) alignmentList.add(getAlnValue(unknownBpRead, ALIGNMENT_MATCH));
		} else if (difference>0) {
			alignmentList.add(getAlnValue(difference, ALIGNMENT_INSERTION));
			if(unknownBpRef>0) alignmentList.add(getAlnValue(unknownBpRef, ALIGNMENT_MATCH));
		} else {
			alignmentList.add(getAlnValue(-difference, ALIGNMENT_DELETION));
			if(unknownBpRead>0) alignmentList.add(getAlnValue(unknownBpRead, ALIGNMENT_MATCH));
		}
		alignmentList.add(getAlnValue(finalMatchLength, ALIGNMENT_MATCH));
		alignment = NumberArrays.toIntArray(alignmentList);
		last = finalMatchRefStart + finalMatchLength -1;
		alleleCallsUpdated = false;
	}

	public boolean hasIndelCalls(int referenceFirst, int referenceLast) {
		Map<Integer, GenomicVariant> calls = getIndelCalls();
		if (calls==null) return false;
		for(int i: calls.keySet()) {
			if(i>=referenceFirst && i<=referenceLast) return true;
		}
		return false;
	}
	public static LinkedList<Integer> encodePairwiseAlignment(String [] queryToSubjectAln ) {
		LinkedList<Integer> answer = new LinkedList<Integer>();
		byte nextOperator = 0;
		int nextLength = 0;
		String queryAln = queryToSubjectAln[0];
		String subjectAln = queryToSubjectAln[1];
		if(queryAln.length()!=subjectAln.length()) throw new IllegalArgumentException("Inconsistent length of alignment. Query length: "+queryAln.length()+" subject length: "+subjectAln.length());
		for(int i=0;i<subjectAln.length();i++) {
			char subjectChar = subjectAln.charAt(i);
			char queryChar = queryAln.charAt(i);
			byte op = ReadAlignment.ALIGNMENT_MATCH;
			if(subjectChar == AlignmentConstantGap.GAP_CHARACTER) {
				op = ReadAlignment.ALIGNMENT_INSERTION;
			} else if(queryChar == AlignmentConstantGap.GAP_CHARACTER) {
				op = ReadAlignment.ALIGNMENT_DELETION;
			}
			if(op != nextOperator) {
				if(nextLength>0) {
					answer.add(getAlnValue(nextLength, nextOperator));
				}
				nextOperator = op;
				nextLength = 0;
			}
			nextLength++;
		}
		if(nextLength>0) {
			answer.add(getAlnValue(nextLength, nextOperator));
		}
		return answer;
	}
	public boolean clipBorders(int minMatchLength) {
		LinkedList<Integer> alnList = new LinkedList<Integer>();
		int removedBpSequenceStart=0;
		int removedBpQueryStart=0;
		int i=0;
		for (;i<alignment.length;i++) {
			int code = alignment[i];
			int length = getOperationLength(code);
			byte op = getOperator(code);
			if(length>=minMatchLength && op== ALIGNMENT_MATCH) break;
			
			if(consumesReadBases(code)) removedBpQueryStart+=length;
			if(consumesReferenceBases(code)) removedBpSequenceStart+=length;
		}
		
		
		
		int removedBpSequenceEnd=0;
		int removedBpQueryEnd=0;
		int j=alignment.length-1;
		for(;j>=0;j--) {
			int code = alignment[j];
			int length = ReadAlignment.getOperationLength(code);
			byte op = ReadAlignment.getOperator(code);
			if(length>=minMatchLength && op== ReadAlignment.ALIGNMENT_MATCH) break;
			if(consumesReadBases(code)) removedBpQueryEnd+=length;
			if(consumesReferenceBases(code)) removedBpSequenceEnd+=length;
			
		}
		if(i>j) {
			//System.err.println("WARN: Failed to clip ends of read alignment at "+sequenceName+":"+first+" CIGAR: "+getCigarString());
			return false;
		}
		
		if(removedBpSequenceStart>0) first+=removedBpSequenceStart;
		if(removedBpSequenceEnd>0) last-=removedBpSequenceEnd;
		if(removedBpQueryStart>0) alnList.add(ReadAlignment.getAlnValue(removedBpQueryStart, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
		for(;i<=j;i++) {
			alnList.add(alignment[i]);
		}
		if(removedBpQueryEnd>0) alnList.add(ReadAlignment.getAlnValue(removedBpQueryEnd, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
		setAlignment(alnList);
		return true;
	}
}
