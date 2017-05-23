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
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ngsep.genome.GenomicRegion;
import ngsep.math.NumberArrays;
import ngsep.variants.GenomicVariant;
import ngsep.variants.GenomicVariantImpl;

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
	
	public static byte ALIGNMENT_DELETION = 1; //D operation in SAM format
	public static byte ALIGNMENT_MATCH = 3; //M and = operations in SAM format
	public static byte ALIGNMENT_INSERTION = 2; //I  operation in SAM format
	public static byte ALIGNMENT_MISMATCH = 7; //X operation in SAM format
	public static byte ALIGNMENT_SKIPFROMREAD = 6; //S operation in SAM format
	public static byte ALIGNMENT_SKIPFROMREF = 5; //N operation in SAM format
	public static byte ALIGNMENT_HARDCLIP = 0; //H operation in SAM format
	public static byte ALIGNMENT_PADDING = 4; //P operation in SAM format
	
	public static String ALIGNMENT_CHAR_CODES = "HDIMPNSX";
	
	
	
	private String sequenceName;
	private int first;
	private int last;
	private int flags = 0;
	private short alignmentQuality = 0;
    private String mateSequenceName = null;
    private int mateFirst = 0;
    private int inferredInsertSize = 0;
    
    //Replaces the CIGAR. value / 8 is the length and value % 8 is the operation. See constants above
    private int [] alignment;
    private int readLength = 0; //Read length provided in the constructor or consistent with the alignment
    
    //Quick retrieval data
    private short [] alleleCallLength; //Length of the allele call at each read position. 0 for skipped bases
    private Map<Integer,GenomicVariant> indelCalls; //Indel calls indexed by the last reference position before the event. Null for alignments without indels 
    private boolean alleleCallsUpdated = false;
    
	
	//Read information. Object not used to save memory 
	private char [] readName=null;
	private CharSequence readCharacters=null;
	private char [] qualityScores=null;
	
	//Attributes to ignore bases
	private byte basesToIgnoreCloseToIndel = 2;
	private short basesToIgnoreStart = 0;
	private short basesToIgnoreEnd = 0;
	
	//Optional information stored
	private String readGroup = null;
	
	public ReadAlignment(String sequenceName, int first, int last, int readLength, int flags) {
		this.sequenceName = sequenceName;
		this.first = first;
		this.last = last;
		this.readLength = readLength;
		this.setFlags(flags);
		updateAlleleCallsInfo();
	}

	@Override
	public String getSequenceName() {
		return sequenceName;
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
	
	
	
	public int getFlags() {
		return flags;
	}

	public void setFlags(int flags) {
		//TODO: Check consistency
		this.flags = flags;
		if(isReadUnmapped()) {
			alignment = null;
			last = 0;
		}
		alleleCallsUpdated = false;
	}
	
	public short getAlignmentQuality() {
		return alignmentQuality;
	}

	public void setAlignmentQuality(short alignmentQuality) {
		this.alignmentQuality = alignmentQuality;
	}

	public String getMateSequenceName() {
		return mateSequenceName;
	}

	public void setMateSequenceName(String mateSequenceName) {
		this.mateSequenceName = mateSequenceName;
	}

	public int getMateFirst() {
		return mateFirst;
	}

	public void setMateFirst(int mateFirst) {
		this.mateFirst = mateFirst;
	}

	public int getInferredInsertSize() {
		return inferredInsertSize;
	}

	public void setInferredInsertSize(int inferredInsertSize) {
		this.inferredInsertSize = inferredInsertSize;
	}

	public boolean isPaired() {
		return (flags & FLAG_PAIRED)!=0;
	}

	public boolean isProperPair() {
		return (flags & FLAG_PROPER)!=0;
	}
	
	public boolean isReadUnmapped() {
		return (flags & FLAG_READ_UNMAPPED)!=0;
	}
	
	public boolean isMateUnmapped() {
		return (flags & FLAG_MATE_UNMAPPED)!=0;
	}
	
	public boolean isMatePositiveStrand() {
		return (flags & FLAG_MATE_REVERSE_STRAND)==0;
	}
	public boolean isMateNegativeStrand() {
		return (flags & FLAG_MATE_REVERSE_STRAND)!=0;
	}
	
	public void setMateNegativeStrand (boolean mateNegativeStrand) {
		if(mateNegativeStrand) flags = flags | FLAG_MATE_REVERSE_STRAND;
		else flags = flags & ~FLAG_MATE_REVERSE_STRAND;
	}
	public boolean isFirstOfPair () {
		return (flags & FLAG_FIRST_OF_PAIR)!=0;
	}
	public boolean isSecondOfPair () {
		return (flags & FLAG_SECOND_OF_PAIR)!=0;
	}
	public boolean isSecondary () {
		return (flags & FLAG_SECONDARY)!=0;
	}
	public boolean failsQC () {
		return (flags & FLAG_FAILS_QC)!=0;
	}
	public boolean isPCRDuplicate () {
		return (flags & FLAG_PCR_DUP)!=0;
	}
	public boolean isSupplementary () {
		return (flags & FLAG_SUPPLEMENTARY)!=0;
	}
	public boolean isUnique () {
		//TODO: Make rules to be unique part of the constructor
		return (flags & FLAG_MULTIPLE_ALN)==0;
	}
	public boolean isMateDifferentSequence () {
		return (flags & FLAG_MATE_DIFFERENT_SEQUENCE)!=0;
	}
	public boolean isMateSameSequence () {
		return (flags & FLAG_MATE_DIFFERENT_SEQUENCE)==0;
	}
	public String getReadName() {
		if(readName == null) return null;
		return new String (readName);
	}

	public void setReadName(String readName) {
		if(readName == null) this.readName = null;
		this.readName = readName.toCharArray();
	}

	public CharSequence getReadCharacters() {
		return readCharacters;
	}

	public void setReadCharacters(CharSequence readCharacters) {
		if(readCharacters!=null && readLength>0 &&readCharacters.length()!=readLength) throw new IllegalArgumentException("Input Read length: "+readCharacters.length()+" inconsistent with the expected length "+readLength);		 
		this.readCharacters = readCharacters;
		if(readCharacters == null) this.qualityScores = null;
	}
	
	public int getReadLength() {
		return readLength;
	}

	public String getQualityScores() {
		if(qualityScores == null) return null;
		return new String(qualityScores);
	}

	public void setQualityScores(String qualityScores) {
		if(qualityScores==null) this.qualityScores = null;
		this.qualityScores = qualityScores.toCharArray();
		if(readCharacters!=null) fillQualityScores();
	}
	
	public void setBasesToIgnoreCloseToIndel(byte basesToIgnoreCloseToIndel) {
		if(basesToIgnoreCloseToIndel<1) throw new IllegalArgumentException("Bases to ignore close to indel must be at least 1");
		this.basesToIgnoreCloseToIndel = basesToIgnoreCloseToIndel;
		this.alleleCallsUpdated = false;
	}
	

	public void setBasesToIgnore5P(short basesToIgnore5P) {
		if(isNegativeStrand()) this.basesToIgnoreEnd = basesToIgnore5P;
		else this.basesToIgnoreStart = basesToIgnore5P;
		this.alleleCallsUpdated = false;
	}

	public void setBasesToIgnore3P(short basesToIgnore3P) {
		if(isNegativeStrand()) this.basesToIgnoreStart = basesToIgnore3P;
		else this.basesToIgnoreEnd = basesToIgnore3P;
		this.alleleCallsUpdated = false;
	}
	
	public short getBasesToIgnoreStart() {
		return basesToIgnoreStart;
	}

	public void setBasesToIgnoreStart(short basesToIgnoreStart) {
		this.basesToIgnoreStart = basesToIgnoreStart;
		this.alleleCallsUpdated = false;
	}

	public short getBasesToIgnoreEnd() {
		return basesToIgnoreEnd;
	}
	
	public void setBasesToIgnoreEnd(short basesToIgnoreEnd) {
		this.basesToIgnoreEnd = basesToIgnoreEnd;
		this.alleleCallsUpdated = false;
	}
	
	public String getReadGroup() {
		return readGroup;
	}

	public void setReadGroup(String readGroup) {
		this.readGroup = readGroup;
	}
	
	/**
	 * Tells if the alignment is partial (if it has soft clipped basepairs)
	 * @param minClipLength Minimum length of the clip to be considered 
	 * @return boolean true if the alignment starts or ends with softclipped basepairs 
	 */
	public boolean isPartialAlignment(int minClipLength) {
	
		if(alignment == null || alignment.length<2) return false;
		int s1 = 0;
		if(getOperator(alignment[0])==ALIGNMENT_SKIPFROMREAD) s1 = getOperationLength(alignment[0]);
		int s2 = 0;
		int l = alignment.length-1;
		if(getOperator(alignment[l])==ALIGNMENT_SKIPFROMREAD) s2 = getOperationLength(alignment[l]);
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
	 * @return 0-based position of the read aligning with the given position. -1 if the reference position does not align with the read
	 */
	public int getReadPosition (int referencePos) {
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
					return currentReadPos + referencePos-currentRefPos; 
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
	 * @param readPos 0-based position of the read
	 * @return int 1-based reference coordinate aligning with the given position. -1 if the read position does not align with the reference
	 */
	public int getReferencePosition (int readPos) {
		failIfReadUnmappedOrInconsistentAlignment();
		int currentRefPos = first;
		int currentReadPos = 0;
		for(int i=0;i<alignment.length;i++) {
			int length = getOperationLength(alignment[i]);
			boolean cRef = consumesReferenceBases(alignment[i]);
			boolean cRead = consumesReadBases(alignment[i]);
			if(cRef && cRead) {
				if(readPos<currentReadPos) return -1;
				else if(currentReadPos+length>readPos) {
					return currentRefPos + readPos-currentReadPos; 
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
	public CharSequence getAlleleCall (int referencePos) {
		if(readCharacters ==null) return null;
		int readPos = getReadPosition(referencePos);
		if(readPos<0) return null;
		updateAlleleCallsInfo();
		//if(readPos>=alleleCallLength.length) System.err.println("Invalid read position: "+readPos+" Read: "+getReadName()+" Located at "+sequenceName+":"+first+" last:"+last+". Length: "+readCharacters.length()+" CIGAR: "+getCigarString()+" expected length: "+expectedReadLength+" reference pos: "+referencePos);
		//if(readPos>=alleleCallLength.length) System.err.println("Alignment: "+alignment[0]+" "+alignment[1]+" "+alignment[2]);
		int length = alleleCallLength[readPos];
		if(length == 0) return null;
		//if(referencePos==-1) System.out.println("ReadAlignment. Read id: "+getReadName()+". Sequence: "+readCharacters.toString()+". readpos: "+readPos+". end: "+(readPos+length)+" subseq: "+readCharacters.subSequence(readPos, readPos+length)+". length subseq: "+readCharacters.subSequence(readPos, readPos+length).length());
		return readCharacters.subSequence(readPos, readPos+length);
	}
	public CharSequence getAlleleCall (int referenceFirst, int referenceLast) {
		if(readCharacters == null) return null;
		updateAlleleCallsInfo();
		int readFirst = getReadPosition(referenceFirst);
		int readLast = getReadPosition(referenceLast);
		if(readFirst<0 || readLast<0 || readLast < readFirst) return null;
		if(withinIgnoreRegions(readFirst, readLast)) return null;
		return readCharacters.subSequence(readFirst, readLast+1);
	}
	public char getBaseQualityScore (int referencePos) {
		if(qualityScores == null) return 33;
		int readPos = getReadPosition(referencePos);
		if(readPos<0) return 33;
		return qualityScores[readPos];
	}
	
	public String getBaseQualityScores (int referenceFirst, int referenceLast) {
		if(qualityScores == null) return null;
		int readFirst = getReadPosition(referenceFirst);
		int readLast = getReadPosition(referenceLast);
		if(readFirst<0 || readLast<0  || readLast < readFirst) return null;
		if(withinIgnoreRegions(readFirst, readLast)) return null;
		return (new String (qualityScores)).substring(readFirst,readLast+1);
	}
	private boolean withinIgnoreRegions (int readFirst, int readLast) {
		return readFirst<basesToIgnoreStart || readLength - readLast <= basesToIgnoreEnd;
	}
	public Map<Integer,GenomicVariant> getIndelCalls () {
		failIfReadUnmappedOrInconsistentAlignment();
		updateAlleleCallsInfo();
		return indelCalls;
	}
	public GenomicVariant getIndelCall(int referencePos) {
		failIfReadUnmappedOrInconsistentAlignment();
		updateAlleleCallsInfo();
		if(indelCalls == null) return null;
		return indelCalls.get(referencePos);
	}
	
	public void collapseIndelEvents() {
		failIfReadUnmappedOrInconsistentAlignment();
		List<Integer> newAlignmentList = new ArrayList<Integer>();
		int i=0;
		int n = alignment.length;
		//String cigar = getCigarString();
		if(n<5) return;
		int changeNextEvent = 0;
		while(i<n) {
			int length = getOperationLength(alignment[i]);
			byte operator = getOperator(alignment[i]);
			//if(aln.getAlignmentStart()==77467) System.out.println("WARN. Processing read: "+aln.getReadName()+" at "+aln.getReferenceName()+":"+aln.getAlignmentStart()+". Current CIGAR: "+aln.getCigarString()+" next elem: "+length+operator.name()+" change next event: "+changeNextEvent+" is Indel: "+isIndel(elem));
			if(isIndel(operator)) {
				for(int j=i+1;j<n;j++) {
					int l2 = getOperationLength(alignment[j]);
					byte op2 = getOperator(alignment[j]);
					
					if(isIndel(op2)) {
						if(operator == op2) {
							length+=l2;
						} else {
							int newLength = length - l2;
							changeNextEvent+=Math.min(length, l2);
							if(newLength<0) {
								//Switch operator to make length positive
								if(operator == ALIGNMENT_INSERTION) operator = ALIGNMENT_DELETION;
								else operator = ALIGNMENT_INSERTION;
								newLength = -newLength;
							}
							length = newLength;
						}
						if(j==n-1) {
							if(changeNextEvent == 0) System.err.println("WARN. A modified CIGAR will be finished with an indel for the alignment of read: "+getReadName()+" at "+getSequenceName()+":"+getFirst()+". Current CIGAR: "+getCigarString());
							i=n;
						}
					} else if(j==n-1) {
						//Last event. Do not try to merge it
						i=j;
						break;
					} else {
						//Look ahead for next indel
						int l3 = getOperationLength(alignment[j+1]);
						byte op3 = getOperator(alignment[j+1]);
						int length2 = l2+changeNextEvent;
						if (isMatchMismatch(op2) && isIndel(op3) && (length2<=5 || length2<length+l3)) {
							//Short M event in the middle of two indels. Collapse with event after
							changeNextEvent+=l2;
							
						} else {
							i=j;
							break;
						}
					}
				}
			} else {
				if(changeNextEvent>0) {
					length+=changeNextEvent;
					changeNextEvent = 0;
				}
				i++;
			}
			if(length>0) {
				newAlignmentList.add(getAlnValue(length, operator));
			}
		}
		if(changeNextEvent>0) {
			//Soft clip if the read ends with an indel
			newAlignmentList.add(getAlnValue(changeNextEvent, ALIGNMENT_SKIPFROMREAD));
			System.err.println("WARN. A modified CIGAR will be finished with soft clip for the alignment of read at "+getSequenceName()+":"+first);
		}
		
		alignment = NumberArrays.toIntArray(newAlignmentList);
		//if(!cigar.equals(getCigarString())) System.out.println("WARN. Changing CIGAR of read: "+getReadName()+" at "+getSequenceName()+":"+first+". Current CIGAR: "+cigar+" new CIGAR: "+getCigarString());
		alleleCallsUpdated = false;
	}
	
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
	//Alignment encoding / decoding methods
	private int getOperationLength(int alnValue) {
		return alnValue/8;
	}

	private byte getOperator(int alnValue) {
		return (byte)(alnValue & 0x7);
	}
	
	private boolean consumesReadBases (int alnValue) {
		return (alnValue & 0x2) != 0;
	}
	
	private boolean consumesReferenceBases (int alnValue) {
		return (alnValue & 0x1) != 0;
	}
	
	private int getAlnValue (int operationLength , byte operator) {
		return 8*operationLength + operator;
	}
	
	private boolean isIndel(byte operator) {
		return operator == ALIGNMENT_DELETION || operator == ALIGNMENT_INSERTION;
	}
	
	private boolean isMatchMismatch(byte operator) {
		return operator == ALIGNMENT_MATCH || operator == ALIGNMENT_MISMATCH;
	}
	
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
		alignment = NumberArrays.toIntArray(alignmentList);
		if(expectedReadLength != this.readLength) {
			System.out.println("WARN. New CIGAR changes read length for read at "+sequenceName+":"+first);
			this.readLength = expectedReadLength;
		}
		if(expectedEnd -1 != this.last) {
			System.out.println("WARN. New CIGAR changes last alignment position for read at "+sequenceName+":"+first);
			this.last = expectedEnd - 1;
		}
		if(first == posPrint && alignment.length >=3) System.out.println("Alignment codes: "+alignment[0]+" "+alignment[1]+" "+alignment[2]+" Expected read length: "+expectedReadLength);
		alleleCallsUpdated = false;
	}
	
	private void fillQualityScores() {
		int l = readCharacters.length(); 
		if(l>qualityScores.length) {
			char [] fakeLowScores = new char[l];
			for(int i=0;i<l;i++) {
				if(i<qualityScores.length) fakeLowScores[i] = qualityScores[i];
				else fakeLowScores[i] = '$';
			}
			qualityScores = fakeLowScores;
		}
	}

}
