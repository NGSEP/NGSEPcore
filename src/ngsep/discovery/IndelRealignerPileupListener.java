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
package ngsep.discovery;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.math.NumberArrays;
import ngsep.sequences.DNASequence;
import ngsep.sequences.HammingSequenceDistanceMeasure;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.GenomicVariant;


public class IndelRealignerPileupListener implements PileupListener {

	private static final int DEF_REGION_BOUNDARY = 100;
	private GenomicRegionSortedCollection<? extends GenomicVariant> inputVariants = null;
	private ReferenceGenome genome;
	private int minBPForGoodRefAln = 5;
	private int maxBPRealignmentEnd = 50;
	private HammingSequenceDistanceMeasure hammingMeasure = new HammingSequenceDistanceMeasure();
	
	
	private List<? extends GenomicVariant> seqInputVariants;
	private int idxNextVariant = 0;
	
	//DEBUG
	private int posPrint = -1;

	public GenomicRegionSortedCollection<? extends GenomicRegion> getInputVariants() {
		return inputVariants;
	}

	public void setInputVariants(GenomicRegionSortedCollection<? extends GenomicVariant> inputVariants) {
		this.inputVariants = inputVariants;
	}
	
	
	
	public ReferenceGenome getGenome() {
		return genome;
	}

	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}

	@Override
	public void onPileup(PileupRecord pileup) {
		int currentPos = pileup.getPosition();
		int referenceSpan = 1;
		int predictedEventEnd=currentPos;
		GenomicVariant var = intersectWithVariants (pileup);
		if(var!=null) {
			if(var.getFirst()==pileup.getPosition()) {
				if(var.getType()==GenomicVariant.TYPE_STR) pileup.setSTR(true);
				referenceSpan = var.getLast() - var.getFirst() + 1;
				predictedEventEnd = var.getLast();
				if(currentPos==posPrint)System.out.println("IndelRealigner. Var span: "+referenceSpan+". Predicted end: "+predictedEventEnd+" Var type "+var.getType());
			} else {
				pileup.setEmbedded(true);
			}
		}
		List<ReadAlignment> alignments = pileup.getAlignments();
		if(var==null) {
			//Look for the start of a new indel event from this position
			int maxIndelLength = 0;
			int maxIndelSpan = 0;
			for(ReadAlignment aln:alignments) {
				GenomicVariant indel = aln.getIndelCall(currentPos);
				if(currentPos==posPrint) System.out.println("Read name: "+aln.getReadName()+". Alignment start: "+aln.getFirst()+" CIGAR: "+aln.getCigarString()+" Indels: "+aln.getIndelCalls());
				if(indel!=null) {
					int indelSpan = indel.getLast() - indel.getFirst() + 1;
					int indelLength = indel.length();
					if(currentPos==posPrint) System.out.println("Read name: "+aln.getReadName()+". Alignment start: "+aln.getFirst()+" CIGAR: "+aln.getCigarString()+" Indel length: "+indelLength);
					if(indelLength>maxIndelLength) maxIndelLength = indelLength;
					if(indelSpan > maxIndelSpan) maxIndelSpan = indelSpan;
				}
			}
			if(maxIndelLength>0) predictedEventEnd = currentPos+Math.max(maxIndelLength, maxIndelSpan)+1;
			if(currentPos==posPrint)System.out.println("IndelRealigner. Max indel length: "+maxIndelLength+". Max indel Span: "+maxIndelSpan+" predicted end: "+predictedEventEnd+" time: "+System.currentTimeMillis());
		}
		if(predictedEventEnd>currentPos) {
			int conciliatedSpan = conciliateIndels(pileup,alignments, predictedEventEnd,var);
			if(conciliatedSpan > 0) referenceSpan = conciliatedSpan;
			if(currentPos==posPrint)System.out.println("IndelRealigner. New reference span: "+referenceSpan+". STR: "+pileup.isSTR()+" time: "+System.currentTimeMillis());
		}
		pileup.setReferenceSpan(referenceSpan);
		if(currentPos==posPrint)System.out.println("IndelRealigner. Processed pileup at  "+pileup.getSequenceName()+":"+currentPos+" Reference span: "+referenceSpan+" STR: "+pileup.isSTR()+". time: "+System.currentTimeMillis());
	}

	@Override
	public void onSequenceStart(QualifiedSequence sequence) {
		if(inputVariants!=null) {
			seqInputVariants = inputVariants.getSequenceRegions(sequence.getName()).asList();
			idxNextVariant = 0;
		}
	}

	@Override
	public void onSequenceEnd(QualifiedSequence sequence) {

	}

	private GenomicVariant intersectWithVariants(PileupRecord pileup) {
		int currentPos = pileup.getPosition();
		if(inputVariants!=null) {
			if(currentPos==posPrint)System.out.println("Idx next variant. "+idxNextVariant+" inputVariants: "+inputVariants.size());
			while (idxNextVariant<seqInputVariants.size()) {
				GenomicVariant var = seqInputVariants.get(idxNextVariant);
				int first = var.getFirst();
				int last = var.getLast();
				if(currentPos==posPrint)System.out.println("Next variant. "+var.getSequenceName()+":"+first+"-"+last);
				if (currentPos<first) break;
				else if(first<=currentPos && currentPos<=last) return var;
				//System.out.println("Passed STR at "+currentReferenceName+":"+first+"-"+last+" current reference pos: "+currentReferencePos);
				idxNextVariant++;
			}
		}
		return null;
	}
	
	/**
	 * @param pileup record with alignments
	 * @param eventEnd predicted event end
	 * @param varG variant to genotype
	 * @return int The reference span if indels are called (zero otherwise).
	 */
	private int conciliateIndels(PileupRecord pileup, List<ReadAlignment> alignments, int eventEnd, GenomicVariant varG) {
		int answer = 0;
		int currentPos = pileup.getPosition();
		boolean fixedEvent = (varG!=null);
		
		Set<Integer> lengths = new TreeSet<Integer>();
		List<ReadAlignment> indelAlns = new ArrayList<ReadAlignment>();
		//Vote for possible indel starts
		int [] votes = new int [eventEnd-currentPos+1];
		
		int maxLength = analyzeIndels(alignments,currentPos,eventEnd,varG, lengths, indelAlns, votes);
		
		if(currentPos==posPrint)System.out.println("ConciliateIndels. Current pos: "+currentPos+" Max length: "+maxLength+" lengths: "+ lengths.size());
		
		if(lengths.size()==0) return answer;
		
		
		int maxI = 0;
		if(!fixedEvent) {
			maxI = NumberArrays.getIndexMaximum(votes);
			if(lengths.size()>1) {
				//Look for new tandem repeats
				if(currentPos==posPrint) System.out.println("ConciliateIndels. Checking tandem repeats.");
				int newSpan = lookForNewSTR(pileup, indelAlns, maxLength);
				if(newSpan>1) {
					//eventEnd = currentPos+newSpan-1;
					maxI = 0;
					answer = newSpan;
					eventEnd = currentPos+answer-1;
					fixedEvent = true;
					pileup.setSTR(true);
					pileup.setNewSTR(true);
					if(currentPos==posPrint) System.out.println("New tandem repeat identified from reads at "+pileup.getSequenceName()+":"+currentPos+" New max span: "+answer);
				}
			}
		}
		if(currentPos==posPrint)System.out.println("ConciliateIndels. Current pos: "+currentPos+" Start max votes: "+(currentPos+maxI)+" total alns: "+alignments.size()+" indelAlns: "+ indelAlns.size()+" event end: "+eventEnd+" max length: "+maxLength+" maxI: "+maxI);
		
		//Move indel starts to the position with the largest number of votes
		int newPredictedEventEnd = moveIndelStarts(indelAlns, currentPos, eventEnd, maxLength, maxI);
		if(maxI > 0) return answer;
		if(!fixedEvent && newPredictedEventEnd!=eventEnd) {
			//Update the predicted end the new indel event after realignment
			eventEnd = newPredictedEventEnd;
			answer = eventEnd-currentPos+1;
			if(currentPos==posPrint) System.out.println("ConciliateIndels. Predicted span for new indel: "+answer);	
		}
		
		//Try to realign ends of alignments falling within the event
		processEndsOfAlignments(alignments,pileup.getSequenceName(), currentPos, eventEnd);
		return answer;
	}

	/**
	 * Calculates max length, set of indel lengths, alignments with indels and votes for indel start position
	 * @param alignments to analyze
	 * @param eventStart First genomic coordinate of the event
	 * @param eventEnd End genomic coordinate of the event
	 * @param varG Variant to genotype
	 * @param lengths Output set of event lengths
	 * @param indelAlns Alignments with indels within the event
	 * @param votes Array of votes for start positions. Indexes are relative to eventStart
	 * @return int Maximum length of an indel event
	 */
	private int analyzeIndels(List<ReadAlignment> alignments, int eventStart, int eventEnd, GenomicVariant varG, Set<Integer> lengths, List<ReadAlignment> indelAlns, int[] votes) {	
		
		int firstOverlap = eventStart;
		int lastOverlap = eventEnd;
		/*if(varG!=null && varG.getLast()>varG.getFirst()) {
			//Small increase in search range to genotype indels and STRs
			firstOverlap-=3;
			lastOverlap+=3;
		}*/
		int maxLength = 0;
		lengths.clear();
		indelAlns.clear();
		if(votes!=null)Arrays.fill(votes,0);
		for(ReadAlignment aln:alignments) {
			boolean indelFound = false;
			Map<Integer,GenomicVariant> indels = aln.getIndelCalls();
			if(indels!=null) {
				for(int start:indels.keySet()) {
					GenomicVariant indel = indels.get(start);
					if(eventStart==posPrint) System.out.println("Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+" CIGAR: "+aln.getCigarString()+" Next indel start "+start+" event limits "+eventStart+"-"+eventEnd);
					
					if(indel.getLast() >= firstOverlap && start <=lastOverlap) {
						indelFound = true;
						int length = indel.length();
						lengths.add(length);
						if(maxLength<length) maxLength = length;
						int i = start-eventStart;
						if(eventStart==posPrint) System.out.println("Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+" CIGAR: "+aln.getCigarString()+" Insertion start: "+start+" vote: "+i);
						if(votes!=null && i>=0 && i<votes.length)votes[i]++;
						break;
					}
				}
			}
			if(indelFound) indelAlns.add(aln);
		}
		return maxLength;
	}
	/**
	 * Realigns the indels of the given alignments located within 
	 * @param alignments
	 * @param currentPos
	 * @param eventEnd
	 * @param offset
	 * @return int New predicted end of the event
	 */
	public int moveIndelStarts(List<ReadAlignment> alignments, int first, int last, int maxLength, int offset) {
		int firstOverlap = first;
		int lastOverlap = last;
		/*if(varG!=null && varG.getLast()>varG.getFirst()) {
			//Small increase in search range to genotype indels and STRs
			firstOverlap-=3;
			lastOverlap+=3;
		}*/
		int answer = first+1;
		for(ReadAlignment aln:alignments) {
			Map<Integer,GenomicVariant> indels = aln.getIndelCalls();
			if(indels!=null) {
				for(int start:indels.keySet()) {
					GenomicVariant indel = indels.get(start);
					//if(eventStart==posPrint) System.out.println("Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+" CIGAR: "+aln.getCigarString()+" Next indel start "+start+" event limits "+eventStart+"-"+eventEnd);
					
					if(indel.getLast() >= firstOverlap && start <=lastOverlap) {
					//if(start >= currentPos && start <=eventEnd) {		
						//if(aln.getFirst()==1291016) System.out.println("Trying to move indel start for alignment of read "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" indel reference pos "+start+" offset: "+offset);
						boolean moved = aln.moveIndelStart(start,first+offset);
						if(first==posPrint && moved == false) System.err.println("WARN: Failed attempt to move indel start for alignment of read "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" indel reference pos "+start+" current pileup pos: "+first+" new indel start "+(first+offset)+" read pos: "+aln.getAlignedReadPosition(first));
						break;
					}
				}
			}
			indels = aln.getIndelCalls();
			if(indels!=null) {
				int alnRefLast = first;
				for(int start:indels.keySet()) {
					//This takes into account several close indel events within the same alignment
					if(start >= first && start <=alnRefLast+maxLength) {
						alnRefLast = indels.get(start).getLast();
					}
				}
				if(alnRefLast>answer) answer = alnRefLast;
			}
			if(first==posPrint) System.err.println("Ref end after trying to  move indel start for alignment of read "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" is "+answer);
		}
		return answer;
	}
	
	private int lookForNewSTR(PileupRecord pileup, List<ReadAlignment> alns, int maxLength) {
		int newSpan = 0;
		int currentPos = pileup.getPosition();
		//Look first in the reference
		if(alns.size()==0) return 0;
		CharSequence seq = genome.getReference(pileup.getSequenceName(), currentPos+1, alns.get(alns.size()-1).getLast());
		int lengthRef = 0;
		if(seq!=null) {
			String reference = seq.toString().toUpperCase();
			lengthRef = checkTandemRepeat(reference);
		}
		if(lengthRef>0) {
			//System.out.println("IndelRealigner. Found new STR from reference at "+pileup.getSequenceName()+":"+currentPos+" sequence: "+reference+" length: "+lengthRef);
			return lengthRef+2;
		}
		for(ReadAlignment aln:alns) {
			GenomicVariant indel = aln.getIndelCall(currentPos);
			if(indel!=null && indel.length() == maxLength) {
				int indelSpan = indel.getLast() - indel.getFirst() + 1;
				int indelLength = indel.length();
				int lengthTR = checkTandemRepeat(aln, currentPos);
				if(lengthTR > 0) {
					//SAMRecord alnRecord = aln.getSAMRecord();
					//System.out.println("IndelRealigner. Found new STR candidate from reads at "+pileup.getSequenceName()+":"+currentPos+" read: "+alnRecord.getReadString()+" readPos: "+aln.getReadPosition(currentPos)+" CIGAR: "+alnRecord.getCigarString()+" indel length: "+length+" reference span: "+span+" length tandem repeat: "+lengthTR);
					//For new STRs not present in the reference (hence insertions), although the span should be 2, it is better to leave a large span to accommodate alternative alignments
					if(indelLength>=indelSpan) newSpan = indelLength+2;
					else {
						newSpan = lengthTR+indelSpan;
					}
					break;
				}
			}
		}
		return newSpan;
	}

	private int checkTandemRepeat(ReadAlignment aln, int currentPos) {
		int readFirst = aln.getAlignedReadPosition(currentPos);
		if(readFirst < 0) return 0;
		String seq = aln.getReadCharacters().toString();
		seq = seq.substring(readFirst+1).toUpperCase();
		return checkTandemRepeat(seq);
	}

	private int checkTandemRepeat(String seq) {
		int length = checkMonoNucleotide(seq);
		if(length==0) length = checkDinucleotide(seq);
		return length;
	}

	private int checkMonoNucleotide(String seq) {
		int [] counts = new int [4];
		int minLength = 5;
		int i=0;
		while(i<seq.length() && i<minLength) {
			int j = DNASequence.BASES_STRING.indexOf(seq.charAt(i));
			if(j>=0) counts[j]++;
			i++;
		}
		int baseIdx = -1;
		for(int j=0;j<counts.length;j++) {
			if(counts[j]>=i-1) {
				baseIdx = j;
				break;
			}
		}
		if(baseIdx == -1 || i<minLength) return 0;
		while(i<seq.length() && counts[baseIdx]>=i-1) {
			int j = DNASequence.BASES_STRING.indexOf(seq.charAt(i));
			if(j>=0) counts[j]++;
			i++;
		}
		i--;
		//System.out.println("Tandem repeat from sequence: "+seq+" I: "+i);
		if(DNASequence.BASES_STRING.indexOf(seq.charAt(i-1))!=baseIdx) return i-1;
		return i;
	}

	private int checkDinucleotide(String seq) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	
	
	private void processEndsOfAlignments(List<ReadAlignment> alignments, String sequenceName, int eventFirst, int eventLast) {
		//Calculate reference alleles before and after
		CharSequence seqBefore = genome.getReference(sequenceName, eventFirst-DEF_REGION_BOUNDARY, eventFirst);
		CharSequence seqAfter = genome.getReference(sequenceName, eventLast, eventLast+DEF_REGION_BOUNDARY);
		CharSequence seqWithin = null;
		if(eventFirst!=eventLast-1) seqWithin = genome.getReference(sequenceName, eventFirst+1, eventLast-1);
		if(seqWithin!=null) seqWithin = seqWithin.toString().toUpperCase();
		String refAlleleBefore = seqBefore!=null?seqBefore.toString().toUpperCase():null;
		String refAlleleAfter = seqAfter!=null?seqAfter.toString().toUpperCase():null;
		if(refAlleleBefore!=null && seqWithin!=null) refAlleleBefore+=seqWithin;
		if(refAlleleAfter!=null && seqWithin!=null) refAlleleAfter=seqWithin+refAlleleAfter;
		
		//Calculate alternative alleles before and after
		String insertedConsensusSequence = calculateInsertedConsensusSequence (alignments, eventFirst);
		if(eventFirst==posPrint) System.out.println("Inserted consensus sequence: "+insertedConsensusSequence);
		
		String altAlleleBefore = seqBefore!=null?seqBefore.toString().toUpperCase():null;
		String altAlleleAfter = seqAfter!=null?seqAfter.toString().toUpperCase():null;
		int offset = 0;
		if(insertedConsensusSequence!=null) {
			offset = insertedConsensusSequence.length();
			if(altAlleleBefore!=null) altAlleleBefore+=insertedConsensusSequence;
			if(altAlleleBefore!=null && seqWithin!=null) altAlleleBefore+=seqWithin;
			if(altAlleleAfter!=null && seqWithin!=null) altAlleleAfter=seqWithin+altAlleleAfter;
			if(altAlleleAfter!=null) altAlleleAfter=insertedConsensusSequence+altAlleleAfter;
		} else {
			int deletionLength = calculateDeletionConsensusLength (alignments, eventFirst);
			if(deletionLength>eventLast-eventFirst-1) deletionLength=eventLast-eventFirst-1;
			offset = -deletionLength;
			if(altAlleleBefore!=null && seqWithin!=null) {
				if(deletionLength==0) altAlleleBefore+=seqWithin;
				else if(deletionLength<seqWithin.length()) {
					String remainder = seqWithin.toString().substring(deletionLength);
					altAlleleBefore+=remainder;
				}
			}
			if(altAlleleAfter!=null && seqWithin!=null) {
				if(deletionLength==0) altAlleleAfter=seqWithin+altAlleleAfter;
				else if(deletionLength<seqWithin.length()) {
					String remainder = seqWithin.toString().substring(deletionLength);
					altAlleleAfter=remainder+altAlleleAfter;
				}
			}
		}
		if(eventFirst==posPrint) System.out.println("Alt allele after: "+altAlleleAfter);
		for(ReadAlignment aln:alignments) {
			int alnFirst = aln.getFirst();
			int alnLast = aln.getLast();
			
			String cigarStr = aln.getCigarString();
			boolean hasIndelCallsBefore = aln.hasIndelCalls(alnFirst,eventFirst-1);
			boolean hasIndelCallsAfter = aln.hasIndelCalls(eventLast+1,alnLast);
			
			int bpForGoodRefAln = Math.max(offset, minBPForGoodRefAln);
			boolean trimStart=eventFirst-alnFirst<bpForGoodRefAln && !hasIndelCallsBefore;
			int readPosAfter = aln.getAlignedReadPosition(eventLast);
			if(eventFirst==posPrint) System.out.println("IndelRealigner. realignStarts. Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+". CIGAR: "+aln.getCigarString()+". Event limits: "+eventFirst+"-"+eventLast+" readPosAfter: "+readPosAfter+" offset:"+offset );
			if(!hasIndelCallsBefore && refAlleleBefore!=null && altAlleleBefore!=null && readPosAfter>=bpForGoodRefAln && readPosAfter-offset<=maxBPRealignmentEnd && readPosAfter<refAlleleBefore.length() && readPosAfter<altAlleleBefore.length() && aln.getIndelCall(eventFirst)==null) {
				CharSequence readPrefix = aln.getReadCharacters().subSequence(0, readPosAfter);
				CharSequence refSuffix = refAlleleBefore.substring(refAlleleBefore.length()-readPosAfter);
				if(eventFirst==posPrint) System.out.println(readPrefix);
				if(eventFirst==posPrint) System.out.println(refSuffix);
				double referenceDistance = hammingMeasure.calculateDistance(refSuffix, readPrefix);
				CharSequence altSuffix = altAlleleBefore.substring(altAlleleBefore.length()-readPosAfter);
				if(eventFirst==posPrint) System.out.println(altSuffix);
				double alternativeDistance = hammingMeasure.calculateDistance(altSuffix, readPrefix);
				int newAlnFirst = eventLast-readPosAfter+1+offset;
				int firstMatchLength = eventFirst-newAlnFirst+1;
				if(eventFirst == posPrint) System.out.println("IndelRealigner. realignStarts. Reference distance: "+referenceDistance+" alt distance: "+alternativeDistance+" new aln first: "+newAlnFirst+"");
				if(alternativeDistance<referenceDistance && alternativeDistance<3 && firstMatchLength>=minBPForGoodRefAln) {
					aln.realignStart(newAlnFirst,firstMatchLength,eventLast, readPosAfter);
					trimStart = false;
					if(eventFirst == posPrint) System.out.println("IndelRealigner. realignEnds. Realigned start of alignment with original coordinates: "+alnFirst+"-"+alnLast+" old CIGAR: "+cigarStr+" new start: "+aln.getFirst()+" new CIGAR: "+aln.getCigarString());
				}
			}	
			if(trimStart) {
				int ignoreBP = eventLast-alnFirst+1;
				ignoreBP+=aln.getSoftClipStart();
				byte bpToIgnoreStart = (byte)Math.max(aln.getBasesToIgnoreStart(), ignoreBP);
				
				//if(bpToIgnoreStart>10)System.err.println("WARN: Ignoring "+bpToIgnoreStart+" base pairs at the start of alignment of read "+aln.getSAMRecord().getReadName()+" at "+aln.getReferenceName()+":"+aln.getAlignmentStart()+ " Current CIGAR: "+aln.getSAMRecord().getCigarString()+" indel alns: "+numIndelAlns+" non indel alns: "+alns.size()+" event first: "+eventFirst+" event last: "+eventLast);
				aln.setBasesToIgnoreStart(bpToIgnoreStart);
				if(eventFirst == posPrint) System.out.println("IndelRealigner. Trimmed "+bpToIgnoreStart+" at the start of alignment with coordinates: "+alnFirst+"-"+alnLast);
			}
			
			boolean trimEnd = alnLast-eventLast<bpForGoodRefAln && !hasIndelCallsAfter;
			int readPosBefore = aln.getAlignedReadPosition(eventFirst);
			//
			int readSuffixLength = 0;
			if(readPosBefore>=0) { 
				readSuffixLength = aln.getReadLength()-readPosBefore-1;
			}
			if(eventFirst==posPrint) System.out.println("IndelRealigner. realignEnd. Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+". CIGAR: "+aln.getCigarString()+" Event first: "+eventFirst+" readPosBefore: "+readPosBefore+" suffix length: "+readSuffixLength);
			if(!hasIndelCallsAfter && refAlleleAfter!=null && altAlleleAfter!=null && readSuffixLength>=bpForGoodRefAln && readSuffixLength-offset<=maxBPRealignmentEnd && readSuffixLength<refAlleleAfter.length() && readSuffixLength<altAlleleAfter.length() && (aln.getIndelCall(eventFirst)==null || readPosAfter<0)) {
				CharSequence readSuffix = aln.getReadCharacters().subSequence(readPosBefore+1,aln.getReadLength());
				CharSequence refPrefix = refAlleleAfter.substring(0, readSuffixLength);
				if(eventFirst==posPrint) System.out.println(readSuffix);
				if(eventFirst==posPrint) System.out.println(refPrefix);
				double referenceDistance = hammingMeasure.calculateDistance(refPrefix, readSuffix);
				CharSequence altPrefix = altAlleleAfter.substring(0, readSuffixLength);
				if(eventFirst==posPrint) System.out.println(altPrefix);
				double alternativeDistance = hammingMeasure.calculateDistance(altPrefix, readSuffix);
				int finalMatchLength = readSuffixLength;
				if(offset>0)finalMatchLength-=offset;
				int newEventLast = eventFirst+1;
				if(offset<0) newEventLast-=offset;
				
				
				if(eventFirst == posPrint) System.out.println("IndelRealigner. realignEnds. Reference distance: "+referenceDistance+" alt distance: "+alternativeDistance+" event coords: "+eventFirst+"-"+eventLast+" offset: "+offset+" final match length: "+finalMatchLength);
				if(alternativeDistance<referenceDistance && alternativeDistance<3 && finalMatchLength>=minBPForGoodRefAln) {
					aln.realignEnd(eventFirst, readPosBefore, newEventLast, finalMatchLength);
					trimEnd = false;
					if(eventFirst == posPrint) System.out.println("IndelRealigner. realignEnds. Realigned end of alignment with original coordinates: "+alnFirst+"-"+alnLast+" old CIGAR: "+cigarStr+" new end: "+aln.getLast()+" new CIGAR: "+aln.getCigarString());
				}
				
				
			}
			if(trimEnd) {
				int ignoreBP = alnLast-eventFirst+1;
				ignoreBP+=aln.getSoftClipEnd();
				byte bpToIgnoreEnd = (byte)Math.max(aln.getBasesToIgnoreEnd(), ignoreBP);
				//if(bpToIgnoreEnd>10)System.err.println("WARN: Ignoring "+bpToIgnoreEnd+" base pairs at the end of alignment of read "+aln.getSAMRecord().getReadName()+" at "+aln.getReferenceName()+":"+aln.getAlignmentStart()+ " Current CIGAR: "+aln.getSAMRecord().getCigarString()+" indel alns: "+numIndelAlns+" non indel alns: "+alns.size()+" event first: "+eventFirst+" event last: "+eventLast);
				aln.setBasesToIgnoreEnd(bpToIgnoreEnd);
				if(eventFirst == posPrint) System.out.println("IndelRealigner. Trimmed "+bpToIgnoreEnd+" at the end of alignment with coordinates: "+alnFirst+"-"+alnLast);
			}
		}
	}

	private String calculateInsertedConsensusSequence(List<ReadAlignment> alignments, int eventFirst) {
		Map<Integer, List<String>> insertionAllelesByLength = new HashMap<>();
		for(ReadAlignment aln:alignments) {
			GenomicVariant call = aln.getIndelCall(eventFirst);
			if(call==null) continue;
			CharSequence allele = aln.getAlleleCall(eventFirst);
			if(allele==null) continue;
			allele = allele.subSequence(1, allele.length()-1);
			if(allele.length()==0) continue;
			//if(eventFirst==posPrint) System.out.println("Calculating insertion consensus allele. Next insertion allele: "+allele+ " read: "+aln.getReadName());
			List<String> insertionAllelesLength = insertionAllelesByLength.get(allele.length());
			if(insertionAllelesLength==null) {
				insertionAllelesLength = new ArrayList<>();
				insertionAllelesByLength.put(allele.length(), insertionAllelesLength);
			}
			insertionAllelesLength.add(allele.toString());
		}
		int numAllelesConsensus = 0;
		List<String> allelesConsensus = null;
		for(List<String> alleles:insertionAllelesByLength.values()) {
			if(alleles.size()>numAllelesConsensus) {
				numAllelesConsensus = alleles.size();
				allelesConsensus = alleles;
			}
		}
		if(allelesConsensus==null) return null;
		return HammingSequenceDistanceMeasure.makeHammingConsensus(allelesConsensus);
	}

	private int calculateDeletionConsensusLength(List<ReadAlignment> alignments, int eventFirst) {
		Map<Integer, Integer> counts = new HashMap<>();
		for(ReadAlignment aln:alignments) {
			GenomicVariant call = aln.getIndelCall(eventFirst);
			if(call==null) continue;
			int innerLength = call.getLast()-call.getFirst()-1;
			Integer count = counts.get(innerLength);
			if(count == null) count = 0;
			count++;
			counts.put(innerLength, count);
		}
		int max = 0;
		int answer = 0;
		for(int length:counts.keySet()) {
			int count = counts.get(length);
			if(max<count) {
				answer = length;
				max = count;
			}
		}
		return answer;
	}
}
