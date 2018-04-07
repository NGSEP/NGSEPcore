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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionSortedCollection;
import ngsep.genome.ReferenceGenome;
import ngsep.math.NumberArrays;
import ngsep.sequences.AbstractLimitedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.variants.GenomicVariant;


public class IndelRealignerPileupListener implements PileupListener {

	private GenomicRegionSortedCollection<? extends GenomicVariant> inputVariants = null;
	private ReferenceGenome genome;
	private int bpForGoodRefAln = 5;
	
	
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
			int conciliatedSpan = conciliateIndels(pileup,predictedEventEnd,var);
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
	private int conciliateIndels(PileupRecord pileup, int eventEnd, GenomicVariant varG) {
		int answer = 0;
		int currentPos = pileup.getPosition();
		boolean fixedEvent = (varG!=null);
		
		Set<Integer> lengths = new TreeSet<Integer>();
		List<ReadAlignment> alignments = pileup.getAlignments();
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
		if(currentPos==posPrint)System.out.println("ConciliateIndels. Current pos: "+currentPos+" Start max votes: "+(currentPos+maxI)+" indelAlns: "+ indelAlns.size());
		
		//Move indel starts to the position with the largest number of votes
		int newPredictedEventEnd = moveIndelStarts(indelAlns, currentPos, eventEnd, maxLength, maxI);
		if(maxI > 0) return answer;
		if(!fixedEvent && newPredictedEventEnd!=eventEnd) {
			//Update the predicted end the new indel event after realignment
			eventEnd = newPredictedEventEnd;
			answer = eventEnd-currentPos+1;
			if(currentPos==posPrint) System.out.println("ConciliateIndels. Predicted span for new indel: "+answer);	
		}
		/*if(eventChange) {
			if(currentPos==posPrint) System.out.println("Recalculating alignments with indels");
			maxLength = analyzeIndels(alignments,currentPos,eventEnd, varG, lengths,indelAlns,null);
			if(currentPos==posPrint)System.out.println("ConciliateIndels. Current pos: "+currentPos+" indelAlns: "+ indelAlns.size());
		}*/
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
						//if(aln.getFirst()==787293) System.out.println("Trying to move indel start for alignment of read "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" indel reference pos "+start+" current pileup pos: "+currentPos+" new indel start "+(currentPos+maxI)+" read pos: "+aln.getReadPosition(currentPos));
						boolean moved = aln.moveIndelStart(start,first+offset);
						if(first==posPrint && moved == false) System.err.println("WARN: Failed attempt to move indel start for alignment of read "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" indel reference pos "+start+" current pileup pos: "+first+" new indel start "+(first+offset)+" read pos: "+aln.getReadPosition(first));
						break;
					}
				}
				int alnRefLast = first;
				for(int start:indels.keySet()) {
					//This takes into account several close indel events within the same alignment
					if(start >= first && start <=alnRefLast+maxLength) {
						alnRefLast = indels.get(start).getLast();
					}
				}
				if(alnRefLast>answer) answer = alnRefLast;
			}
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
			String reference = seq.toString();
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
		int readFirst = aln.getReadPosition(currentPos);
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
		CharSequence sequenceBefore = genome.getReference(sequenceName, eventFirst-100, eventFirst);
		CharSequence sequenceAfter = genome.getReference(sequenceName, eventLast, eventLast+100);
		for(ReadAlignment aln:alignments) {
			int alnFirst = aln.getFirst();
			int alnLast = aln.getLast();
			String cigarStr = aln.getCigarString();
			
			if(eventFirst-alnFirst<bpForGoodRefAln) {
				boolean trimStart=true;
				int readPosAfter = aln.getReadPosition(eventLast);
				if(eventFirst==posPrint) System.out.println("IndelRealigner. realignEnds. Start falls within event. Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+". Event last: "+eventLast+" readPosAfter: "+readPosAfter);
				if(sequenceBefore!=null && readPosAfter>=bpForGoodRefAln) {
					CharSequence prefix = aln.getReadCharacters().subSequence(0, readPosAfter);
					int overlapLength = AbstractLimitedSequence.getOverlapLength (sequenceBefore,prefix);
					if(eventFirst == posPrint) System.out.println("IndelRealigner. realignEnds. overlap length: "+overlapLength);
					if(overlapLength>=bpForGoodRefAln) {
						int newAlnFirst = eventFirst-overlapLength+1; 
						aln.realignStart(newAlnFirst,overlapLength,readPosAfter);
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
				
			}
			if(alnLast-eventLast<bpForGoodRefAln) {
				boolean trimEnd = true;
				int length = aln.getReadLength();
				int readPosBefore = aln.getReadPosition(eventFirst);
				if(eventFirst==posPrint) System.out.println("IndelRealigner. realignEnds. End falls within event. Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+". Event last: "+eventLast+" readPosBefore: "+readPosBefore+" length: "+length);
				int lastReadBp = length - readPosBefore;
				if(sequenceAfter!=null && lastReadBp>=bpForGoodRefAln) {
					CharSequence suffix = aln.getReadCharacters().subSequence(readPosBefore+1,length);
					int overlapLength = AbstractLimitedSequence.getOverlapLength (suffix,sequenceAfter);
					if(eventFirst == posPrint) System.out.println("IndelRealigner. realignEnds. overlap length: "+overlapLength);
					if(overlapLength>=bpForGoodRefAln) {
						aln.realignEnd(readPosBefore, eventLast, overlapLength);
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
	}
}
