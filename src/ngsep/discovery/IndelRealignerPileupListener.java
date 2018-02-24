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
			int conciliatedSpan = conciliateIndels(pileup,predictedEventEnd,var!=null);
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
				if (currentPos<first) break;
				else if(first<=currentPos && currentPos<=last) return var;
				//System.out.println("Passed STR at "+currentReferenceName+":"+first+"-"+last+" current reference pos: "+currentReferencePos);
				idxNextVariant++;
			}
		}
		return null;
	}
	/**
	 * 
	 * @param predictedEventEnd
	 * @return int The reference span if indels are called (zero otherwise).
	 */
	private int conciliateIndels(PileupRecord pileup, int eventEnd, boolean fixedEvent) {
		int answer = 0;
		boolean eventChange = false;
		int currentPos = pileup.getPosition();
		
		
		Set<Integer> lengths = new TreeSet<Integer>();
		List<ReadAlignment> alignments = pileup.getAlignments();
		List<ReadAlignment> indelAlns = new ArrayList<ReadAlignment>();
		List<ReadAlignment> nonIndelGoodSpanAlns = new ArrayList<ReadAlignment>();
		List<ReadAlignment> nonIndelBadSpanAlns = new ArrayList<ReadAlignment>();
		//Vote for possible indel starts
		int [] votes = new int [eventEnd-currentPos+1];
		
		int maxLength = analyzeIndels(alignments,currentPos,eventEnd, lengths,indelAlns,nonIndelGoodSpanAlns,nonIndelBadSpanAlns,votes);
		
		if(currentPos==posPrint)System.out.println("ConciliateIndels. Current pos: "+currentPos+" Max length: "+maxLength+" lengths: "+ lengths.size());
		
		if(lengths.size()==0) return answer;
		
		
		int maxI = 0;
		if(!fixedEvent) {
			maxI = AlignmentsPileupGenerator.getMaxIndex(votes,-1);
			if(lengths.size()>1) {
				//Look for new tandem repeats
				if(currentPos==posPrint) System.out.println("ConciliateIndels. Checking tandem repeats.");
				int newSpan = lookForNewSTR(pileup, indelAlns, maxLength);
				if(newSpan>1) {
					//eventEnd = currentPos+newSpan-1;
					maxI = 0;
					answer = newSpan;
					eventChange = true;
					fixedEvent = true;
					pileup.setSTR(true);
					pileup.setNewSTR(true);
					if(currentPos==posPrint) System.out.println("New tandem repeat identified from reads at "+pileup.getSequenceName()+":"+currentPos+" New max span: "+answer);
				}
			}
		}
		if(currentPos==posPrint)System.out.println("ConciliateIndels. Current pos: "+currentPos+" Start max votes: "+(currentPos+maxI)+" indelAlns: "+ indelAlns.size()+" non indel good alignments "+nonIndelGoodSpanAlns.size()+" nonIndel partialAlns: "+nonIndelBadSpanAlns.size());
		
		
		
		//Move indel starts to the position with the largest number of votes. Use initial event end to be consistent with votes
		for(ReadAlignment aln:indelAlns) {
			Map<Integer,GenomicVariant> indels = aln.getIndelCalls();
			if(indels!=null) {
				for(int start:indels.keySet()) {
					if(start >= currentPos && start <=eventEnd) {		
						//if(aln.getFirst()==265419) System.out.println("Trying to move indel start for alignment of read "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" indel reference pos "+start+" current pileup pos: "+currentPos+" new indel start "+(currentPos+maxI)+" read pos: "+aln.getReadPosition(currentPos));
						boolean moved = aln.moveIndelStart(start,currentPos+maxI);
						if(currentPos==posPrint && moved == false) System.err.println("WARN: Failed attempt to move indel start for alignment of read "+aln.getReadName()+" at "+aln.getSequenceName()+":"+aln.getFirst()+" indel reference pos "+start+" current pileup pos: "+currentPos+" new indel start "+(currentPos+maxI)+" read pos: "+aln.getReadPosition(currentPos));
						break;
					}
				}
			}
		}
		if(maxI > 0) return answer;
		if(eventChange) eventEnd = currentPos+answer-1;
		if(currentPos==posPrint) System.out.println("ConciliateIndels. Calculating span. Number of lengths: "+lengths.size()+ ". Max indel length: "+maxLength);
		
		if(!fixedEvent) {
			//Update the predicted end the new indel event
			int newPredictedEventEnd = calculateEndOfIndel(currentPos, indelAlns, maxLength);
			if(newPredictedEventEnd!=eventEnd) {
				eventChange = true;
				eventEnd = newPredictedEventEnd;
				answer = eventEnd-currentPos+1;
				if(currentPos==posPrint) System.out.println("ConciliateIndels. Predicted span for new indel: "+answer);
			}	
		}
		if(eventChange) {
			if(currentPos==posPrint) System.out.println("Recalculating alignments with indels");
			maxLength = analyzeIndels(alignments,currentPos,eventEnd, lengths,indelAlns,nonIndelGoodSpanAlns,nonIndelBadSpanAlns,null);
			if(currentPos==posPrint)System.out.println("ConciliateIndels. Current pos: "+currentPos+" indelAlns: "+ indelAlns.size()+" non indel good alignments "+nonIndelGoodSpanAlns.size()+" nonIndel partialAlns: "+nonIndelBadSpanAlns.size());
		}
		
		//Soft clip not good reference alignments
		trimAlignmentEndsWithinEvent(currentPos, eventEnd, nonIndelBadSpanAlns);
		return answer;
	}
	
	private int analyzeIndels(List<ReadAlignment> alignments, int eventStart, int eventEnd, Set<Integer> lengths, List<ReadAlignment> indelAlns,
			List<ReadAlignment> nonIndelGoodSpanAlns, List<ReadAlignment> nonIndelBadSpanAlns, int[] votes) {
		
		int maxLength = 0;
		lengths.clear();
		indelAlns.clear();
		nonIndelGoodSpanAlns.clear();
		nonIndelBadSpanAlns.clear();
		if(votes!=null)Arrays.fill(votes,0);
		for(ReadAlignment aln:alignments) {
			boolean indelFound = false;
			Map<Integer,GenomicVariant> indels = aln.getIndelCalls();
			if(indels!=null) {
				for(int start:indels.keySet()) {
					if(start >= eventStart && start <=eventEnd) {
						indelFound = true;
						int length = indels.get(start).length();
						lengths.add(length);
						if(maxLength<length) maxLength = length;
						int i = start-eventStart;
						if(eventStart==posPrint) System.out.println("Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+" CIGAR: "+aln.getCigarString()+" Insertion start: "+start+" vote: "+i);
						//Test: Only primary alignments can vote for indel start
						//if(!aln.getSAMRecord().getNotPrimaryAlignmentFlag()) 
						if(votes!=null)votes[i]++;
						break;
					}
				}
			}
			if(indelFound) indelAlns.add(aln);
			else {
				if(eventStart - aln.getFirst()>=bpForGoodRefAln && aln.getLast()-eventEnd>=bpForGoodRefAln ) {
					if(eventStart==posPrint) System.out.println("Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+" reference spans. Event end: "+eventEnd);
					nonIndelGoodSpanAlns.add(aln);
				} else {
					if(eventStart==posPrint) System.out.println("Read name: "+aln.getReadName()+". Aln limits: "+aln.getFirst()+"-"+aln.getLast()+" partial reference. Event end: "+eventEnd);
					nonIndelBadSpanAlns.add(aln);
				}
			}
		}
		return maxLength;
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
	
	private int calculateEndOfIndel(int currentPos, List<ReadAlignment> alignments, int maxLength) {
		int predictedIndelEnd = currentPos+1;
		for(ReadAlignment aln:alignments) {
			int alnRefLast = currentPos;
			Map<Integer,GenomicVariant> indels = aln.getIndelCalls();
			if(indels!=null) {
				for(int start:indels.keySet()) {
					//This takes into account several close indel events within the same alignment
					if(start >= currentPos && start <=alnRefLast+maxLength) {
						alnRefLast = indels.get(start).getLast();
					}
				}
			}
			if(alnRefLast>predictedIndelEnd) predictedIndelEnd = alnRefLast;
		}
		return predictedIndelEnd;
	}
	
	private void trimAlignmentEndsWithinEvent(int eventFirst, int eventLast, List<ReadAlignment> alns) {
		if(eventFirst == posPrint) System.out.println("IndelRealigner. Trimming ends for "+alns.size()+" alignments");
		for(ReadAlignment aln:alns) {
			int alnFirst = aln.getFirst();
			int alnLast = aln.getLast();
			
			if(eventFirst-alnFirst<bpForGoodRefAln) {
				int ignoreBP = eventLast-alnFirst+1;
				ignoreBP+=aln.getSoftClipStart();
				byte bpToIgnoreStart = (byte)Math.max(aln.getBasesToIgnoreStart(), ignoreBP);
				
				//if(bpToIgnoreStart>10)System.err.println("WARN: Ignoring "+bpToIgnoreStart+" base pairs at the start of alignment of read "+aln.getSAMRecord().getReadName()+" at "+aln.getReferenceName()+":"+aln.getAlignmentStart()+ " Current CIGAR: "+aln.getSAMRecord().getCigarString()+" indel alns: "+numIndelAlns+" non indel alns: "+alns.size()+" event first: "+eventFirst+" event last: "+eventLast);
				aln.setBasesToIgnoreStart(bpToIgnoreStart);
				if(eventFirst == posPrint) System.out.println("IndelRealigner. Trimmed "+bpToIgnoreStart+" at the start of alignment with coordinates: "+alnFirst+"-"+alnLast);
			}
			if(alnLast-eventLast<bpForGoodRefAln) {
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
