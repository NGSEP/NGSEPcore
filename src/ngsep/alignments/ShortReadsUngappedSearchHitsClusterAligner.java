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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.genome.GenomicRegion;
import ngsep.genome.GenomicRegionImpl;
import ngsep.genome.io.SimpleGenomicRegionFileHandler;
import ngsep.sequences.LimitedSequence;

/**
 * @author German Andrade
 * @author Jorge Duitama
 *
 */
public class ShortReadsUngappedSearchHitsClusterAligner  implements UngappedSearchHitsClusterAligner {
	private int minMatchLength = 15;
	// Statistics
	private int fewMismatchesAlns = 0;
	private int completeAlns = 0;
	
	private Map<String, List<GenomicRegion>> knownSTRs;
	public Map<String, List<GenomicRegion>> getKnownSTRs() {
		return knownSTRs;
	}
	public void setKnownSTRs(Map<String, List<GenomicRegion>> knownSTRs) {
		this.knownSTRs = knownSTRs;
	}
	
	
	
	public int getFewMismatchesAlns() {
		return fewMismatchesAlns;
	}
	public int getCompleteAlns() {
		return completeAlns;
	}
	public ReadAlignment buildAlignment(CharSequence query, CharSequence subject,  UngappedSearchHitsCluster kmerHitsCluster) {
		int subjectIdx = kmerHitsCluster.getSubjectIdx();
		int first = kmerHitsCluster.getSubjectPredictedStart()+1;
		int last = kmerHitsCluster.getSubjectPredictedEnd();
		int lastPerfect = first+query.length()-1;
		List<Integer> alignment = new ArrayList<Integer>(1);
		alignment.add(ReadAlignment.getAlnValue(query.length(), ReadAlignment.ALIGNMENT_MATCH));
		ReadAlignment aln = buildAln(query, subjectIdx, subject, first, lastPerfect, alignment);
		//System.out.println("Built ungapped alignment  "+aln+" predicted limits: "+first+" "+last+" diff kmer hits: "+kmerHitsCluster.getNumDifferentKmers()+" allConsistent: "+kmerHitsCluster.isAllConsistent());
		if(aln!=null) {
			//TODO: use sequence name to find tandem repeats
			GenomicRegion region =findTandemRepeat("",first,last);
			if(region!=null) {
				ReadAlignment newaln=verifyShortTandemRepeats(subjectIdx,aln.getFirst(), aln.getLast(),query,subject,region);
				//System.out.println("Found overlapping tandem repeat at "+region.getSequenceName()+":"+region.getFirst()+"-"+region.getLast()+" new aln: "+newaln);
				if(newaln!=null) {
					return newaln;
				}
			}
			if(kmerHitsCluster.getCountKmerHitsCluster()>2 && kmerHitsCluster.isAllConsistent()) {
				int [] mismatches = countMismatches (query, subject, aln);
				if(mismatches !=null && mismatches[0]<0.05*query.length() && mismatches[1]+mismatches[2] < 0.1*query.length()) {
					int ends = mismatches[1]+mismatches[2]; 
					//if (ends > mismatches[0]) System.err.println("Problem counting mismatches for "+sequenceName+":"+first+" read: "+query+" mismatches: "+mismatches[0]+" "+mismatches[1]+" "+mismatches[2]);
					if (ends>0) aln = buildAln(query, subjectIdx, subject, first+mismatches[1], lastPerfect-mismatches[2], encodeAlignment(query.length(),mismatches));
					if(aln!=null) {
						aln.setAlignmentQuality((byte) Math.round(100-5*mismatches[0]));
						aln.setNumMismatches((short) mismatches[0]);
						fewMismatchesAlns++;
						//System.out.println("Few mismatches alignment: "+aln);
						return aln;
					}
				}
			}
		}
		//Perform smith waterman
		first = Math.max(1, first-3);
		last = Math.min(subject.length(), last+3);
		int d = last-first+1;
		//System.out.println("Building gapped alignment. Predicted limits: "+first+" "+last+" distance: "+d+" query length: "+query.length());
		if(d>1.5*query.length()) return null;
		if(d<0.5*query.length()) return null;
		CharSequence refSeq = subject.subSequence(first-1, last);
		
		//System.out.println("Aligning reference from "+first+" to "+last+ " to query. length: "+refSeq.length());
		completeAlns++;
		PairwiseAlignerAffineGap alignerFullRead = createAlignerFullRead(Math.max(query.length(), refSeq.length()));
		String [] rawAln = alignerFullRead.calculateAlignment(query, refSeq.toString());
		int mismatches = countMismatches(rawAln);
		if(mismatches>0.1*query.length()) return null;
		LinkedList<Integer> alnCodes = ReadAlignment.encodePairwiseAlignment(rawAln);
		aln = buildAln(query, subjectIdx, subject, first, last, alnCodes);
		if(aln==null) return null;
		if (!aln.clipBorders(minMatchLength)) return null;
		//System.out.println("New genomic coordinates : "+first+"-"+last+" CIGAR:" +aln.getCigarString()+" mismatches: "+mismatches);
		aln.setAlignmentQuality((byte) Math.round(100-5*mismatches));
		aln.setNumMismatches((short)mismatches);
		
		return aln;
	}
	public List<Integer> encodeAlignment(int length, int[] mismatches) {
		List<Integer> answer = new LinkedList<Integer>();
		int l2 = length-mismatches[1]-mismatches[2];
		if(mismatches[1]>0) {
			answer.add(ReadAlignment.getAlnValue(mismatches[1], ReadAlignment.ALIGNMENT_SKIPFROMREAD));
		}
		answer.add(ReadAlignment.getAlnValue(l2, ReadAlignment.ALIGNMENT_MATCH));
		if(mismatches[2]>0) {
			answer.add(ReadAlignment.getAlnValue(mismatches[2], ReadAlignment.ALIGNMENT_SKIPFROMREAD));
		}
		return answer;
	}
	private PairwiseAlignerAffineGap createAlignerFullRead(int capacity) {
		PairwiseAlignerAffineGap alignerFullRead = new PairwiseAlignerAffineGap(capacity);
		alignerFullRead.setForceStart2(false);
		alignerFullRead.setForceEnd2(false);
		return alignerFullRead;
	}
	public int countMismatches(String[] alignedSequences) {
		int answer = 0;
		boolean lastIsGap = true;
		for(int i=0;i<alignedSequences[0].length();i++) {
			char c1 = alignedSequences[0].charAt(i);
			char c2 = alignedSequences[1].charAt(i);
			if(c1==LimitedSequence.GAP_CHARACTER || c2 == LimitedSequence.GAP_CHARACTER) {
				if(!lastIsGap) answer+=2;
				lastIsGap = true;
			} else {
				if(c1!=c2) answer++;
				lastIsGap = false;
			}
		}
		if(lastIsGap) answer-=2;
		return answer;
	}
	public int [] countMismatches(CharSequence query, CharSequence subject, ReadAlignment aln) {
		int [] answer = {0,0,0};
		String refSeq = subject.subSequence(aln.getFirst()-1, aln.getLast()).toString();
		int lastMismatch = -1;
		boolean startAssigned = false;
		for (int i=0;i<query.length() && i<refSeq.length();i++ ) {
			if(query.charAt(i)!=refSeq.charAt(i)) {
				answer[0]++;
				lastMismatch=i;
			} else if (startAssigned==false && answer[0]+3<i) {
				answer[1]=lastMismatch+1;
				startAssigned = true;
			}
		}
		if (query.length()!=refSeq.length()) {
			answer[0]+=Math.abs(query.length()-refSeq.length());
			answer[2]=Math.max(0, query.length()-refSeq.length());
		} else {
			lastMismatch=refSeq.length();
			int numM =0;
			for (int i=query.length()-1;i>=0;i-- ) {
				if (query.charAt(i)!=refSeq.charAt(i)) {
					lastMismatch = i;
					numM++;
				} else {
					int revIdx = refSeq.length()-1-i;
					if (numM+3<revIdx) {
						answer[2]=refSeq.length()-lastMismatch;
						break;
					}
				}
			}
		}
		
		return answer;
	}
	
	public GenomicRegion findTandemRepeat(String sequenceName, int first, int last) {
		GenomicRegionImpl region = new GenomicRegionImpl(sequenceName, first, last);
		if(knownSTRs==null) return null;
		List<GenomicRegion> l =knownSTRs.get(sequenceName);
		if(l==null) return null;
		return binaryContains(l, 0, l.size()-1, region);
	}

	private GenomicRegion binaryContains(List<GenomicRegion> l, int left, int rigth, GenomicRegion element) {
		if(rigth>=left) {
			int middle = left +(rigth-left)/2;
			GenomicRegion actual = l.get(middle);
			if(isOverlappgingSorted(actual, element)||isOverlappgingSorted(element, actual	)) {
				return actual;
			}
			if(l.get(middle).getFirst()>element.getFirst()) {
				return binaryContains(l, left, middle-1, element);
			}
			return binaryContains(l, middle+1,rigth, element);
		}
		return null;
	}
	
	/**
	 * Creates an alignment taking into account that the region overlaps with a tandem repeat
	 * @param aln
	 * @param read
	 * @param qualityScores
	 * @param region
	 * @return
	 */
	public ReadAlignment verifyShortTandemRepeats(int subjectIdx, int first, int last, CharSequence read, CharSequence subject, GenomicRegion region) {
		int firstLeftPart = Math.max(first,1);
		int softClipLeft = 0;
		int softClipRight = 0;
		LinkedList<Integer> encodedLeftAln = null;
		LinkedList<Integer> encodedRightAln = null;
		int leftMismatches = 0;
		int rightMismatches = 0;
		
		if(first<region.getFirst()-5) {
			CharSequence refSeq = subject.subSequence(firstLeftPart-1, region.getFirst()-1).toString();
			if(refSeq!=null) {
				int endReadSegment= Math.min(read.length(), region.getFirst()-first+5);
				String readSegment = read.toString().substring(0,endReadSegment);
				//System.out.println(refSeq);
				//System.out.println(readSegment);
				PairwiseAlignerAffineGap alignerSTRsLeft = createAlignerLeftTR(Math.max(readSegment.length(), refSeq.length()));
				String [] alignmentLeft = alignerSTRsLeft.calculateAlignment(readSegment, refSeq.toString());
				leftMismatches = countMismatches(alignmentLeft);
				encodedLeftAln = ReadAlignment.encodePairwiseAlignment(alignmentLeft);
				int lastCode = encodedLeftAln.getLast();
				if (leftMismatches<=readSegment.length()/10 && ReadAlignment.getOperator(lastCode)==ReadAlignment.ALIGNMENT_INSERTION) {
					softClipLeft = ReadAlignment.getOperationLength(lastCode);
					encodedLeftAln.removeLast();
				} else {
					encodedLeftAln = null;
				}
				softClipLeft+=(read.length()-endReadSegment);
			}	
		}
		if(last>region.getLast()+5) {
			CharSequence refSeq = subject.subSequence(region.getLast(), last);
			if(refSeq!=null) {
				int startReadSegment= Math.max(0, read.length()-(last-region.getLast())-5);
				String readSegment = read.toString().substring(startReadSegment);
				//System.out.println(refSeq);
				//System.out.println(readSegment);
				PairwiseAlignerAffineGap alignerSTRsRight = createAlignerRightTR(Math.max(readSegment.length(), refSeq.length()));
				String [] alignmentRight = alignerSTRsRight.calculateAlignment(readSegment, refSeq.toString());
				rightMismatches = countMismatches(alignmentRight);
				encodedRightAln = ReadAlignment.encodePairwiseAlignment(alignmentRight);
				int firstCode = encodedRightAln.getFirst();
				if (rightMismatches<=readSegment.length()/10 && ReadAlignment.getOperator(firstCode)==ReadAlignment.ALIGNMENT_INSERTION) {
					softClipRight = ReadAlignment.getOperationLength(firstCode);
					encodedRightAln.removeFirst();
				} else {
					encodedRightAln=null;
				}
				softClipRight+=startReadSegment;
				
			}	
		}
		if(encodedLeftAln==null && encodedRightAln ==null) {
			return null;
		}
		if(encodedRightAln==null) {
			//Left alignment with right soft clip
			if(softClipLeft>0) encodedLeftAln.add(ReadAlignment.getAlnValue(softClipLeft, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
			ReadAlignment aln = buildAln(read, subjectIdx, subject, firstLeftPart, region.getFirst()-1, encodedLeftAln);
			if(aln==null) return null;
			//System.out.println("Left alignment new genomic coordinates : "+first+"-"+last+" CIGAR:" +aln.getCigarString());
			if (!aln.clipBorders(minMatchLength)) return null;
			aln.setAlignmentQuality((byte)(90-5*leftMismatches));
			aln.setNumMismatches((short) leftMismatches);
			return aln;
		}
		if(encodedLeftAln==null) {
			//Right alignment with left soft clip
			first = region.getLast()+1;
			if(softClipRight>0) encodedRightAln.addFirst(ReadAlignment.getAlnValue(softClipRight, ReadAlignment.ALIGNMENT_SKIPFROMREAD));
			ReadAlignment aln = buildAln(read, subjectIdx, subject, region.getLast()+1, last, encodedRightAln);
			if(aln==null) return null;
			//System.out.println("Right alignment new genomic coordinates : "+first+"-"+last+" CIGAR:" +aln.getCigarString());
			if (!aln.clipBorders(minMatchLength)) return null;
			aln.setAlignmentQuality((byte)(90-5*rightMismatches));
			aln.setNumMismatches((short) rightMismatches);
			return aln;
		}
		first = firstLeftPart;
		int alignedLeft = read.length()-softClipLeft;
		int alignedRight = read.length()-softClipRight;
		int middleLength = read.length()-alignedLeft-alignedRight;
		if(middleLength<0) return null;
		int difference = region.length()-middleLength;
		//System.out.println("Aligned left: "+alignedLeft+" aligned right: "+alignedRight+" middle length: "+middleLength+" difference: "+difference);
		LinkedList<Integer> alignmentList = new LinkedList<Integer>();
		alignmentList.addAll(encodedLeftAln);
		if(difference>0) {
			// Region length > middle length. Add deletion
			alignmentList.add(ReadAlignment.getAlnValue(difference, ReadAlignment.ALIGNMENT_DELETION));
			if(middleLength>0) alignmentList.add(ReadAlignment.getAlnValue(middleLength, ReadAlignment.ALIGNMENT_MATCH));
		} else if (difference<0) {
			// Region length < middle length. Add insertion
			alignmentList.add(ReadAlignment.getAlnValue(-difference, ReadAlignment.ALIGNMENT_INSERTION));
			if(region.length()>0) alignmentList.add(ReadAlignment.getAlnValue(region.length(), ReadAlignment.ALIGNMENT_MATCH));
		} else {
			if(middleLength>0) alignmentList.add(ReadAlignment.getAlnValue(middleLength, ReadAlignment.ALIGNMENT_MATCH));
		}
		
		alignmentList.addAll(encodedRightAln);
		short mismatches = (short) (leftMismatches+rightMismatches);
		ReadAlignment aln = buildAln(read, subjectIdx, subject, first, last, alignmentList);
		if(aln==null) return null;
		//System.out.println("Building alignment from first "+first+" last: "+last+" softClipLeft: "+softClipLeft+" softClip right "+softClipRight+" cigar "+aln.getCigarString());
		if (!aln.clipBorders(minMatchLength)) return null;
		aln.setAlignmentQuality((byte)(100-5*mismatches));
		aln.setNumMismatches(mismatches);
		return aln;
	}
	public PairwiseAlignerAffineGap createAlignerLeftTR (int capacity) {
		PairwiseAlignerAffineGap alignerSTRsLeft = new PairwiseAlignerAffineGap(capacity);
		alignerSTRsLeft.setForceEnd1(false);
		return alignerSTRsLeft;
	}
	
	public PairwiseAlignerAffineGap createAlignerRightTR (int capacity) {
		PairwiseAlignerAffineGap alignerSTRsRight = new PairwiseAlignerAffineGap(capacity);
		alignerSTRsRight.setForceStart1(false);
		return alignerSTRsRight;
	}
	
	
	public ReadAlignment buildAln(CharSequence query, int sequenceIdx, CharSequence subject, int first, int last, List<Integer> alignmentCodes) {
		if(first <=0) return null;
		ReadAlignment aln = new ReadAlignment(sequenceIdx, first, last, query.length(), 0);
		aln.setReadCharacters(query);
		if(alignmentCodes!=null)aln.setAlignment(alignmentCodes);
		if(aln.getLast()>subject.length()) return null;
		return aln;
	}
	public void loadSTRsFile(String strsFile) throws IOException {
		SimpleGenomicRegionFileHandler handler = new SimpleGenomicRegionFileHandler();
		knownSTRs=handler.loadRegionsAsMap(strsFile);
		flat();
	}

	private void flat() {
		Set<String>keys=knownSTRs.keySet();
		Iterator<String>it=keys.iterator();
		while(it.hasNext()) {
			String key = it.next();
			List<GenomicRegion> l =knownSTRs.get(key);
			List<GenomicRegion> newList = flat(l);
			while(isOverlappging(newList)) {
				newList=flat(newList);
			}
			knownSTRs.put(key, newList);
		}
	}
	private boolean isOverlappging(List<GenomicRegion>l) {
		for (int i = 0; i < l.size()-1; i++) {
			GenomicRegion current = l.get(i);
			GenomicRegion next = l.get(i+1);
			if(isOverlappgingSorted(current, next)) {
				return true;
			}
		}
		return false;
	}

	private boolean isOverlappgingSorted(GenomicRegion a,GenomicRegion b) {
		return b.getFirst()>=a.getFirst()&&b.getFirst()<=a.getLast();
	}

	private List<GenomicRegion> flat(List<GenomicRegion> l) {
		List<GenomicRegion> newList = new ArrayList<GenomicRegion>();
		for (int i = 0; i < l.size(); i++) {
			if(i+1==l.size())newList.add(l.get(i));
			else {
				GenomicRegion current =l.get(i);
				GenomicRegion next =l.get(i+1);
				if(isOverlappgingSorted(current, next)) {
					newList.add(new GenomicRegionImpl(current.getSequenceName(), current.getFirst(), Math.max(current.getLast(), next.getLast())));
					i++;
				}
				else {
					newList.add(current);
				}
			}	
		}

		return newList;
	}
	
}
