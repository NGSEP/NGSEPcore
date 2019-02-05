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
package ngsep.sequences;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import ngsep.alignments.ReadAlignment;

/**
 * Class able to build combined FM indexes for multiple small sequences
 * @author German Andrade
 * @author Jorge Duitama
 */
public class FMIndex implements Serializable
{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 6155320979838428304L;
	private QualifiedSequenceList sequencesWithNames;
	private List<Integer> sequenceLengths = new ArrayList<>();
	private List<FMIndexSingleSequence> internalIndexes = new ArrayList<>();
	private List<CombinedMultisequenceFMIndexMetadata> internalMetadata = new ArrayList<>();

	/**
	 * Loads the sequences in the given list to allow searches from these sequences
	 * @param sequences to add to the index. Each QualifiedSequence object in the list should have a name and its characters
	 */
	public void loadQualifiedSequenceList (QualifiedSequenceList sequences) {
		sequencesWithNames = sequences;
		StringBuffer internalSequence = new StringBuffer();
		CombinedMultisequenceFMIndexMetadata internalIdxMetadata = new CombinedMultisequenceFMIndexMetadata();
		int nI=0;
		int i=0;
		for(QualifiedSequence seq:sequences) {
			String next = seq.getCharacters().toString();
			if(internalSequence.length() + next.length() > 1000000000) {
				System.out.println("Building index for "+nI+" sequences. Total sequence length: "+internalSequence.length());
				long time = System.currentTimeMillis();
				FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
				System.out.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
				internalIndexes.add(index);
				internalMetadata.add(internalIdxMetadata);
				internalSequence = new StringBuffer();
				internalIdxMetadata = new CombinedMultisequenceFMIndexMetadata();
				nI=0;
			}
			internalSequence.append(next);
			internalIdxMetadata.addInputSequence(i, seq.getLength());
			sequenceLengths.add(seq.getLength());
			nI++;
			i++;
		}
		System.out.println("Building index for "+nI+" sequences. Total sequence length: "+internalSequence.length());
		long time = System.currentTimeMillis();
		FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
		System.out.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
		internalIndexes.add(index);
		internalMetadata.add(internalIdxMetadata);
	}
	public void loadUnnamedSequences (List<? extends CharSequence> sequences, int tally, int indexl) {
		int n = sequences.size();
		//TODO: Change to a LimitedSequence to make it space efficient
		StringBuffer internalSequence = new StringBuffer();
		CombinedMultisequenceFMIndexMetadata internalIdxMetadata = new CombinedMultisequenceFMIndexMetadata();
		int nI=0;
		for(int i=0;i<n;i++) {
			String next = sequences.get(i).toString();
			if(internalSequence.length() + next.length() > (100000000) ) {
				System.out.println("		Building index for "+nI+" sequences. Total sequence length: "+internalSequence.length());
				long time = System.currentTimeMillis();
				FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence,tally,indexl);
				System.out.println("		Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
				System.out.println("		-----------------");
				internalIndexes.add(index);
				internalMetadata.add(internalIdxMetadata);
				internalSequence = new StringBuffer();
				internalIdxMetadata = new CombinedMultisequenceFMIndexMetadata();
				nI=0;
			}
			internalSequence.append(next);
			internalIdxMetadata.addInputSequence(i, next.length());
			sequenceLengths.add(next.length());
			nI++;
		}
		System.out.println("		Building index for "+nI+" sequences. Total sequence length: "+internalSequence.length());
		long time = System.currentTimeMillis();
		FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence,tally,indexl);
		System.out.println("		Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
		internalIndexes.add(index);
		internalMetadata.add(internalIdxMetadata);
	}
	public List<ReadAlignment> search (String searchSequence) {
		return search(searchSequence, 0, sequenceLengths.size());
	}
	public List<ReadAlignment> search (String searchSequence, int firstIndex, int lastIndex) {
		List<ReadAlignment> alignments = new ArrayList<>();
		String searchUp = searchSequence.toUpperCase();
		for (int i=0;i<internalIndexes.size();i++) 
		{
			FMIndexSingleSequence idxSeq = internalIndexes.get(i);
			CombinedMultisequenceFMIndexMetadata metadata = internalMetadata.get(i);
			if(!metadata.overlapWithIndexes(firstIndex, lastIndex)) continue;
			Set<Integer> matches = idxSeq.search(searchUp);
			for (int internalPosMatch:matches) 
			{
				int [] realData = metadata.getSequenceIdxAndStart(internalPosMatch);
				if(realData==null) continue;
				if(realData[0]>=sequenceLengths.size()) throw new RuntimeException("Problem with internal index answer: "+realData[0]+"-"+realData[1]+". Absolute: "+internalPosMatch+" total length: "+metadata.getTotalLength()+" first idx: "+metadata.getFirstInputSequenceIdx()+" last idx: "+metadata.getLastInputSequenceIdx());
				int first = internalPosMatch-realData[1];
				int l = searchSequence.length();
				int last = first + l - 1;
				if(last>=sequenceLengths.get(realData[0])) continue;
				String seqName = ""+realData[0];
				if(sequencesWithNames!=null) seqName = sequencesWithNames.get(realData[0]).getName();
				ReadAlignment alignment = new ReadAlignment(seqName, first, last, l, 0);
				alignments.add(alignment);
			}
		}
		return alignments;
	}	
}
class CombinedMultisequenceFMIndexMetadata implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -4208212205942777492L;
	
	
	private int firstInputSequenceIdx=-1;
	private int lastInputSequenceIdx=-1;
	private List<Integer> sequenceStarts=new ArrayList<>();
	private int totalLength = 0;
	
	public void addInputSequence(int idx, int length) {
		if(firstInputSequenceIdx == -1) firstInputSequenceIdx = idx;
		lastInputSequenceIdx = idx;
		sequenceStarts.add(totalLength);
		totalLength += length;
	}
	
	public int [] getSequenceIdxAndStart (int absolutePosition) {
		if(absolutePosition>=totalLength) return null;
		int [] answer = new int [2];
		int i = Collections.binarySearch(sequenceStarts, absolutePosition);
		if(i>=0) {
			answer[0] = firstInputSequenceIdx+i;
			answer[1] = sequenceStarts.get(i);
		} else {
			//throw new RuntimeException("Position to search: "+absolutePosition+" i: "+i+ " start: "+sequenceStarts.get(-i-2));
			i = -i - 2;
			
			answer[0] = firstInputSequenceIdx+i;
			if(i<sequenceStarts.size()) answer[1] = sequenceStarts.get(i);
		}
		return answer;
	}
	public boolean overlapWithIndexes (int idxFirst, int idxLast) {
		return firstInputSequenceIdx<=idxLast && idxFirst<=lastInputSequenceIdx;
	}

	/**
	 * @return the firstInputSequenceIdx
	 */
	public int getFirstInputSequenceIdx() {
		return firstInputSequenceIdx;
	}

	/**
	 * @return the lastInputSequenceIdx
	 */
	public int getLastInputSequenceIdx() {
		return lastInputSequenceIdx;
	}

	/**
	 * @return the totalLength
	 */
	public int getTotalLength() {
		return totalLength;
	}
	
	
	
}
