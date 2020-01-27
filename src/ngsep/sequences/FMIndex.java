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
			if(nI>0 && internalSequence.length() + next.length() > 1000000000) {
				System.err.println("Building index for "+nI+" sequences. Total sequence length: "+internalSequence.length());
				long time = System.currentTimeMillis();
				FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
				System.err.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
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
		if(nI>0) {
			System.err.println("Building index for "+nI+" sequences. Total sequence length: "+internalSequence.length());
			long time = System.currentTimeMillis();
			FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
			System.err.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
			internalIndexes.add(index);
			internalMetadata.add(internalIdxMetadata);
		}
	}
	public void loadUnnamedSequences (List<? extends CharSequence> sequences, int tally, int suffixFraction) {
		int n = sequences.size();
		StringBuffer internalSequence = new StringBuffer();
		CombinedMultisequenceFMIndexMetadata internalIdxMetadata = new CombinedMultisequenceFMIndexMetadata();
		int nI=0;
		for(int i=0;i<n;i++) {
			String next = sequences.get(i).toString();
			if(nI>0 && internalSequence.length() + next.length() > (100000000) ) {
				System.err.println("Building index for "+nI+" sequences. Total sequence length: "+internalSequence.length());
				long time = System.currentTimeMillis();
				FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence,tally,suffixFraction);
				System.err.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
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
		if(nI>0) {
			System.err.println("Building index for "+nI+" sequences. Total sequence length: "+internalSequence.length());
			long time = System.currentTimeMillis();
			FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence,tally,suffixFraction);
			System.err.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
			internalIndexes.add(index);
			internalMetadata.add(internalIdxMetadata);
		}
	}
	/**
	 * Searches the given sequence against this FMindex.
	 * This search is case sensitive.
	 * @param query Sequence to search
	 * @return List<FMIndexUngappedSearchHit> exact hits to the sequences indexed by this FMIndex 
	 */
	public List<FMIndexUngappedSearchHit> exactSearch (String query) {
		return exactSearch(query, 0, sequenceLengths.size());
	}
	/**
	 * Searches the given sequence against the given subsequences of this FMindex.
	 * This search is case sensitive.
	 * @param query Sequence to search
	 * @param firstIndex of the subject sequence to look for
	 * @param lastIndex of the subject sequence to look for
	 * @return
	 */
	public List<FMIndexUngappedSearchHit> exactSearch (String query, int firstIndex, int lastIndex) {
		List<FMIndexUngappedSearchHit> hits = new ArrayList<>();
		for (int i=0;i<internalIndexes.size();i++) 
		{
			FMIndexSingleSequence idxSeq = internalIndexes.get(i);
			CombinedMultisequenceFMIndexMetadata metadata = internalMetadata.get(i);
			if(!metadata.overlapWithIndexes(firstIndex, lastIndex)) continue;
			Set<Integer> matches = idxSeq.exactSearch(query);
			for (int internalPosMatch:matches) 
			{
				int [] realData = metadata.getSequenceIdxAndStart(internalPosMatch);
				if(realData==null) continue;
				int sequenceIdx = realData[0];
				int sequenceStart = realData[1];
				if(sequenceIdx>=sequenceLengths.size()) throw new RuntimeException("Problem with internal index answer: "+realData[0]+"-"+realData[1]+". Absolute: "+internalPosMatch+" total length: "+metadata.getTotalLength()+" first idx: "+metadata.getFirstInputSequenceIdx()+" last idx: "+metadata.getLastInputSequenceIdx());
				//Match to other sequences sharing internal index with queried sequence
				if(sequenceIdx<firstIndex) continue;
				if(sequenceIdx>lastIndex) continue;
				int start = internalPosMatch-sequenceStart;
				int sequenceLength = sequenceLengths.get(sequenceIdx); 
				int queryLength = query.length();
				int last = start + queryLength - 1;
				//Match with artificial concatenation between sequences
				if(last>=sequenceLength) continue;
				String seqName = ""+sequenceIdx;
				if(sequencesWithNames!=null) seqName = sequencesWithNames.get(sequenceIdx).getName();
				//ReadAlignment alignment = new ReadAlignment(seqName, first, last, searchLength, 0);
				FMIndexUngappedSearchHit hit = new FMIndexUngappedSearchHit(query, sequenceIdx, seqName, start);
				hits.add(hit);
			}
		}
		return hits;
	}
	public CharSequence getSequence(String sequenceName, int first, int last) {
		QualifiedSequence seq = sequencesWithNames.get(sequenceName);
		if(seq==null) return null;
		CharSequence characters = seq.getCharacters();
		if(first <=0 || first>characters.length()) return null;
		if(last <=0 || last>characters.length()) return null;
		return characters.subSequence(first-1, last);
		
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
