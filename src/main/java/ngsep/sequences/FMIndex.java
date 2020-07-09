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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPOutputStream;

import ngsep.main.io.ConcatGZIPInputStream;

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
	private int maxHitsQuery = 100000;
		
	public int getMaxHitsQuery() {
		return maxHitsQuery;
	}
	public void setMaxHitsQuery(int maxHitsQuery) {
		this.maxHitsQuery = maxHitsQuery;
		for(FMIndexSingleSequence internalIndex:internalIndexes) internalIndex.setMaxHitsQuery(maxHitsQuery);
	}
	/**
	 * Loads the sequences in the given list to allow searches from these sequences
	 * @param sequences to add to the index. Each QualifiedSequence object in the list should have a name and its characters
	 */
	public void loadQualifiedSequences (List<QualifiedSequence> sequences) {
		if(sequences instanceof QualifiedSequenceList) sequencesWithNames = (QualifiedSequenceList)sequences;
		else {
			sequencesWithNames = new QualifiedSequenceList();
			sequencesWithNames.addAll(sequences);
		}
		StringBuffer internalSequence = new StringBuffer();
		CombinedMultisequenceFMIndexMetadata internalIdxMetadata = new CombinedMultisequenceFMIndexMetadata();
		int nI=0;
		int i=0;
		for(QualifiedSequence seq:sequences) {
			String next = seq.getCharacters().toString();
			if(nI>0 && internalSequence.length() + next.length() > 100000000) {
				System.err.println("Building index for "+nI+" sequences. Total sequence length: "+internalSequence.length());
				long time = System.currentTimeMillis();
				FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
				index.setMaxHitsQuery(maxHitsQuery);
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
			index.setMaxHitsQuery(maxHitsQuery);
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
	public List<UngappedSearchHit> exactSearch (String query) {
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
	public List<UngappedSearchHit> exactSearch (String query, int firstIndex, int lastIndex) {
		List<UngappedSearchHit> hits = new ArrayList<>();
		for (int i=0;i<internalIndexes.size() && hits.size()<maxHitsQuery;i++) 
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
				
				//ReadAlignment alignment = new ReadAlignment(seqName, first, last, searchLength, 0);
				UngappedSearchHit hit = new UngappedSearchHit(query, sequenceIdx, start);
				if(sequencesWithNames!=null) hit.setSequenceName(sequencesWithNames.get(sequenceIdx).getName());
				hits.add(hit);
				if(hits.size()>=maxHitsQuery) break;
			}
		}
		for(UngappedSearchHit hit: hits) hit.setTotalHitsQuery(hits.size());
		return hits;
	}
	/**
	 * Return the sequence with the given name
	 * @param sequenceName Name of the sequence to search
	 * @return CharSequence sequence with the given name. Null if the name is not found
	 */
	public CharSequence getSequence(String sequenceName) {
		QualifiedSequence seq = sequencesWithNames.get(sequenceName);
		if(seq==null) return null;
		return seq.getCharacters();
	}
	/**
	 * Return the subsequence of the indexed sequence between the given genomic coordinates
	 * @param sequenceName Name of the sequence to search
	 * @param first position of the sequence (1-based, included)
	 * @param last position of the sequence (1-based, included)
	 * @return CharSequence segment of the given sequence between the given coordinates
	 */
	public CharSequence getSequence(String sequenceName, int first, int last) {
		QualifiedSequence seq = sequencesWithNames.get(sequenceName);
		if(seq==null) return null;
		CharSequence characters = seq.getCharacters();
		if(first <=0 || first>characters.length()) return null;
		if(last <=0 || last>characters.length()) return null;
		return characters.subSequence(first-1, last);
		
	}	
	public void save (String filename) throws IOException {
		try(OutputStream os = new GZIPOutputStream(new FileOutputStream(filename));
			PrintStream out = new PrintStream(os)) {
			save(out);
		}
	}
	public void save (PrintStream out) {
		out.println("#COMPOUNDINDEX\t"+maxHitsQuery);
		for (CombinedMultisequenceFMIndexMetadata metadata:internalMetadata) {
			metadata.save(out);
		}
		out.println("#INTERNALINDEXES");
		int i=0;
		for(FMIndexSingleSequence index:internalIndexes) {
			System.out.println("Saving internal index: "+i);
			index.save(out);
			System.out.println("Saved internal index: "+i);
			i++;
		}
	}
	public static FMIndex load (QualifiedSequenceList sequences, String indexFile) throws IOException {
		FMIndex index = new FMIndex();
		index.sequencesWithNames = sequences;
		for(QualifiedSequence seq:sequences) index.sequenceLengths.add(seq.getLength());
		try (FileInputStream fis = new FileInputStream(indexFile);
			 ConcatGZIPInputStream gzis = new ConcatGZIPInputStream(fis);
			 InputStreamReader isr = new InputStreamReader(gzis);
			 BufferedReader reader = new BufferedReader(isr)) {
			String line = reader.readLine();
			if(line==null) throw new IOException("Empty index file");
			if(!line.startsWith("#COMPOUNDINDEX")) throw new IOException("#COMPOUNDINDEX section not found. Line: "+line);
			String [] items = line.split("\t");
			index.maxHitsQuery = Integer.parseInt(items[1]);
			line = reader.readLine();
			while (line!=null && !line.equals("#INTERNALINDEXES")) {
				items = line.split("\t");
				if(!"#METADATA".equals(items[0])) throw new IOException("Unexpected line reading metadata. Line: "+line);
				CombinedMultisequenceFMIndexMetadata metadata = new CombinedMultisequenceFMIndexMetadata();
				for(int i=1;i<items.length;i+=2) {
					metadata.addInputSequence(Integer.parseInt(items[i]), Integer.parseInt(items[i+1]));
				}
				index.internalMetadata.add(metadata);
				line = reader.readLine();
			}
			if(line == null) throw new IOException("Unexpected end of file reading metadata.");
			while(true) {
				FMIndexSingleSequence internalIndex = FMIndexSingleSequence.load(reader);
				if(internalIndex==null) break;
				System.out.println("Loaded internal index: "+index.internalIndexes.size());
				index.internalIndexes.add(internalIndex);			
			}
			if(index.internalMetadata.size()!=index.internalIndexes.size())  throw new IOException("Inconsistent metadata and internal indexes. Metadata entries: "+index.internalMetadata.size()+" indexes: "+index.internalIndexes.size());
		}
		return index;
	}
}
class CombinedMultisequenceFMIndexMetadata implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -4208212205942777492L;
	
	
	private int firstInputSequenceIdx=-1;
	private int lastInputSequenceIdx=-1;
	private List<Integer> idxs=new ArrayList<>();
	private List<Integer> lengths=new ArrayList<>();
	private List<Integer> sequenceStarts=new ArrayList<>();
	private int totalLength = 0;
	
	public void addInputSequence(int idx, int length) {
		if(firstInputSequenceIdx == -1) firstInputSequenceIdx = idx;
		lastInputSequenceIdx = idx;
		idxs.add(idx);
		lengths.add(length);
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
	
	public void save (PrintStream out) {
		out.print("#METADATA");
		for(int i=0;i<idxs.size();i++) out.print("\t"+idxs.get(i)+"\t"+lengths.get(i));
		out.println();
	}
	
}
