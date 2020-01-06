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
package ngsep.genome;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.sequences.FMIndexSingleSequence;
import ngsep.sequences.FMIndexUngappedSearchHit;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

/**
 * FMIndex for reference genomes
 * @author German Andrade
 * @author Jorge Duitama
 */
public class ReferenceGenomeFMIndex implements Serializable {
	/**
	 * Serial number
	 */
	private static final long serialVersionUID = 5577026857894649939L;
	private QualifiedSequenceList sequencesMetadata;
	private Map<Integer,FMIndexSingleSequence> internalIndexes = new HashMap<>();
	
	public ReferenceGenomeFMIndex (ReferenceGenome genome) {
		sequencesMetadata = genome.getSequencesMetadata();
		int n = genome.getNumSequences();
		
		for (int i = 0; i < n; i++) 
		{
			QualifiedSequence q = genome.getSequenceByIndex(i);
			CharSequence seqChars = q.getCharacters();
			FMIndexSingleSequence idx = new FMIndexSingleSequence(seqChars.toString().toUpperCase());
			internalIndexes.put(i,idx);
		}
	}
	
	/**
	 * Loads an instance of the FMIndex from a serialized binary file
	 * @param filename Binary file with the serialization of an FMIndex
	 * @return FMIndex serialized in the given file
	 * @throws IOException If there were errors reading the file
	 */
	public static ReferenceGenomeFMIndex loadFromBinaries(String filename) throws IOException
	{
		ReferenceGenomeFMIndex fmIndex;
		try (FileInputStream fis = new FileInputStream(filename);
			 ObjectInputStream ois = new ObjectInputStream(fis);) {
			fmIndex = (ReferenceGenomeFMIndex) ois.readObject();
		} catch (ClassNotFoundException e) {
			throw new RuntimeException("FMIndex class not found",e);
		}
		return fmIndex;
	}
	
	/**
	 * Saves this FM-Index as a serializable object
	 * @param filename
	 * @throws IOException
	 */
	public void save (String filename) throws IOException 
	{
		try (FileOutputStream fos = new FileOutputStream( filename );
			 ObjectOutputStream oos = new ObjectOutputStream( fos );) {
			oos.writeObject( this );
		}
	}
	/**
	 * @return The list of sequences and lengths related to the reference genome
	 */
	public QualifiedSequenceList getSequencesMetadata() {
		return sequencesMetadata;
	}
	/**
	 * Returns the length of the reference sequence with the given name
	 * @param sequenceName Name of the sequence
	 * @return int length of the sequence. Zero if there are no reference sequences with the given name 
	 */
	public int getReferenceLength(String sequenceName) {
		QualifiedSequence seq = sequencesMetadata.get(sequenceName);
		if(seq==null) return 0;
		return seq.getLength();
	}
	/**
	 * Searches the given sequence in the index
	 * This search is case sensitive
	 * @param searchSequence sequence to search
	 * @return List<ReadAlignment> Alignments of the given sequence to segments of sequences in this index
	 */
	public List<FMIndexUngappedSearchHit> exactSearch (String searchSequence) {
		List<FMIndexUngappedSearchHit> hits = new ArrayList<>();
		for (int i=0;i<sequencesMetadata.size();i++) {
			QualifiedSequence seq = sequencesMetadata.get(i);
			FMIndexSingleSequence idxSeq = internalIndexes.get(i);
			Set<Integer> matches = idxSeq.exactSearch(searchSequence);
			for (int internalPosMatch:matches) {
				FMIndexUngappedSearchHit hit = new FMIndexUngappedSearchHit(searchSequence, i, seq.getName(), internalPosMatch);
				hits.add(hit);
			}
		}
		return hits;
	}
	
	/**
	 * Return the subsequence of the indexed sequence between the given genomic coordinates
	 * @param sequenceName Name of the sequence to search
	 * @param first position of the sequence (1-based, included)
	 * @param last position of the sequence (1-based, included)
	 * @return CharSequence segment of the given sequence between the given coordinates
	 */
	public CharSequence getSequence (String sequenceName, int first, int last) {
		FMIndexSingleSequence idxSeq = internalIndexes.get(sequenceName);
		if(idxSeq==null) return null;
		CharSequence seq = idxSeq.getSequence(first-1, last);
		return seq;
	}

	public boolean isValidAlignment(String sequenceName,int last) {
		FMIndexSingleSequence internalIndex = internalIndexes.get(sequenceName);
		if(internalIndex == null) return false;
		return internalIndex.getSequenceLength()>=last;
	}
	
}
