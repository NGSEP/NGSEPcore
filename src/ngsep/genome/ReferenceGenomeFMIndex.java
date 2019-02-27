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

import ngsep.alignments.ReadAlignment;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndexSingleSequence;
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
	private Map<String,FMIndexSingleSequence> internalIndexes = new HashMap<>();
	
	public ReferenceGenomeFMIndex (ReferenceGenome genome) {
		sequencesMetadata = genome.getSequencesMetadata();
		int n = genome.getNumSequences();
		
		for (int i = 0; i < n; i++) 
		{
			QualifiedSequence q = genome.getSequenceByIndex(i);
			CharSequence seqChars = q.getCharacters();
			FMIndexSingleSequence idxForward = new FMIndexSingleSequence(seqChars.toString().toUpperCase());
			internalIndexes.put(q.getName(),idxForward);
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
	 * @param searchSequence sequence to search
	 * @return List<ReadAlignment> Alignments of the given sequence to segments of sequences in this index
	 */
	public List<ReadAlignment> search (String searchSequence) {
		return search(searchSequence, false);
	}
	/**
	 * Searches the given sequence in the index
	 * @param searchSequence sequence to search
	 * @param searchReverseComplement If true, the reverse complement of the sequence will be calculated and searched in the index
	 * @return List<ReadAlignment> Alignments of the given sequence to segments of sequences in this index
	 */
	public List<ReadAlignment> search (String searchSequence, boolean searchReverseComplement)
	{
		List<ReadAlignment> alignments = new ArrayList<>();
		String searchUp = searchSequence.toUpperCase();
		int lq = searchUp.length();
		for (String seqName:internalIndexes.keySet()) 
		{
			FMIndexSingleSequence idxSeq = internalIndexes.get(seqName);
			Set<Integer> matches = idxSeq.search(searchUp);
			//System.out.println("Search: "+searchUp+" matches: "+matches);
			for (int internalPosMatch:matches) 
			{
				ReadAlignment alignment = new ReadAlignment(seqName, internalPosMatch+1, internalPosMatch+lq, lq, 0);
				alignment.setAlignmentQuality((short) 100);
				alignments.add(alignment);
			}
		}
		if(searchReverseComplement) {
			//Search the reverse complement
			searchUp = DNAMaskedSequence.getReverseComplement(searchUp);
			for (String seqName:internalIndexes.keySet()) 
			{
				FMIndexSingleSequence idxSeq = internalIndexes.get(seqName);
				Set<Integer> matches = idxSeq.search(searchUp);
				for (int internalPosMatch:matches) 
				{
					ReadAlignment alignment = new ReadAlignment(seqName, internalPosMatch+1, internalPosMatch+lq, lq, ReadAlignment.FLAG_READ_REVERSE_STRAND);
					alignment.setAlignmentQuality((short) 100);
					alignments.add(alignment);
				}
			}
		}
		return alignments;
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
	
	public static void main(String[] args) throws IOException
	{
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		ReferenceGenomeFMIndex f = new ReferenceGenomeFMIndex(genome);
		f.save(args[1]);
		f = loadFromBinaries(args[1]);
		List<ReadAlignment> a = f.search("ata");
		for (int i = 0; i < a.size(); i++) {
			System.out.println(a.get(i).getSequenceName()+" pos:"+a.get(i).getFirst()+" to: "+a.get(i).getLast()+" flags:"+a.get(i).getFlags());
		}
	}
	
}
