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
import ngsep.genome.ReferenceGenome;

/**
 * Class able to build FM indexes for multiple sequences
 * @author German Andrade
 * @author Jorge Duitama
 */
public class FMIndex implements Serializable 
{
	/**
	 * Serial number
	 */
	private static final long serialVersionUID = 5577026857894649939L;
	private List<FMIndexSingleSequence> internalIndexes = new ArrayList<>();
	private List<List<SequenceMetadata>> realSequencesMap = new ArrayList<>();
	private Map<String, SequenceMetadata> realSequencesByKey = new HashMap<>();

	/**
	 * Loads the given genome in this FMIndex
	 * @param filename with the fasta file with genome to load
	 * @throws IOException If the genome file can not be read
	 */
	public void loadGenome(String genomeFilename) throws IOException 
	{
		ReferenceGenome referenceGenome = new ReferenceGenome(genomeFilename);
		int n = referenceGenome.getNumSequences();
		
		for (int i = 0; i < n; i++) 
		{
			QualifiedSequence q = referenceGenome.getSequenceByIndex(i);
			DNAMaskedSequence seqChars = (DNAMaskedSequence)q.getCharacters();
			FMIndexSingleSequence seqForward = new FMIndexSingleSequence(seqChars.toString().toUpperCase());
			internalIndexes.add(seqForward);
			
			SequenceMetadata seqM = new SequenceMetadata(q.getName(), q.getLength(), false, realSequencesMap.size(),0);
			realSequencesByKey.put(seqM.getUniqueKey(), seqM);
			List<SequenceMetadata> metadata = new ArrayList<>();
			metadata.add(seqM);
			realSequencesMap.add(metadata);
			seqChars = seqChars.getReverseComplement();
			FMIndexSingleSequence seqReverse = new FMIndexSingleSequence(seqChars.toString().toUpperCase());
			internalIndexes.add(seqReverse);
			seqM = new SequenceMetadata(q.getName(), q.getLength(), true, realSequencesMap.size(),0);
			realSequencesByKey.put(seqM.getUniqueKey(), seqM);
			metadata = new ArrayList<>();
			metadata.add(seqM);
			realSequencesMap.add(metadata);
		}
	}
	/**
	 * Loads the sequences in the given list to allow searches from these sequences
	 * @param sequences to add to the index. Each QualifiedSequence object in the list should have a name and its characters
	 */
	public void loadQualifiedSequenceList (QualifiedSequenceList sequences) {
		StringBuffer internalSequence = new StringBuffer();
		List<SequenceMetadata> internalSequenceMetadata = new ArrayList<>();
		for(QualifiedSequence seq:sequences) {
			String next = seq.getCharacters().toString();
			if(internalSequence.length() + next.length() > 1000000000) {
				System.out.println("Building index for "+internalSequenceMetadata.size()+" sequences. Total sequence length: "+internalSequence.length());
				long time = System.currentTimeMillis();
				FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
				System.out.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
				internalIndexes.add(index);
				realSequencesMap.add(internalSequenceMetadata);
				internalSequence = new StringBuffer();
				internalSequenceMetadata = new ArrayList<>();
			}
			SequenceMetadata seqM = new SequenceMetadata(seq.getName(), next.length(), false,realSequencesMap.size(),internalSequence.length());
			realSequencesByKey.put(seqM.getUniqueKey(), seqM);
			internalSequence.append(next);
			internalSequenceMetadata.add(seqM);
			
		}
		System.out.println("Building index for "+internalSequenceMetadata.size()+" sequences. Total sequence length: "+internalSequence.length());
		long time = System.currentTimeMillis();
		FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
		System.out.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
		internalIndexes.add(index);
		realSequencesMap.add(internalSequenceMetadata);
	}
	public void loadUnnamedSequences (String groupName, List<? extends CharSequence> sequences) {
		int n = sequences.size();
		StringBuffer internalSequence = new StringBuffer();
		List<SequenceMetadata> internalSequenceMetadata = new ArrayList<>();
		for(int i=0;i<n;i++) {
			String next = sequences.get(i).toString();
			if(internalSequence.length() + next.length() > (100000000) ) {
				System.out.println("Building index for "+internalSequenceMetadata.size()+" sequences. Total sequence length: "+internalSequence.length());
				long time = System.currentTimeMillis();
				FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
				System.out.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
				internalIndexes.add(index);
				realSequencesMap.add(internalSequenceMetadata);
				internalSequence = new StringBuffer();
				internalSequenceMetadata = new ArrayList<>();
			}
			SequenceMetadata seqM = new SequenceMetadata(groupName+"_"+i, next.length(), false,realSequencesMap.size(),internalSequence.length());
			realSequencesByKey.put(seqM.getUniqueKey(), seqM);
			internalSequence.append(next);
			internalSequenceMetadata.add(seqM);
			
		}
		System.out.println("Building index for "+internalSequenceMetadata.size()+" sequences. Total sequence length: "+internalSequence.length());
		long time = System.currentTimeMillis();
		FMIndexSingleSequence index = new FMIndexSingleSequence(internalSequence);
		System.out.println("Built index in "+(System.currentTimeMillis()-time)+" milliseconds");
		internalIndexes.add(index);
		realSequencesMap.add(internalSequenceMetadata);
	}
	
	/**
	 * Loads an instance of the FMIndex from a serialized binary file
	 * @param filename Binary file with the serialization of an FMIndex
	 * @return FMIndex serialized in the given file
	 * @throws IOException If there were errors reading the file
	 */
	public static FMIndex loadFromBinaries(String filename) throws IOException
	{
		FMIndex fmIndex = new FMIndex();
		try (FileInputStream fis = new FileInputStream(filename);
			 ObjectInputStream ois = new ObjectInputStream(fis);) {
			fmIndex = (FMIndex) ois.readObject();
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
	
	public List<ReadAlignment> search (String searchSequence)
	{
		
		List<ReadAlignment> alignments = new ArrayList<>();
		String searchUp = searchSequence.toUpperCase();
		for (int i=0;i<internalIndexes.size();i++) 
		{
			FMIndexSingleSequence idxSeq = internalIndexes.get(i);
			Set<Integer> matches = idxSeq.search(searchUp);
			for (int internalPosMatch:matches) 
			{
				ReadAlignment alignment = buildAlignmentFromMetadata(searchUp, i,internalPosMatch);
				if(alignment!=null) alignments.add(alignment);
			}
		}
		return alignments;
	}
	
	/**
	 * Return the subsequence of the indexed sequence between the given coordinates
	 * @param sequenceName Name of the sequence to search
	 * @param start position of the sequence (0-based, included)
	 * @param end position of the sequence (0-based, excluded)
	 * @param negativeStrand if true, the negative strand is requested
	 * @return CharSequence segment of the given sequence between the given coordinates
	 */
	public CharSequence getSequence (String sequenceName, int start, int end, boolean negativeStrand) {
		String key = (new SequenceMetadata(sequenceName, 0, negativeStrand, 0, 0)).getUniqueKey();
		SequenceMetadata seqM = realSequencesByKey.get(key);
		if(seqM==null) return null;
		//System.out.println("Internal i: "+seqM.getInternalIndex());
		FMIndexSingleSequence internalFM = internalIndexes.get(seqM.getInternalIndex());
		if(internalFM==null) throw new RuntimeException("Index not found for sequence: "+sequenceName);
		int internalSequenceStart = seqM.getInternalStart(); 
		int internalQueryStart = start;
		int internalQueryEnd = end;
		if(negativeStrand) {
			internalQueryEnd = seqM.getLength()-start;
			internalQueryStart = seqM.getLength()-end;
		}
		internalQueryStart += internalSequenceStart;
		internalQueryEnd += internalSequenceStart;
		return internalFM.getSequence(internalQueryStart, internalQueryEnd);
	}
	
	private ReadAlignment buildAlignmentFromMetadata(String searchSequence, int i, int internalPosMatch) 
	{
		List<SequenceMetadata> metadata = realSequencesMap.get(i);
		int internalPostStartSeq = 0;
		for(SequenceMetadata seqM:metadata) {
			if(internalPostStartSeq + seqM.getLength() > internalPosMatch) {
				String seqName = seqM.getSeqName();
				int first = internalPosMatch-internalPostStartSeq;
				int l = searchSequence.length();
				int last = first + l - 1;
				int flags = 0;
				if(seqM.isNegativeStrand()) {
					last = seqM.getLength()-1-first;
					first = last - l + 1 ;
					flags+=ReadAlignment.FLAG_READ_REVERSE_STRAND;
				}
				if(first < 0) return null;
				if(last>=seqM.getLength()) return null;
				return new ReadAlignment(seqName, first+1, last+1,l,flags);
			}
			internalPostStartSeq += seqM.getLength();
		}
		return null;
	}

	/*public List<GenomicRegion> search (String sequenceName, String searchSequence) 
	{
		List<GenomicRegion> alignments = new ArrayList<>();
		for (FMIndexSingleSequence idxSeq:singleSequenceIndexes) 
		{
			if(idxSeq.getSequenceName().equals(searchSequence))
			{
				return idxSeq.search(searchSequence);
			}
		}
		return alignments;
	}*/

	
	public static void main(String[] args) throws IOException
	{
		FMIndex f = new FMIndex();
		f.loadGenome(args[0]);
		//f.save(args[1]);
		//f.loadGenome(args[1]);
		System.out.println("Internal sequences: "+f.realSequencesByKey.size());
		SequenceMetadata seqM = f.realSequencesByKey.get("seq1_false");
		int l = seqM.getLength();
		System.out.println("Length: "+l);
		System.out.println(f.getSequence("seq1", 0, l, false));
		System.out.println(f.getSequence("seq1", 0, l, true));
		System.out.println("Positive 3-10: "+f.getSequence("seq1", 3, 10, false));
		System.out.println("Negative 3-10: "+f.getSequence("seq1", 3, 10, true));
		List<ReadAlignment> a = f.search("ata");
		for (int i = 0; i < a.size(); i++) {
			System.out.println(a.get(i).getSequenceName()+" pos:"+a.get(i).getFirst()+" to: "+a.get(i).getLast()+" flags:"+a.get(i).getFlags());
		}
	}
}
class SequenceMetadata implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -4208212205942777492L;
	private String seqName;
	private int length;
	private boolean negativeStrand;
	private int internalIndex;
	private int internalStart;
	public SequenceMetadata(String seqName, int length, boolean negativeStrand, int internalIndex, int internalStart) {
		super();
		this.seqName = seqName;
		this.length = length;
		this.negativeStrand = negativeStrand;
		this.internalIndex = internalIndex;
		this.internalStart = internalStart;
	}
	/**
	 * @return the seqName
	 */
	public String getSeqName() {
		return seqName;
	}
	/**
	 * @return the length
	 */
	public int getLength() {
		return length;
	}
	/**
	 * @return the negativeStrand
	 */
	public boolean isNegativeStrand() {
		return negativeStrand;
	}
	/**
	 * @return the internalIndex
	 */
	public int getInternalIndex() {
		return internalIndex;
	}
	/**
	 * @return the internalStart
	 */
	public int getInternalStart() {
		return internalStart;
	}
	
	public String getUniqueKey() {
		return seqName+"_"+negativeStrand;
	}
	
	
	
}
