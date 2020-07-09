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

import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.LimitedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;


/**
 * Implementation of a reference genome storing the sequences as DNAMaskedSequences to
 * minimize memory consumption 
 * @author Jorge Duitama
 */
public class ReferenceGenome { 
	private QualifiedSequenceList sequences;
	private String filename;
	/**
	 * Creates a new ReferenceGenome with the given data
	 * @param filename Name of the fasta file with the reference genome 
	 * @throws IOException If the file can not be read
	 */
	public ReferenceGenome (String filename) throws IOException {
		this.filename = filename;
		FastaSequencesHandler handler = new FastaSequencesHandler();
		handler.setSequenceType(DNAMaskedSequence.class);
		sequences = handler.loadSequences(filename);
		sequences.setAllowChanges(false);
	}
	/**
	 * Creates a reference genome sequence with the given sequence
	 * @param refQS Unique sequence in this reference. refQS!=null
	 */
	public ReferenceGenome(QualifiedSequence refQS) {
		sequences = new QualifiedSequenceList();
		sequences.add(refQS);
	}
	
	/**
	 * @return String the path of the file from which this genome was loaded
	 */
	public String getFilename() {
		return filename;
	}
	/**
	 * Returns the reference base pair at the given coordinate
	 * @param sequenceName Name of the sequence to search
	 * @param absolutePosition Position within the sequence
	 * @return char Basepair at the given coordinate. 0 if the coordinate is outside the genome
	 */
	public char getReferenceBase(String sequenceName, int absolutePosition) {
		QualifiedSequence qS = sequences.get(sequenceName);
		if(qS == null) return 0;
		CharSequence seq = qS.getCharacters();
		if(seq!=null && absolutePosition>=1 && absolutePosition<=seq.length()) {
			return seq.charAt(absolutePosition-1);
		}
		return 0;
	}
	/**
	 * Returns an unmodifiable list with objects storing the metadata of the sequences within the genome
	 * @return QualifiedSequenceList List of QualifiedSequence objects including the sequences metadata
	 */
	public QualifiedSequenceList getSequencesMetadata() {
		QualifiedSequenceList list = new QualifiedSequenceList();
		for(QualifiedSequence seqGenome:sequences) {
			QualifiedSequence seq = new QualifiedSequence(seqGenome.getName());
			seq.setComments(seqGenome.getComments());
			seq.setLength(seqGenome.getLength());
			list.add(seq);
		}
		list.setAllowChanges(false);
		return list;
	}
	/**
	 * Returns an unmodifiable list with objects storing the complete information of the sequences
	 * @return QualifiedSequenceList List of QualifiedSequence objects including the sequences in this genome
	 */
	public QualifiedSequenceList getSequencesList() {
		QualifiedSequenceList list = new QualifiedSequenceList();
		for(QualifiedSequence seqGenome:sequences) {
			QualifiedSequence seq = new QualifiedSequence(seqGenome.getName());
			seq.setComments(seqGenome.getComments());
			seq.setLength(seqGenome.getLength());
			seq.setCharacters(seqGenome.getCharacters());
			list.add(seq);
		}
		list.setAllowChanges(false);
		return list;
	}

	/**
	 * Prints this reference genome in the given stream
	 * @param out Stream to print the genome in fasta format
	 * @param lineLength Number of bases per line
	 */
	public void saveGenome(PrintStream out, int lineLength) {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		handler.saveSequences(sequences, out, lineLength);
	}

	/**
	 * Changes the reference base at the given genomic coordinate
	 * @param sequenceName Name of the sequence to look for
	 * @param absolutePosition One based position relative to the start of the sequence with
	 * the given name  
	 * @param base New base to set.
	 */
	public void setReferenceBase(String sequenceName, int absolutePosition, char base) {
		QualifiedSequence qS = sequences.get(sequenceName);
		if(qS == null) return;
		LimitedSequence seq = (LimitedSequence)qS.getCharacters();
		int pos = absolutePosition-1;
		if(seq==null || pos <0 || pos >= seq.length() || !seq.isInAlphabet(base)) return;
		char current = seq.charAt(pos);
		if (Character.toUpperCase(current) != Character.toUpperCase(base)) {
			if(Character.isUpperCase(current)) {
				base = Character.toUpperCase(base);
			} else {
				base = Character.toLowerCase(base);
			}
		}
		seq.setCharAt(pos, base);
	}
	/**
	 * Finds the sequence with the given name
	 * @param sequenceName Name to look for
	 * @return QualifiedSequence Object with the sequence information
	 */
	public QualifiedSequence getSequenceByName(String sequenceName) {
		return sequences.get(sequenceName);
	}
	/**
	 * Finds the sequence with the given index
	 * @param index Sequence index to look for
	 * @return QualifiedSequence Object with the sequence information
	 */
	public QualifiedSequence getSequenceByIndex(int index) {
		return sequences.get(index);
	}
	public CharSequence getSequenceCharacters (String sequenceName) {
		QualifiedSequence qs = sequences.get(sequenceName);
		if(qs==null) return null;
		return qs.getCharacters();
	}
	/**
	 * Returns the complete characters of the sequence with the given index
	 * @param index 0-based index to look for
	 * @return CharSequence sequence at the given index
	 */
	public CharSequence getSequenceCharacters (int index) {
		QualifiedSequence qs = sequences.get(index);
		if(qs==null) return null;
		return qs.getCharacters();
	}
	/**
	 * Returns the sequence at the given coordinates
	 * @param r Genomic region to search
	 * @return CharSequence object with the sequence within the given region
	 */
	public CharSequence getReference(GenomicRegion r) {
		return getReference(r.getSequenceName(),r.getFirst(),r.getLast());
	}
	/**
	 * Returns the sequence at the given coordinates
	 * @param sequenceName Name of the sequence to search
	 * @param first 1-based first genomic position within the sequence with the given name
	 * @param last 1-based last genomic position within the sequence with the given name
	 * @return CharSequence object with the sequence within the given region
	 */
	public CharSequence getReference(String sequenceName, int first, int last) {
		CharSequence seq = getSequenceCharacters(sequenceName);
		return getSubsequenceGenomicCoordinates(first, last, seq);
	}
	/**
	 * Returns the sequence at the given coordinates
	 * @param index 0-ased index in the sequences list
	 * @param first 1-based first genomic position within the sequence with the given name
	 * @param last 1-based last genomic position within the sequence with the given name
	 * @return CharSequence object with the sequence within the given region
	 */
	public CharSequence getReference(int index, int first, int last) {
		CharSequence seq = getSequenceCharacters(index);
		return getSubsequenceGenomicCoordinates(first, last, seq);
	}
	private CharSequence getSubsequenceGenomicCoordinates(int first, int last, CharSequence seq) {
		if(seq!=null && first >=1 && last <=seq.length()) {
			return seq.subSequence(first-1, last);
		}
		return null;
	}
	/**
	 * @return long total genome length
	 */
	public long getTotalLength() {
		return sequences.getTotalLength();
	}
	/**
	 * @return int the number of sequences in the genome
	 */
	public int getNumSequences() {
		return sequences.size();
	}
	/**
	 * @return List<String> List of names of the sequences in the genome
	 */
	public List<String> getSequenceNamesStringList() {
		return sequences.getNamesStringList();
	}

}
