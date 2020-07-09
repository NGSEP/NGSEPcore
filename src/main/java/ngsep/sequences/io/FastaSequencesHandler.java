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
package ngsep.sequences.io;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.List;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;



/**
 * General handler for sequences in fasta format. Implements methods to load and
 * write entire sequence sets
 * @author Jorge Duitama
 *
 */
public class FastaSequencesHandler {
	
	private Class<? extends CharSequence> sequenceType = null;
	
	public FastaSequencesHandler () {
		setSequenceType(DNAMaskedSequence.class);
	}
	/**
	 * Loads the sequences present in the given filename
	 * @param filename Name of the fasta file where sequences must be loaded 
	 * @throws IOException If the file can not be read
	 */
	public QualifiedSequenceList loadSequences(String filename) throws IOException {
		QualifiedSequenceList answer = new QualifiedSequenceList();
		try (FastaFileReader reader = new FastaFileReader(filename)) {
			reader.setSequenceType(sequenceType);
			Iterator<QualifiedSequence> it = reader.iterator();
			while(it.hasNext()) {
				answer.add(it.next());
			}
		}
		return answer;
	}
	/**
	 * Loads the sequences present in the given stream
	 * @param is stream to load sequences in fasta format 
	 * @throws IOException If the file can not be read
	 */
	public QualifiedSequenceList loadSequences(InputStream is) throws IOException {
		QualifiedSequenceList answer = new QualifiedSequenceList();
		try (FastaFileReader reader = new FastaFileReader(is)) {
			reader.setSequenceType(sequenceType);
			Iterator<QualifiedSequence> it = reader.iterator();
			while(it.hasNext()) {
				answer.add(it.next());
			}
		}
		return answer;
	}
	/**
	 * @return Class<?> Type of the sequences that will be loaded
	 */
	public Class<?> getSequenceType() {
		return sequenceType;
	}
	/**
	 * New type for sequences that will be loaded
	 * @param sequenceType New sequence type
	 */
	public void setSequenceType(Class<? extends CharSequence> sequenceType) {
		this.sequenceType = sequenceType;
	}

	/**
	 * Dump all sequences in the given print stream
	 * @param out Stream to print the sequences
	 * @param lineLength Number of bases per line
	 */
	public void saveSequences(List<QualifiedSequence> sequences, PrintStream out,int lineLength) {
		for(QualifiedSequence seq: sequences) {
			out.print(">");
			out.print(seq.getName());
			if(seq.getComments()!=null) {
				out.print(" ");
				out.print(seq.getComments());
			}
			out.println();
			CharSequence characters = seq.getCharacters();
			int l = characters.length();
			for(int j=0;j<l;j+=lineLength) {
				out.println(characters.subSequence(j, Math.min(l, j+lineLength)));
			}
		}
	}
}
