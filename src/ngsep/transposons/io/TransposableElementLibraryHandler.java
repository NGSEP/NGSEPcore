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
package ngsep.transposons.io;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.transposons.TransposableElement;

/**
 * @author Jorge Duitama
 */
public class TransposableElementLibraryHandler {
	/**
	 * Load transposable elements from a given file
	 * @param filename Name of the file with the TE library
	 * @return List<TransposableElement> Transposable elements from the given file
	 * @throws IOException If the file can not be read
	 */
	public List<TransposableElement> load(String filename) throws IOException{
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> sequences = handler.loadSequences(filename);
		List<TransposableElement> tes = new ArrayList<TransposableElement>();
		for(QualifiedSequence seq:sequences) {
			tes.add(new TransposableElement(seq));
		}
		return tes;
	}
	
	public void save(List<TransposableElement> tes, String outputFile) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> sequences = new ArrayList<QualifiedSequence>();
		for(TransposableElement te:tes) {
			QualifiedSequence seq = new QualifiedSequence(te.getId(),te.getSequence());
			sequences.add(seq);
		}
		try(PrintStream out = new PrintStream(outputFile)) {
			handler.saveSequences(sequences, out, 100);
		}
	}
}
