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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import ngsep.sequences.QualifiedSequenceList;

/**
 * @author Jorge Duitama
 */
public class SimpleSequenceListLoader {
	/**
	 * Loads a list of sequences from a file
	 * @param filename Name of the file to load the sequences
	 * @return QualifiedSequenceList List of sequence names. The list is set as unmodifiable by default
	 * @throws IOException If there is a problem reading the given file
	 */
	public QualifiedSequenceList loadSequences(String filename) throws IOException {
		QualifiedSequenceList answer = new QualifiedSequenceList();
		
		try (FileReader fr = new FileReader(filename);
			BufferedReader in = new BufferedReader(fr)) {
			String line=in.readLine();
			while(line!=null) {
				String [] items = line.split("\t| ");
				answer.addOrLookupName(items[0]);
				line=in.readLine();
			}
			
		}
		answer.setAllowChanges(false);
		return answer;
	}
}
