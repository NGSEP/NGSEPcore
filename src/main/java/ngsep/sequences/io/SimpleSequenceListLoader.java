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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

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
		FileInputStream fis = null;
		BufferedReader in = null;
		try {
			fis = new FileInputStream(filename);
			//Open buffer
			in = new BufferedReader(new InputStreamReader (fis));
			String line=in.readLine();
			while(line!=null) {
				String [] items = line.split("\t| ");
				answer.addOrLookupName(items[0]);
				line=in.readLine();
			}
			
		} finally {
			if(in!=null) in.close();
			if(fis!=null) fis.close();
		}
		answer.setAllowChanges(false);
		return answer;
	}
}
