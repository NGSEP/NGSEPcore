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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
/**
 * Script that sorts, orients and rename sequences in a fasta file 
 * @author Jorge Duitama
 *
 */
public class FastaSortAndRenameSequences {
	/**
	 * Sorts, orients and rename sequences in a fasta file, based on a text file describing the required changes.
	 * Writes a modified fasta file to the standard output
	 * @param args. Program arguments. args[0] is the fasta file to process. args[1] is a text file delimited by tab or space with the following fields:
	 * - Contig name
	 * - Orientation (+,-)
	 * - New name
	 * Sequences in this file should be sorted according to the desired order
	 * @throws Exception if the files can not be read
	 */
	public static void main(String[] args) throws Exception {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList sequences = new QualifiedSequenceList(handler.loadSequences(args[0]));
		String indexFile = args[1];
		Map<String,Boolean> changeOrientations = new LinkedHashMap<>();
		Map<String,String> changeNames = new LinkedHashMap<>();
		try (FileReader reader = new FileReader(indexFile);
			 BufferedReader in = new BufferedReader(reader)) {
			String line = in.readLine();
			while (line != null) {
				String [] items = line.split(" |\t");
				changeOrientations.put(items[0],"-".equals(items[1]));
				changeNames.put(items[0],items[2]);
				line = in.readLine();
			}
		}
		List<QualifiedSequence> sortedSeqs = new ArrayList<>();
		Set<String> sortedSeqIds = new HashSet<>();
		Set<String> sortedSeqNewIds = new HashSet<>();
		for(String seqName:changeOrientations.keySet()) {
			QualifiedSequence seq =  sequences.get(seqName);
			if(seq==null) {
				System.err.println("Sequence name: "+seqName+" not found in fasta");
				continue;
			}
			DNAMaskedSequence characters = (DNAMaskedSequence) seq.getCharacters();
			boolean rc = changeOrientations.get(seqName);
			String changedName = changeNames.get(seqName);
			if(rc) {
				characters = characters.getReverseComplement();
			}
			if(sortedSeqNewIds.contains(changedName)) {
				System.err.println("Duplicated new name : "+changedName);
				continue;
			}
			sortedSeqs.add(new QualifiedSequence(changedName,characters));
			sortedSeqIds.add(seqName);
			sortedSeqNewIds.add(changedName);
		}
		for(QualifiedSequence sequence:sequences) {
			if(!sortedSeqIds.contains(sequence.getName())) {
				System.err.println("Unplaced sequence "+sequence.getName()+" "+sequence.getLength());
				sortedSeqs.add(sequence);
			}
		}
		handler.saveSequences(sortedSeqs, System.out, 100);
	}

}
