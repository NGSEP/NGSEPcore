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
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
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
	Class<? extends CharSequence> sequenceType = null;
	private Constructor<? extends CharSequence> emptyConstructor = null;
	private Constructor<? extends CharSequence> stringConstructor = null;
	private Method appendMethod = null;
	
	
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
		FileInputStream fis = null;
		BufferedReader in = null;
		try {
			fis = new FileInputStream(filename);
			//Open buffer
			in = new BufferedReader(new InputStreamReader (fis));
			//Read sequences
			String id =null;
			String comment = null;
			String line=in.readLine();
			
			StringBuilder buffer = new StringBuilder();
			CharSequence nextSequence = createEmptySequence();
			while(line!=null) {
				if(line.startsWith(">")) {
					addSequence(answer, id, comment, nextSequence, buffer);
					nextSequence = createEmptySequence();
					String idLine = line.substring(1);
					String [] items = idLine.split(" |\t");
					id = items[0];
					if(id.length()+1<idLine.length()) {
						comment = idLine.substring(id.length()+1);
					} else {
						comment = null;
					}	
				} else if (!line.startsWith("#")) {
					line = removeSpaces(line);
					buffer.append(line);
					if(nextSequence!=null && appendMethod!=null && buffer.length() > 10000000) {
						try {
							appendMethod.invoke(nextSequence, buffer.toString());
						} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
							throw new IOException(e);
						}
						buffer.delete(0, buffer.length());
					}
				}
				line=in.readLine();
			}
			if(id!=null) {
				//System.out.println("Adding sequence "+id+". Length: "+nextSequence.length());
				addSequence(answer,id, comment, nextSequence, buffer);
			}
		} finally {
			if(in!=null) in.close();
			if(fis!=null) fis.close();
		}
		return answer;
	}
	private CharSequence createEmptySequence() throws IOException {
		CharSequence nextSequence = null;
		if(emptyConstructor!=null && appendMethod !=null) {	
			try {
				nextSequence = emptyConstructor.newInstance();
			} catch (Exception e) {
				throw new IOException(e);
			}
		}
		return nextSequence;
	}
	private String removeSpaces(String line) {
		//return line.replaceAll("\\s", "");
		StringBuilder answer = new StringBuilder();
		for(int i=0;i<line.length();i++) {
			char c = line.charAt(i);
			if(!Character.isISOControl(c)) {
				answer.append(c);
			}
		}
		return answer.toString();
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
		try {
			emptyConstructor = (Constructor<? extends CharSequence>)sequenceType.getConstructor();
		} catch (NoSuchMethodException | SecurityException e1) {
			emptyConstructor = null;
		}
		try {
			stringConstructor = (Constructor<? extends CharSequence>)sequenceType.getConstructor(CharSequence.class);
		} catch (NoSuchMethodException | SecurityException e1) {
			stringConstructor = null;
		}
		try {
			appendMethod = sequenceType.getMethod("append", CharSequence.class);
		} catch (NoSuchMethodException | SecurityException e) {
			appendMethod = null;
		}
		if(stringConstructor == null && (emptyConstructor == null || appendMethod == null)) throw new RuntimeException("Class for sequences must have an empty constructor plus an append method or a string constructor");
	}
	/**
	 * Adds the sequence with the given information to the given list
	 * @param sequences
	 * @param name
	 * @param comments
	 * @param sequence Sequence to be added if append is supported
	 * @param buffer with the end of the sequence or the complete sequence if append is not supoported
	 * @throws IOException 
	 */
	private void addSequence(QualifiedSequenceList sequences, String name, String comments, CharSequence sequence, StringBuilder buffer) throws IOException {
		CharSequence finalSequence;
		try {
			if(sequence!=null && sequence.length()>0) {
				appendMethod.invoke(sequence, buffer.toString());
				finalSequence = sequence;
			} else {
				finalSequence = stringConstructor.newInstance(buffer.toString());
			}
		} catch (Exception e) {
			throw new IOException(e);
		}
		if(name!=null && finalSequence.length()>0) {
			//System.out.println("Adding sequence "+id+". Length: "+nextSequence.length());
			QualifiedSequence seq = new QualifiedSequence(name,finalSequence);
			if(comments !=null) seq.setComments(comments);
			sequences.add(seq);
		}
		buffer.delete(0, buffer.length());
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
