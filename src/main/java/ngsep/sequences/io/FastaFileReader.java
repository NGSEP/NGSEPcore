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
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.logging.Logger;

import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class FastaFileReader implements Iterable<QualifiedSequence>,Closeable  {
	
	private Logger log = Logger.getLogger(FastaFileReader.class.getName());
	
	/**
	 * Constant to load the read id, read sequence and quality information 
	 */
	public static final int LOAD_MODE_FULL = 0;
	/**
	 * Constant to load only the read sequence of each read
	 */
	public static final int LOAD_MODE_MINIMAL = 1;
	
	private BufferedReader in;
	
	private String currentLine = "";
	
	private FastaFileIterator currentIterator = null;
	
	private int loadMode = LOAD_MODE_FULL;
	
	private boolean keepLowerCase = false;
	
	private Class<? extends CharSequence> sequenceType = null;
	private Constructor<? extends CharSequence> emptyConstructor = null;
	private Method appendMethod = null;
	
	public FastaFileReader (String filename) throws IOException {
		init(null,new File(filename));
	}
	public FastaFileReader (File file) throws IOException {
		init(null,file);
	}
	public FastaFileReader (InputStream stream) throws IOException {
		init(stream,null);
	}
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		if (log == null) throw new NullPointerException("Log can not be null");
		this.log = log;
	}
	
	public int getLoadMode() {
		return loadMode;
	}
	public void setLoadMode(int loadMode) {
		this.loadMode = loadMode;
	}
	
	public boolean isKeepLowerCase() {
		return keepLowerCase;
	}
	/**
	 * Changes the behavior to keep lowercase characters if they exist. By default all characters are converted to upper case
	 * @param keepLowerCase
	 */
	public void setKeepLowerCase(boolean keepLowerCase) {
		this.keepLowerCase = keepLowerCase;
	}
	/**
	 * @return Class datatype for sequences to load
	 */
	public Class<? extends CharSequence> getSequenceType() {
		return sequenceType;
	}
	/**
	 * Changes the datatype for the sequences
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
			appendMethod = sequenceType.getMethod("append", CharSequence.class);
		} catch (NoSuchMethodException | SecurityException e) {
			appendMethod = null;
		}
		if(emptyConstructor == null || appendMethod == null) throw new RuntimeException("Class for sequences must have an empty constructor plus an append method");
	}
	@Override
	public void close() throws IOException {
		in.close();
	}

	@Override
	public Iterator<QualifiedSequence> iterator() {
		if (in == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (currentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        currentIterator = new FastaFileIterator();
		return currentIterator;
	}
	
	private void init (InputStream stream, File file) throws IOException {
		if (stream != null && file != null) throw new IllegalArgumentException("Stream and file are mutually exclusive");
		if(file!=null) {
			stream = new FileInputStream(file);
			if(file.getName().toLowerCase().endsWith(".gz")) {
				stream = new ConcatGZIPInputStream(stream);
			}
		}
		in = new BufferedReader(new InputStreamReader(stream));
		//Sequence name by default
		setSequenceType(DNAMaskedSequence.class);
	}
	/**
	 * Loads a raw read from the given BufferedReader
	 * @param in buffer to read
	 * @return RawRead
	 * @throws IOException if the buffer can not be read
	 */
	private QualifiedSequence load (BufferedReader in) throws IOException {
		// end of file
		if (currentLine == null) return null;
		//Read sequences
		String id = null;
		String comments = null;
		if (currentLine.length()>0) {
			int i = currentLine.indexOf(' ');
			int i2 = currentLine.indexOf('\t');
			if (i<=1) i = i2;
			if (i2>1 && i2<i) i = i2;
			if(i<=1) {
				id = currentLine.substring(1);
			} else {
				id = currentLine.substring(1, i);
				if(i+1<currentLine.length()) {
					comments = currentLine.substring(i+1);
				}
			}
		}
		StringBuilder buffer = new StringBuilder();
		CharSequence nextSequence = createEmptySequence();
		currentLine = in.readLine();
		while(currentLine!=null) {
			currentLine = currentLine.trim();
			if (currentLine.length()==0) {
				currentLine = in.readLine();
				continue;
			}
			char firstChr = currentLine.charAt(0);
			if(firstChr=='>') {
				break;
			} else if (firstChr!='#') {
				String procesedLine = removeSpaces(currentLine);
				buffer.append(procesedLine);
				if(buffer.length() > 10000000) {
					try {
						appendMethod.invoke(nextSequence, buffer.toString());
					} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
						throw new IOException(e);
					}
					buffer.delete(0, buffer.length());
				}
			}
			currentLine = in.readLine();
		}
		if (buffer.length()>0) {
			try {
				appendMethod.invoke(nextSequence, buffer.toString());
			} catch (IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
				throw new IOException(e);
			}
		}
		if (id==null) {
			if (nextSequence.length()>0) log.warning("Sequence with length "+nextSequence.length()+" found before the first id sequence");
			return new QualifiedSequence(null,nextSequence);
		}
		if (nextSequence.length()==0) log.warning("Loaded empty sequence with id: "+id);
		QualifiedSequence seq = new QualifiedSequence(id,nextSequence);
		if(comments!=null) seq.setComments(comments);
		return seq;
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
	
	private boolean passFilters (QualifiedSequence record) {
		//TODO: Implement filters
		return true;
	}
	
	private class FastaFileIterator implements Iterator<QualifiedSequence> {
		private QualifiedSequence nextRecord;
		public FastaFileIterator() {
			nextRecord = loadRecord();
		}
		@Override
		public boolean hasNext() {
			return nextRecord!=null;
		}

		@Override
		public QualifiedSequence next() {
			if(nextRecord==null) throw new NoSuchElementException();
			QualifiedSequence answer = nextRecord;
			nextRecord = loadRecord();
			return answer;
		}

		private QualifiedSequence loadRecord() {
			QualifiedSequence record;
			while(true) {
				try {
					record = load(in);
					//Case to reach the first record
					if(record !=null && record.getName()==null) record = load(in);
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
				if(record==null) return null;
				if(passFilters(record)) return record;
			} 
		}
		@Override
		public void remove() {
			throw new UnsupportedOperationException("Remove not supported by FastqFileIterator");
		}
	}
}
