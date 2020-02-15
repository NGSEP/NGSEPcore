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
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.logging.Logger;

import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.sequences.RawRead;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class FastqFileReader implements Iterable<RawRead>,Closeable  {
	
	private Logger log = Logger.getLogger(FastqFileReader.class.getName());
	
	/**
	 * Constant to load the read id, read sequence and quality information 
	 */
	public static final int LOAD_MODE_FULL = 0;
	/**
	 * Constant to load the read sequence and quality information of each read
	 */
	public static final int LOAD_MODE_QUALITY = 1;
	/**
	 * Constant to load only the read sequence of each read
	 */
	public static final int LOAD_MODE_MINIMAL = 2;
	
	private BufferedReader in;
	
	private FastqFileIterator currentIterator = null;
	
	private int loadMode = LOAD_MODE_FULL;
	
	private boolean keepLowerCase = false;
	
	private Class<? extends CharSequence> sequenceType = null;
	
	private Constructor<? extends CharSequence> charSequenceConstructor = null;
	
	public FastqFileReader (String filename) throws IOException {
		init(null,new File(filename));
	}
	public FastqFileReader (File file) throws IOException {
		init(null,file);
	}
	public FastqFileReader (InputStream stream) throws IOException {
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
			charSequenceConstructor = (Constructor<? extends CharSequence>)sequenceType.getConstructor(CharSequence.class);
		} catch (NoSuchMethodException | SecurityException e1) {
			throw new RuntimeException("The given sequence type does not have a string constructor", e1);
		}
	}
	@Override
	public void close() throws IOException {
		in.close();
	}

	@Override
	public Iterator<RawRead> iterator() {
		if (in == null) {
            throw new IllegalStateException("File reader is closed");
        }
        if (currentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
        currentIterator = new FastqFileIterator();
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
	}
	/**
	 * Loads a raw read from the given BufferedReader
	 * @param in buffer to read
	 * @return RawRead
	 * @throws IOException if the buffer can not be read
	 */
	private RawRead load (BufferedReader in) throws IOException {
		String id = in.readLine();
		if(id==null) return null;
		int i = id.indexOf(" ");
		if(i>0) id = id.substring(0,i);
		CharSequence seq = in.readLine();
		if(seq==null) return null;
		if(!keepLowerCase) seq = seq.toString().toUpperCase();
		String plus = in.readLine();
		if(plus==null) return null;
		String qs = in.readLine();
		if(qs==null) return null;
		
		if(sequenceType!=null) {
			try {
				seq = charSequenceConstructor.newInstance(seq);
			} catch (Exception e) {
				throw new RuntimeException("Can not invoke char sequence constructor",e);
			}
		}
		
		if(loadMode == LOAD_MODE_MINIMAL) return new RawRead(null, seq, null);
		else if (loadMode == LOAD_MODE_QUALITY) return new RawRead(null, seq, qs);
		else return new RawRead(id.substring(1), seq, qs);
	}
	
	private boolean passFilters (RawRead read) {
		//TODO: Implement filters
		return true;
	}
	
	private class FastqFileIterator implements Iterator<RawRead> {
		private RawRead nextRecord;
		public FastqFileIterator() {
			nextRecord = loadRecord();
		}
		@Override
		public boolean hasNext() {
			return nextRecord!=null;
		}

		@Override
		public RawRead next() {
			if(nextRecord==null) throw new NoSuchElementException();
			RawRead answer = nextRecord;
			nextRecord = loadRecord();
			return answer;
		}

		private RawRead loadRecord() {
			RawRead read;
			while(true) {
				try {
					read = load(in);
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
				if(read==null) return null;
				if(passFilters(read)) return read;
			} 
		}
		@Override
		public void remove() {
			throw new UnsupportedOperationException("Remove not supported by FastqFileIterator");
		}
	}
}
