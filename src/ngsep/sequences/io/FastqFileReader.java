package ngsep.sequences.io;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.logging.Logger;

import ngsep.main.io.ConcatGZIPInputStream;
import ngsep.sequences.RawRead;

public class FastqFileReader implements Iterable<RawRead>,Closeable  {

	private Logger log = Logger.getLogger(FastqFileReader.class.getName());
	
	private BufferedReader in;
	
	private FastqFileIterator currentIterator = null;
	
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
			if(file.getName().endsWith(".gz")) {
				stream = new ConcatGZIPInputStream(stream);
			}
		}
		in = new BufferedReader(new InputStreamReader(stream));
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
					read = RawRead.load(in);
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
