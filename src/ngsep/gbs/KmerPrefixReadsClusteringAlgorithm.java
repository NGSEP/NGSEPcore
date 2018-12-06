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
package ngsep.gbs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNASequence;
import ngsep.sequences.DNAShortKmer;
import ngsep.sequences.DNAShortKmerClusterMap;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;
import ngsep.sequencing.ReadsDemultiplex;

/**
 * @author Jorge Gomez
 * @author Jorge Duitama
 */
public class KmerPrefixReadsClusteringAlgorithm {

	private Logger log = Logger.getLogger(ReadsDemultiplex.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	public static final int DEF_KMER_LENGTH = 31;
	public static final int DEF_START = 8;
	public static final String DEF_REGEXP_SINGLE="<S>.fastq.gz";
	public static final String DEF_REGEXP_PAIRED="<S>_<N>.fastq.gz";
	
	private String inputDirectory=".";
	private String outputDirectory=".";
	private int kmerLength = DEF_KMER_LENGTH;
	private Pattern regexp=Pattern.compile(DEF_REGEXP_SINGLE);
	private Map<String, String> filenamesBySampleId1=new HashMap<>();
	private Map<String, String> filenamesBySampleId2=new HashMap<>();
	private DNAShortKmerClusterMap kmersMap;
	
	
	public static void main(String[] args) throws Exception {
		KmerPrefixReadsClusteringAlgorithm instance = new KmerPrefixReadsClusteringAlgorithm();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.inputDirectory = args[i++];
		instance.outputDirectory = args[i++];
		instance.run();
	}
	
	public void setProgressNotifier(ProgressNotifier progressNotifier) { 
		this.progressNotifier = progressNotifier;
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	/**
	 * @return the kmerLength
	 */
	public int getKmerLength() {
		return kmerLength;
	}
	/**
	 * @param kmerLength the kmerLength to set
	 */
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
	}
	
	public void setKmerLength(String value) {
		setKmerLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	public void run() throws IOException {
		loadFilenamesAndSamples();
		log.info("Loaded "+filenamesBySampleId1.size()+" samples");
		buildKmersMap();
		log.info("Built kmers map with "+kmersMap.size()+" clusters");
		clusterReads();
		log.info("Clustered reads");
		callVariants();
		log.info("Called variants");
		printStatistics();
		log.info("Process finished");
	}

	private void loadFilenamesAndSamples() {
		File[] files = (new File(inputDirectory)).listFiles();
		for(File f : files) {
			String filename = f.getName();
			//TODO: Use pattern
			int i = filename.indexOf(".fastq");
			if(i>=0) {
				String sampleId = filename.substring(0, i);
				filenamesBySampleId1.put(sampleId, f.getAbsolutePath());
			}
		}
	}

	public void buildKmersMap() throws IOException {
		kmersMap = new DNAShortKmerClusterMap();
		for(String filename:filenamesBySampleId1.values()) {
			addKmersFromFile(filename);
		}
	}

	private void addKmersFromFile(String filename) throws IOException {
		int readCount = 0;
		try (FastqFileReader openFile = new FastqFileReader(filename);) {
			Iterator<RawRead> reader = openFile.iterator();
			while(reader.hasNext()) {
				RawRead read = reader.next();
				String s = read.getSequenceString();
				if(DEF_START + kmerLength>s.length()) continue;
				String prefix = s.substring(DEF_START,DEF_START + kmerLength);
				if(DNASequence.isDNA(prefix)) {
					kmersMap.addOcurrance(new DNAShortKmer(prefix));
					readCount++;
				}
			}
		}
		log.info("Processed a total of " + readCount + " reads for file: "+filename);
		
	}
	public void clusterReads() throws IOException {
		ClusteredReadsCache clusteredReadsCache = new ClusteredReadsCache(); 
		for(String sampleId:filenamesBySampleId1.keySet()) {
			String filename1 = filenamesBySampleId1.get(sampleId);
			String filename2 = filenamesBySampleId2.get(sampleId);
			if(filename2 == null) {
				clusterReadsSingleFile (filename1, clusteredReadsCache);
			} else {
				clusterReadsPairedEndFiles (filename1, filename2, clusteredReadsCache);
			}
		}
		
	}

	private void clusterReadsSingleFile(String filename, ClusteredReadsCache clusteredReadsCache) throws IOException {
		try (FastqFileReader openFile = new FastqFileReader(filename);) {
			Iterator<RawRead> reader = openFile.iterator();
			while(reader.hasNext()) {
				RawRead read = reader.next();
				String s = read.getSequenceString();
				if(DEF_START + kmerLength>s.length()) continue;
				String prefix = s.substring(DEF_START,DEF_START + kmerLength);
				if(!DNASequence.isDNA(prefix)) continue;
				Integer k = kmersMap.getCluster(new DNAShortKmer(prefix));
				if(k==null) continue;
				clusteredReadsCache.addSingleRead(k, read);
				//TODO: Dump cache if needed
			}
		}
	}




	private void clusterReadsPairedEndFiles(String filename1, String filename2, ClusteredReadsCache clusteredReadsCache) throws IOException {
		// TODO Auto-generated method stub
		
	}




	public void callVariants() throws IOException {
		// TODO Implement. Load in parallel clustered read files. Align reads within clusters build pileups and call variants
		
	}
	private void printStatistics() {
		// TODO Implement. Save to the output directory a file with process statistics
		
	}
	
}
class ClusteredReadsCache {
	private Map<Integer,List<RawRead>> clusteredReadsCache = new TreeMap<>();
	private int totalReads = 0;
	public void addSingleRead(int k, RawRead read) {
		List<RawRead> readsClusterK = clusteredReadsCache.get(k);
		if(readsClusterK==null) {
			readsClusterK = new ArrayList<>();
			clusteredReadsCache.put(k, readsClusterK);
		}
		readsClusterK.add(read);
	}
	/**
	 * @return the totalReads
	 */
	public int getTotalReads() {
		return totalReads;
	}
	/**
	 * Dumps the cache to the file wit the given name and clears this cache
	 * @param filename to dump the cache
	 */
	public void dump(String filename) {
		//TODO: Implement
	}
}
