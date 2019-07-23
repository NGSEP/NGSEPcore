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
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import ngsep.alignments.ReadAlignment;
import ngsep.discovery.PileupRecord;
import ngsep.discovery.VariantPileupListener;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNASequence;
import ngsep.sequences.DNAShortKmer;
import ngsep.sequences.DNAShortKmerClusterMap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;
import ngsep.sequencing.ReadsDemultiplex;
import ngsep.variants.CalledGenomicVariant;

/**
 * @author Jorge Gomez
 * @author Jorge Duitama
 * @author Andrea Parra
 */
public class KmerPrefixReadsClusteringAlgorithm {

	private Logger log = Logger.getLogger(ReadsDemultiplex.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	public static final int DEF_KMER_LENGTH = 31;
	public static final int DEF_START = 8;
	public static final int DEF_MAX_READS_IN_MEMORY = 1000000;
	public static final String DEF_REGEXP_SINGLE="<S>.fastq.gz";
	public static final String DEF_REGEXP_PAIRED="<S>_<N>.fastq.gz";
	
	private String inputDirectory=".";
	private String outPrefix="./output";
	private int kmerLength = DEF_KMER_LENGTH;
	private Pattern regexp=Pattern.compile(DEF_REGEXP_SINGLE);
	private Map<String, String> filenamesBySampleId1=new HashMap<>();
	private Map<String, String> filenamesBySampleId2=new HashMap<>();
	private DNAShortKmerClusterMap kmersMap;
	
	private ProcessInfo processInfo = new ProcessInfo();
	private int numClustersWithVariants = 0;
	private int numClusteredFiles;
	private int numUnclusteredReads = 0;
	private int numTotalReads = 0;
	
	
	public static void main(String[] args) throws Exception {
		KmerPrefixReadsClusteringAlgorithm instance = new KmerPrefixReadsClusteringAlgorithm();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.inputDirectory = args[i++];
		instance.outPrefix = args[i++];
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
		
		processInfo.addTime(System.nanoTime(), "Load files start");
		loadFilenamesAndSamples();
		processInfo.addTime(System.nanoTime(), "Load files end");
		processInfo.addTime(System.nanoTime(), "BuildKmersMap start");
		log.info("Loaded "+filenamesBySampleId1.size()+" samples");
		buildKmersMap();
		processInfo.addTime(System.nanoTime(), "BuildKmersMap end");
		processInfo.addTime(System.nanoTime(), "Cluster reads start");
		log.info("Built kmers map with "+kmersMap.size()+" clusters");
		List<String> clusteredReadsFilenames = clusterReads();
		processInfo.addTime(System.nanoTime(), "Cluster reads end");
		processInfo.addTime(System.nanoTime(), "Variant calling start");
		//List<String> clusteredReadsFilenames = debug();
		this.numClusteredFiles = clusteredReadsFilenames.size();
		log.info("Clustered reads");
		callVariants(clusteredReadsFilenames);
		processInfo.addTime(System.nanoTime(), "Variant calling end");
		log.info("Called variants");
		printStatistics();
		log.info("Process finished");
	}
	
	private List<String> debug() {		
		log.info("Skipping to call variants");
		List<String> clusteredReadsFilenames = new ArrayList<>();
		clusteredReadsFilenames.add("trial_debug_clusteredReads_0.fastq.gz");
		clusteredReadsFilenames.add("trial_debug_clusteredReads_1.fastq.gz");
		return clusteredReadsFilenames;
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
	public List<String> clusterReads() throws IOException {
		ClusteredReadsCache clusteredReadsCache = new ClusteredReadsCache(); 
		for(String sampleId:filenamesBySampleId1.keySet()) {
			String filename1 = filenamesBySampleId1.get(sampleId);
			String filename2 = filenamesBySampleId2.get(sampleId);
			if(filename2 == null) {
				log.info("Clustering reads from " + filename1);
				clusterReadsSingleFile (sampleId, filename1, clusteredReadsCache);
			} else {
				clusterReadsPairedEndFiles (sampleId, filename1, filename2, clusteredReadsCache);
			}
		}
		clusteredReadsCache.dump(outPrefix);
		return clusteredReadsCache.getClusteredReadFiles();
	}

	private void clusterReadsSingleFile(String sampleId, String filename, ClusteredReadsCache clusteredReadsCache) throws IOException {
		int unmatchedReads = 0;
		int count = 1;
		try (FastqFileReader openFile = new FastqFileReader(filename);) {
			Iterator<RawRead> reader = openFile.iterator();
			while(reader.hasNext()) {
				this.numTotalReads++;
				RawRead read = reader.next();
				String s = read.getSequenceString();
				if(DEF_START + kmerLength>s.length()) continue;
				String prefix = s.substring(DEF_START,DEF_START + kmerLength);
				if(!DNASequence.isDNA(prefix)) continue;
				Integer clusterId = kmersMap.getCluster(new DNAShortKmer(prefix));
				if(clusterId==null) {
					this.numUnclusteredReads++;
					unmatchedReads++;
					continue;
				}
				clusteredReadsCache.addSingleRead(clusterId, new RawRead(sampleId+"_"+clusterId+"_"+read.getName(), s, read.getQualityScores()));
				if(clusteredReadsCache.getTotalReads()>=DEF_MAX_READS_IN_MEMORY) {
					log.info("dumping reads");
					clusteredReadsCache.dump(outPrefix);
				}
				count++;
			}
			log.info(Integer.toString(unmatchedReads) + " reads remained unmatched for file: " + filename);
			log.info(Integer.toString(count) + " reads were succesfully matched for file: " + filename);
		}
	}

	private void clusterReadsPairedEndFiles(String sampleId, String filename1, String filename2, ClusteredReadsCache clusteredReadsCache) throws IOException {
		// TODO Auto-generated method stub
		return;
	}
	

	public void callVariants(List<String> clusteredReadsFilenames) throws IOException {
		int numberOfFiles = clusteredReadsFilenames.size();
		
		//process files in parallel
		FastqFileReader [] readers = new FastqFileReader[numberOfFiles];
		RawRead [] currentReads = new RawRead[numberOfFiles];
		Arrays.fill(currentReads, null);
		List<Iterator<RawRead>> iterators = new ArrayList<>();
		try (PrintStream outVariants = new PrintStream(outPrefix+"_variants.txt");) {
			int numNotNull = 0;
			int numCluster = 0;
			for(int i=0; i<numberOfFiles; i++) {
				readers[i] = new FastqFileReader(clusteredReadsFilenames.get(i));
				Iterator<RawRead> it = readers[i].iterator();
				iterators.add(it);
				if(it.hasNext()) {
					currentReads[i] = it.next();
					numNotNull++;
				}
			}
			log.info("Processing a total of " + Integer.toString(numberOfFiles) + " clustered files.");
			outVariants.println("clusterId \t position \t alleles \t variantQS");
			while(numNotNull>0) {
				//gather reads next cluster
				ReadCluster nextCluster = new ReadCluster(numCluster);
				for(int i=0; i<numberOfFiles; i++) {
					addReadsToCluster(nextCluster, iterators.get(i), currentReads, i);
					if(currentReads[i]==null) numNotNull--;
				}
		
				List<CalledGenomicVariant> variants = processCluster(nextCluster);
				for(CalledGenomicVariant variant:variants) {
					String[] alleles = variant.getCalledAlleles();
					String alls = "";
					for(String allel:alleles) {
						alls += allel + " ";
					}
					outVariants.println(
							variant.getSampleId() + "\t" +
							Integer.toString(variant.getFirst()) + "\t" +
							alls + "\t" +
							variant.getVariantQS());	
				}
				
				numCluster++;
			}
			
		} finally {
			for(FastqFileReader reader:readers) {
				if(reader!=null) reader.close();
			}
		}
	}
		
	
	private void addReadsToCluster(ReadCluster nextCluster, Iterator<RawRead> iterator, 
			RawRead[] currentReads, int i) throws IOException {
		
		RawRead currentRead = currentReads[i];
		int readsAdded = 0;
		int numCluster = nextCluster.getClusterNumber();
		while(currentRead!=null) {
			String readIdWithCluster = currentRead.getName();
			String [] items = readIdWithCluster.split("_");
			if(items.length < 2) {
				continue;
			}
			String sampleId = items[0];
			int currentReadCluster = Integer.parseInt(items[1]);
			String readId = items[2];
			//currentRead.setName(readId);
			if (currentReadCluster>numCluster) break;
			else if (currentReadCluster<numCluster) throw new RuntimeException("Disorgainzed file. Current cluster: "+numCluster+" found: "+currentReadCluster+" in read. "+readIdWithCluster);
			nextCluster.addRead(currentRead, sampleId);
			readsAdded++;
			if(iterator.hasNext()) {
				currentReads[i] = iterator.next();
				currentRead = currentReads[i];
				
			} else {
				log.info("Done with file " + Integer.toString(i) + ".");
				currentReads[i] = null;
				currentRead = null;
				
			}
		}
	}

	private List<CalledGenomicVariant> processCluster(ReadCluster readCluster) {
		List<CalledGenomicVariant> variants = new ArrayList<>();
		List<ReadAlignment> readAlignments = new ArrayList<>();
		
		int clusterId = readCluster.getClusterNumber();
		String refSeq = readCluster.getRefSeq();
		String referenceId = Integer.toString(clusterId);
		QualifiedSequence refQS = new QualifiedSequence(referenceId, refSeq);
		ReferenceGenome singleSequenceGenome = new ReferenceGenome (refQS);
		VariantPileupListener variantsDetector = new VariantPileupListener();
		variantsDetector.setGenome(singleSequenceGenome);
		
		// For each read within the cluster create a ReadAlignment. Set characters and quality scores
		for(RawRead read:readCluster.getReads()) {
			int readLength = read.getLength();
			String CIGARString = Integer.toString(readLength) + "M"; 
			ReadAlignment readAlignment = new ReadAlignment(referenceId, 1, readLength, readLength, 0);
			readAlignment.setQualityScores(read.getQualityScores());
			readAlignment.setReadCharacters(read.getCharacters());
			readAlignment.setReadName(read.getName());
			readAlignment.setCigarString(CIGARString);
			readAlignments.add(readAlignment);	
		}
		
		// For each position in the representative sequence create a pileup record with cluster id as sequence name and position =i
		// Start at 1? if not, no readAlignments are added to the pileup CHECK.
		for(int i=0; i<refSeq.length(); i++) {
			PileupRecord clusterPileUp = new PileupRecord(referenceId, i);
			//  Add the alignments to the pileup record
			for(ReadAlignment readAlgn:readAlignments) {
				clusterPileUp.addAlignment(readAlgn);
			}
			//  Use VariantPileuipListener to discover variants from the pileup record for the discovery step variant=null
			CalledGenomicVariant variant = variantsDetector.processPileup(clusterPileUp, null);
						
			if(variant!=null) {
				variant.setSampleId(referenceId);
				variants.add(variant);
			}
		}
		if (variants.size() > 0) {
			this.numClustersWithVariants++;
		}
		return variants;
	}
	
	private void printStatistics() throws IOException {
		// TODO Implement. Create an output file with process statistics
		try(PrintStream processStats = new PrintStream(outPrefix+"_processInfo.txt");){
			processStats.println("Process times:");
			processStats.print(this.processInfo.getTimes());
			processStats.println("Number of Files: " + Integer.toString(this.filenamesBySampleId1.size()));
			processStats.println("Number of Cluster Files: " + Integer.toString(this.numClusteredFiles));
			processStats.println("Number of Clusters: " + Integer.toString(kmersMap.size()));
			processStats.println("Number of Clusters with variants: " + Integer.toString(this.numClustersWithVariants));
			processStats.println("Number of Reads: " + Integer.toString(this.numTotalReads));
			processStats.println("Number of Unclustered Reads: " + Integer.toString(this.numUnclusteredReads));
		}	
	}
	
	/**
	 * 
	 */
	public void readClustersFromFile(String filename) throws IOException {
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
}

class ProcessInfo {
	private List<Long> processTimes = new ArrayList<>(); 
	private List<String> processTimesLabels = new ArrayList<>();
	
	public void addTime(long time, String label) {
		processTimes.add(time);
		processTimesLabels.add(label);
	}
	
	public String getTimes() {
		String timeInfo = "";
		for(int i=0; i<processTimes.size(); i++) {
			timeInfo += processTimesLabels.get(i) + ": " + Long.toString(processTimes.get(i)) + "\n";
		}
		return timeInfo;
	}
}

class ClusteredReadsCache {
	private Map<Integer,List<RawRead>> clusteredReadsCache = new TreeMap<>();
	private Map<Integer, Map<String, String>> clusterInfo = new TreeMap<>();
	private int totalReads = 0;
	private List<String> outFiles = new ArrayList<>();
	
	public void addSingleRead(int k, RawRead read) {
		List<RawRead> readsClusterK = clusteredReadsCache.get(k);
		if(readsClusterK==null) {
			readsClusterK = new ArrayList<>();
			clusteredReadsCache.put(k, readsClusterK);
		}
		readsClusterK.add(read);
		totalReads++;
	}
	
	public void addClusterInfo(int k, String key, String value) {
		Map<String, String> info = clusterInfo.get(k);
		if(clusterInfo == null) {
			info = new HashMap<>();		
		}
		info.put(key, value);
	}
	
	public Map<String, String> getClusterInfo(int k) {
		return clusterInfo.get(k);
	}
	
	public List<String> getClusteredReadFiles() {
		return outFiles;
	}
	
	public List<RawRead> getClusterReads(int k) {
		return clusteredReadsCache.get(k);
	}
	
	/**
	 * @return the totalReads
	 */
	public int getTotalReads() {
		return totalReads;
	}
	
	/**
	 * Dumps the cache to the file with the given name and clears this cache
	 * @param outPrefix prefix of the file to dump the cache
	 */
	public void dump(String outPrefix) throws IOException {
		int number = outFiles.size();
		String singleFilename = outPrefix+"_clusteredReads_"+number+".fastq.gz";
		outFiles.add(singleFilename);
		try (OutputStream os = new FileOutputStream(singleFilename);
			 GZIPOutputStream gos = new GZIPOutputStream(os);
			 PrintStream out = new PrintStream(gos);) {
			for(List<RawRead> readsCluster:clusteredReadsCache.values()) {
				for(RawRead read:readsCluster) {
					read.save(out);
				}
			}
		}
		clusteredReadsCache.clear();
		totalReads = 0;
	}
}
