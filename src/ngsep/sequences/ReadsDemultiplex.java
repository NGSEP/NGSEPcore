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
package ngsep.sequences;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.io.FastqFileReader;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class ReadsDemultiplex {
	
	private Logger log = Logger.getLogger(ReadsDemultiplex.class.getName());
	private Map<String,Map<String,String>> barcodeMap = new TreeMap<String, Map<String,String>>();
	
	private Map<String, String> samplesPerBarcode = new TreeMap<String, String>();
	//Sorted list of barcodes for log time non-exact search
	private String [] sortedBarcodes;
	private int [] [] ntPairFirsts = new int[4][4];
	private int [] [] ntPairLasts = new int[4][4];
	//private Map<String,Pattern> regExps = new TreeMap<String, Pattern>();
	private String outDirectory = ".";
	private String prefix = "";
	private String trimSequence = null;
	private int minReadLength = 40;
	private boolean uncompressedOutput = false;
	private boolean dualBarcode = false;
	private String laneFilesDescriptor = null;
	private String flowcell = null;
	private String lane = null;
	
	private ProgressNotifier progressNotifier = null;
	
	private Pattern trimRegexp = null;
	
	//Statistics
	private int total = 0;
	private int notFound = 0;
	private int tooShort = 0;
	private Map<String,Integer> counts;
	

	public static void main(String[] args) throws Exception {
		ReadsDemultiplex instance = new ReadsDemultiplex();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		
		String indexFile = args[i++];
		instance.loadIndex(indexFile);
		if(instance.laneFilesDescriptor !=null) {
			instance.demultiplexGroup();
			return;
		}
		if(instance.flowcell== null || instance.lane==null) {
			System.err.println("Either a lane files descriptor or a flow cell and a lane should be provided");
			System.exit(1);
		}
		instance.loadLaneInfo(instance.flowcell,instance.lane);
		if(i==args.length || "-".equals(args[i]) ) {
			instance.demultiplex(System.in);
		} else if(i+1==args.length){
			instance.demultiplex(args[i++]);
		} else {
			String f1 = args[i++];
			String f2 = args[i++];
			instance.demultiplex(f1,f2);
		}
	}

	public void demultiplexGroup() throws IOException {
		if(laneFilesDescriptor==null) return;
		FileInputStream fis = null;
		BufferedReader in = null;
		try {
			fis = new FileInputStream(laneFilesDescriptor);
			in = new BufferedReader(new InputStreamReader(fis));
			String line = in.readLine();
			while (line != null) {
				String[] items = line.split("\t");
				loadLaneInfo(items[0], items[1]);
				if(items.length<4) demultiplex(items[2]);
				else demultiplex(items[2],items[3]);
				line = in.readLine();
			}
		} finally {
			if (in != null)
				in.close();
			if (fis != null)
				fis.close();
		}
		
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

	
	public String getOutDirectory() {
		return outDirectory;
	}

	public void setOutDirectory(String outDirectory) {
		this.outDirectory = outDirectory;
	}
	
	public String getLaneFilesDescriptor() {
		return laneFilesDescriptor;
	}

	public void setLaneFilesDescriptor(String laneFilesDescriptor) {
		this.laneFilesDescriptor = laneFilesDescriptor;
	}
	

	public String getFlowcell() {
		return flowcell;
	}

	public void setFlowcell(String flowcell) {
		this.flowcell = flowcell;
	}

	public String getLane() {
		return lane;
	}

	public void setLane(String lane) {
		this.lane = lane;
	}

	public String getPrefix() {
		return prefix;
	}

	public void setPrefix(String prefix) {
		this.prefix = prefix;
	}

	public int getMinReadLength() {
		return minReadLength;
	}

	public void setMinReadLength(int minReadLength) {
		this.minReadLength = minReadLength;
	}
	
	public void setMinReadLength(Integer minReadLength) {
		this.setMinReadLength(minReadLength.intValue());
	}

	public boolean isDualBarcode() {
		return dualBarcode;
	}

	public void setDualBarcode(boolean dualBarcode) {
		this.dualBarcode = dualBarcode;
	}

	public void setDualBarcode(Boolean dualBarcode) {
		this.setDualBarcode(dualBarcode.booleanValue());
	}

	public String getTrimSequence() {
		return trimSequence;
	}

	public void setTrimSequence(String trimSequence) {
		this.trimSequence = trimSequence;
		if(trimSequence!=null) {
			trimRegexp = Pattern.compile(DegenerateSequence.makeRegularExpression(trimSequence));
		} else {
			trimRegexp = null;
		}
 	}
	
	

	public boolean isUncompressedOutput() {
		return uncompressedOutput;
	}

	public void setUncompressedOutput(boolean uncompressedOutput) {
		this.uncompressedOutput = uncompressedOutput;
	}
	
	public void setUncompressedOutput(Boolean uncompressedOutput) {
		this.uncompressedOutput = uncompressedOutput;
	}

	/**
	 * Loads barcodes from the given file
	 * @param indexFile Text file with the id of the sample corresponding to each single barcode, or barcode pair for dual barcoding 
	 * @throws IOException If the file can not be read
	 */
	public void loadIndex(String indexFile) throws IOException {
		try (FileInputStream fis = new FileInputStream(indexFile);
			 BufferedReader in = new BufferedReader(new InputStreamReader(fis));) {
			String line = in.readLine();
			//Skip header
			log.info("Loading index file: "+indexFile);
			line = in.readLine();
			String flowcell= null;
			String lane = null;
			Map<String, String> laneSamplesPerBarcode = null;
			for (int i=1;line != null;i++) {
				String[] items = line.split("\t");
				
				if(!items[0].equals(flowcell) || !items[1].equals(lane)) {
					if(flowcell!=null) barcodeMap.put(makeLaneKey(flowcell, lane), laneSamplesPerBarcode);
					flowcell = items[0];
					lane = items[1];
					laneSamplesPerBarcode = new TreeMap<String, String>();
					log.info("Loading barcode map for flowcell: "+flowcell+" lane: "+lane);
				}
				if (dualBarcode) {
					log.info("Barcode 1: "+items[2]+" barcode 2: "+items[3]+" sample: "+items[4]);
					if(!DNASequence.isDNA(items[2])) throw new IOException("Barcode: "+items[2]+" at line "+i+" is not a DNA sequence");
					if(!DNASequence.isDNA(items[3])) throw new IOException("Barcode: "+items[3]+" at line "+i+" is not a DNA sequence");
					laneSamplesPerBarcode.put(items[2]+"-"+items[3], items[4]);
				} else {
					log.info("Barcode: "+items[2]+" sample: "+items[3]);
					if(!DNASequence.isDNA(items[2])) throw new IOException("Barcode: "+items[2]+" at line "+i+" is not a DNA sequence");
					laneSamplesPerBarcode.put(items[2], items[3]);
				}
				
				line = in.readLine();
			}
			if(flowcell!=null) barcodeMap.put(makeLaneKey(flowcell, lane), laneSamplesPerBarcode);
		}
	}
	public void loadLaneInfo(String flowcell, String lane) throws IOException {
		System.out.println("Using barcode map for flowcell: "+flowcell+" lane: "+lane);
		String key = makeLaneKey(flowcell,lane);
		samplesPerBarcode = barcodeMap.get(key);
		if(samplesPerBarcode == null) throw new IOException("Flowcell: "+flowcell+" and lane: "+lane+ " not found in barcode map");
		Set<String> barcodeKeys = samplesPerBarcode.keySet();
		
		Set<String> barcodesSet;
		if(dualBarcode) {
			barcodesSet = new TreeSet<>();
			for(String barcodeKey:barcodeKeys) {
				int i = barcodeKey.indexOf('-');
				barcodesSet.add(barcodeKey.substring(0, i));
				barcodesSet.add(barcodeKey.substring(i+1));
			}	
		} else {
			barcodesSet = new TreeSet<>(barcodeKeys);
		}
		sortedBarcodes = barcodesSet.toArray(new String[0]);
		for(int i=0;i<ntPairFirsts.length;i++) {
			Arrays.fill(ntPairFirsts[i], -1);
			Arrays.fill(ntPairLasts[i], -1);
		}
		String alphabet = DNASequence.BASES_STRING;
		for (int k=0;k<sortedBarcodes.length;k++) {
			String barcode = sortedBarcodes[k];
			int i = alphabet.indexOf(barcode.charAt(0));
			int j = alphabet.indexOf(barcode.charAt(1));
			if(ntPairFirsts[i][j]==-1)ntPairFirsts[i][j] = k;
			ntPairLasts[i][j] = k;
		}
	}

	private String makeLaneKey(String flowcell, String lane) {
		return flowcell+"_"+lane;
	}

	/**
	 * Demultiplex a single fastq file
	 * @param filename name of the fastq file to demultiplex
	 * @throws IOException If the file can not be read
	 */
	public void demultiplex(String filename) throws IOException {
		log.info("Demultiplexing file: "+filename);
		
		Map<String, PrintStream> outFiles = new TreeMap<String, PrintStream>();
		initCounts();
		
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				processRead(read, outFiles);
				if(!keepRunning()) break;
			}
		} finally {
			for(PrintStream s:outFiles.values()) {
				s.flush();
				s.close();
			}
		}
		printStatistics();
	}
	public void demultiplex(InputStream input) throws IOException {
		Map<String, PrintStream> outFiles = new TreeMap<String, PrintStream>();
		initCounts();
		
		try (FastqFileReader reader = new FastqFileReader(input)) {
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				processRead(read, outFiles);
				if(!keepRunning()) break;
			}
		} finally {
			for(PrintStream s:outFiles.values()) {
				s.flush();
				s.close();
			}
		}
		printStatistics();
	}

	private void initCounts() {
		counts = new TreeMap<String, Integer>();
		total = 0;
		notFound = 0;
		tooShort = 0;
	}
	
	private void printStatistics() {
		log.info("Demultiplexing report");
		log.info("Total fragments\t"+total);
		
		for(String sampleName:counts.keySet()) {
			log.info(""+sampleName+"\t"+counts.get(sampleName));
		}
		log.info("Not found\t"+notFound);
		log.info("Fragment too short\t"+tooShort);
	}

	/**
	 * Demultiplex one single read
	 * @param read to demultiplex
	 * @param outFiles Map of files to choose
	 * @return boolean true if the process should keep running
	 * @throws IOException 
	 */
	private void processRead(RawRead read, Map<String, PrintStream> outFiles) throws IOException {
		total++;
		String barcode = findBarcode(read.getSequenceString());
		if (barcode == null) {
			notFound++;
		} else {
			String sampleId = samplesPerBarcode.get(barcode);
			PrintStream out = outFiles.get(sampleId);
			if(out==null) {
				out = getOutputStream(sampleId);
				outFiles.put(sampleId, out);
			}
			int barcodeLength = barcode.length();
			//Trim barcode
			read.trimFirstNucleotides(barcodeLength);
			//Trim end if sequence appears
			read.trimFromSequence(trimRegexp);
			if(read.getLength()>=minReadLength) {
				read.save(out);
				Integer count = counts.get(sampleId);
				if(count==null) counts.put(sampleId, 1);
				else counts.put(sampleId, count+1);
			} else {
				tooShort++;
			}
		}
	}

	/**
	 * Reports progress and asks the progress notifier if the process should keep running
	 * @return boolean true if the process should keep running
	 */
	private boolean keepRunning() {
		if (progressNotifier!=null && total%10000==0) {
			int progress = (int) (total/10000);
			return progressNotifier.keepRunning(progress);
		}
		return true;
	}

	private PrintStream getOutputStream(String filePrefix) throws IOException {
		PrintStream out;
		String filename =  filePrefix + ".fastq";
		if(!uncompressedOutput) filename+=".gz";
		if(prefix!=null && Character.isDigit(filename.charAt(0))) {
			filename = prefix+filename;
		}
		if(outDirectory!=null) {
			filename = outDirectory + File.separator + filename;
		}
		File f = new File(filename);
		OutputStream os = new FileOutputStream(filename,f.exists());
		if(!uncompressedOutput) os = new GZIPOutputStream(os); 
		out = new PrintStream(os);
		return out;
	}
	public void demultiplex(String filename1, String filename2) throws IOException {
		Map<String, PrintStream> outFiles1 = new TreeMap<String, PrintStream>();
		Map<String, PrintStream> outFiles2 = new TreeMap<String, PrintStream>();
		initCounts();
		try (FastqFileReader in1 = new FastqFileReader(filename1);
			 FastqFileReader in2 = new FastqFileReader(filename2)) {
			Iterator<RawRead> it1 = in1.iterator();
			Iterator<RawRead> it2 = in2.iterator();
			while (it1.hasNext() && it2.hasNext()) {
				RawRead read1 = it1.next();
				RawRead read2 = it2.next();
				demultiplexReadPair(read1, read2, outFiles1, outFiles2);
				if(!keepRunning()) break;
			}
		} finally {
			for(PrintStream s:outFiles1.values()) {
				s.flush();
				s.close();
			}
			for(PrintStream s:outFiles2.values()) {
				s.flush();
				s.close();
			}
		}
		printStatistics();
	}
	public void demultiplex(InputStream i1, InputStream i2) throws IOException {
		Map<String, PrintStream> outFiles1 = new TreeMap<String, PrintStream>();
		Map<String, PrintStream> outFiles2 = new TreeMap<String, PrintStream>();
		initCounts();
		try (FastqFileReader in1 = new FastqFileReader(i1);
			 FastqFileReader in2 = new FastqFileReader(i2)) {
			Iterator<RawRead> it1 = in1.iterator();
			Iterator<RawRead> it2 = in2.iterator();
			while (it1.hasNext() && it2.hasNext()) {
				RawRead read1 = it1.next();
				RawRead read2 = it2.next();
				demultiplexReadPair(read1, read2, outFiles1, outFiles2);
				if(!keepRunning()) break;
			}
		} finally {
			for(PrintStream s:outFiles1.values()) {
				s.flush();
				s.close();
			}
			for(PrintStream s:outFiles2.values()) {
				s.flush();
				s.close();
			}
		}
		printStatistics();
	}

	/**
	 * Demultiplex two reads as a pair
	 * @param read1 First read of the pair
	 * @param read2 Second read of the pair
	 * @param outFiles1 List of out files to choose
	 * @param outFiles2 List of out files to choose
	 * @throws IOException If the out files can not be written
	 */
	private void demultiplexReadPair(RawRead read1, RawRead read2, Map<String, PrintStream> outFiles1, Map<String, PrintStream> outFiles2) throws IOException {
		total++;
		String barcode1;
		String barcode2=null;
		String barcodeM=null;
		if(dualBarcode) {
			barcode1 = findBarcode(read1.getSequenceString());
			barcode2 = findBarcode(read2.getSequenceString());
			if(barcode1!=null && barcode2!=null) {
				barcodeM = barcode1+"-"+barcode2;
				if(!samplesPerBarcode.containsKey(barcodeM)) barcodeM = null;
			}
		} else {
			barcode1 = findBarcode(read1.getSequenceString());
			barcodeM = barcode1;
		}
		
		if (barcodeM == null) {
			notFound++;
		} else {
			String sampleId = samplesPerBarcode.get(barcodeM);
			if (sampleId == null) System.err.println("Sample not found with barcode key: "+barcodeM);
			PrintStream out1 = outFiles1.get(sampleId);
			PrintStream out2 = outFiles2.get(sampleId);
			if(out1==null) {
				out1 = getOutputStream(sampleId+"_1");
				out2 = getOutputStream(sampleId+"_2");
				outFiles1.put(sampleId, out1);
				outFiles2.put(sampleId, out2);
			}
			
			//Trim barcode
			read1.trimFirstNucleotides(barcode1.length());
			if(dualBarcode) {
				read2.trimFirstNucleotides(barcode2.length());
			}
			//Trim end if sequence appears
			int l1 = read1.getLength();
			read1.trimFromSequence(trimRegexp);
			if(read1.getLength()!=l1) {
				read2.trimToLength(read1.getLength());
			}
			if(read1.getLength()>=minReadLength) {
				read1.save(out1);
				read2.save(out2);
				Integer count = counts.get(sampleId);
				if(count==null) counts.put(sampleId, 1);
				else counts.put(sampleId, count+1);
			} else {
				tooShort++;
			}
		}
	}

	/**
	 * Finds the barcode matching the given sequence
	 * @param readSeq Sequence to barcode
	 * @return String barcode matching the given sequence. Null if the barcode was not found
	 */
	private String findBarcode(String readSeq) {
		String alphabet = DNASequence.BASES_STRING;
		int firstIndex = alphabet.indexOf(readSeq.charAt(0));
		int secondIndex = alphabet.indexOf(readSeq.charAt(1));
		if(firstIndex <0 || secondIndex < 0 ) return null;
		int first = ntPairFirsts[firstIndex][secondIndex];
		int last = ntPairLasts[firstIndex][secondIndex];
		if(first==-1) return null;
		String selectedBarcode = null;
		for(int i=first;i<=last;i++) {
			String barcode = sortedBarcodes[i];
			boolean barcodeFound = false;
			if(readSeq.startsWith(barcode)) barcodeFound = true;
			//The last condition is to taking into account cases where one barcode is substring of another
			//and pick the longest matching barcode
			if(barcodeFound && (selectedBarcode==null || selectedBarcode.length()<barcode.length()) ) {
				selectedBarcode = barcode;
			}
		}
		return selectedBarcode;
	}
	
}
