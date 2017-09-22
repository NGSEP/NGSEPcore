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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.main.io.ConcatGZIPInputStream;


public class ReadsDemultiplex {
	
	private Logger log = Logger.getLogger(ReadsDemultiplex.class.getName());
	private Map<String,Map<String,String>> barcodeMap = new TreeMap<String, Map<String,String>>();
	
	private Map<String, String> samplesPerBarcode = new TreeMap<String, String>();
	//Sorted list of barcodes for log time non-exact search
	private String [] sortedBarcodes;
	private int [] nucleotideStarts = new int[4];
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

	public void loadIndex(String indexFile) throws IOException {
		FileInputStream fis = null;
		BufferedReader in = null;
		try {
			fis = new FileInputStream(indexFile);
			in = new BufferedReader(new InputStreamReader(fis));
			String line = in.readLine();
			//Skip header
			log.info("Loading index file: "+indexFile+". Skipping header");
			line = in.readLine();
			String flowcell= null;
			String lane = null;
			Map<String, String> laneSamplesPerBarcode = null;
			while (line != null) {
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
					laneSamplesPerBarcode.put(items[2]+"-"+items[3], items[4]);
				} else {
					log.info("Barcode: "+items[2]+" sample: "+items[3]);
					laneSamplesPerBarcode.put(items[2], items[3]);
				}
				
				line = in.readLine();
			}
			if(flowcell!=null) barcodeMap.put(makeLaneKey(flowcell, lane), laneSamplesPerBarcode);
		} finally {
			if(in!=null) in.close();
			if (fis!=null)fis.close();
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
			barcodesSet = barcodeKeys;
		}
		sortedBarcodes =  barcodesSet.toArray(new String[0]);
		Arrays.sort(sortedBarcodes);
		nucleotideStarts[0] = 0;
		int j=1;
		for (int i=0;i<sortedBarcodes.length-1;i++) {
			char c1 = sortedBarcodes[i].charAt(0);
			char c2 = sortedBarcodes[i+1].charAt(0);
			if(c1!=c2) {
				nucleotideStarts[j] = i+1;
				j++;
			}
		}
	}

	private String makeLaneKey(String flowcell, String lane) {
		return flowcell+"_"+lane;
	}

	public void demultiplex(String filename) throws IOException {
		log.info("Demultiplexing file: "+filename);
		FileInputStream fis = null;
		try {
			fis = new FileInputStream(filename);
			if(filename.toLowerCase().endsWith(".gz")) {
				demultiplex (new ConcatGZIPInputStream(fis));
			} else {
				demultiplex(fis);
			}
		} finally {
			if (fis!=null)fis.close();
		}
	}
	public void demultiplex(InputStream input) throws IOException {
		Map<String, PrintStream> outFiles = new TreeMap<String, PrintStream>();
		Map<String,Integer> counts = new TreeMap<String, Integer>();
		Pattern trimRegexp = null;
		if(trimSequence!=null) {
			trimRegexp = Pattern.compile(DegenerateSequence.makeRegularExpression(trimSequence));
		}
		BufferedReader in = new BufferedReader(new InputStreamReader(input));
		int total = 0;
		int notFound = 0;
		int tooShort = 0;
		RawRead read = RawRead.load(in);
		while (read != null) {
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
			if (progressNotifier!=null && total%10000==0) {
				int progress = (int) (total/10000);
				if (!progressNotifier.keepRunning(progress)) break;
			}
			read = RawRead.load(in);
		}
		for(PrintStream s:outFiles.values()) {
			s.flush();
			s.close();
		}
		log.info("Demultiplexing report");
		log.info("Total reads\t"+total);
		
		for(String sampleName:counts.keySet()) {
			log.info(""+sampleName+"\t"+counts.get(sampleName));
		}
		log.info("Not found\t"+notFound);
		log.info("Fragment too short\t"+tooShort);
	}

	private PrintStream getOutputStream(String filePrefix)
			throws FileNotFoundException, IOException {
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
		log.info("Demultiplexing paired end files "+filename1+" "+filename2);
		FileInputStream fis1 = null;
		FileInputStream fis2 = null;
		try {
			fis1 = new FileInputStream(filename1);
			InputStream is1= fis1;
			if(filename1.endsWith(".gz")) {
				is1 = new ConcatGZIPInputStream(fis1);
			}
			fis2 = new FileInputStream(filename2);
			InputStream is2= fis2;
			if(filename2.endsWith(".gz")) {
				is2 = new ConcatGZIPInputStream(fis2);
			}
			demultiplex(is1, is2);
		} finally {
			if (fis1!=null)fis1.close();
			if (fis2!=null)fis2.close();
		}
	}
	public void demultiplex(InputStream i1, InputStream i2) throws IOException {
		Map<String, PrintStream> outFiles1 = new TreeMap<String, PrintStream>();
		Map<String, PrintStream> outFiles2 = new TreeMap<String, PrintStream>();
		Map<String,Integer> counts = new TreeMap<String, Integer>();
		Pattern trimRegexp = null;
		if(trimSequence!=null) {
			trimRegexp = Pattern.compile(DegenerateSequence.makeRegularExpression(trimSequence));
		}
		BufferedReader in1 = new BufferedReader(new InputStreamReader(i1));
		BufferedReader in2 = new BufferedReader(new InputStreamReader(i2));
		int total = 0;
		int notFound = 0;
		int tooShort = 0;
		RawRead read1 = RawRead.load(in1);
		RawRead read2 = RawRead.load(in2);
		while (read1 != null && read2 != null) {
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
			if (progressNotifier!=null && total%10000==0) {
				int progress = (int) (total/10000);
				if (!progressNotifier.keepRunning(progress)) break;
			}
			read1 = RawRead.load(in1);
			read2 = RawRead.load(in2);
		}
		for(PrintStream s:outFiles1.values()) {
			s.flush();
			s.close();
		}
		for(PrintStream s:outFiles2.values()) {
			s.flush();
			s.close();
		}
		log.info("Demultiplex report");
		log.info("Total fragments\t"+total);
		
		for(String sampleName:counts.keySet()) {
			log.info(""+sampleName+"\t"+counts.get(sampleName));
		}
		log.info("Not found\t"+notFound);
		log.info("Fragment too short\t"+tooShort);
	}

	private String findBarcode(String readSeq) {
		char firstC = readSeq.charAt(0);
		int start=0;
		int end = 0;
		if(firstC == 'A') {
			start = nucleotideStarts[0];
			end = nucleotideStarts[1];
		}
		else if(firstC == 'C') {
			start = nucleotideStarts[1];
			end = nucleotideStarts[2];
		}
		else if(firstC == 'G') {
			start = nucleotideStarts[2];
			end = nucleotideStarts[3];
		}
		else if(firstC == 'T') {
			start = nucleotideStarts[3];
			end = sortedBarcodes.length;
		}
		else return null;
		String selectedBarcode = null;
		for(int i=start;i<end;i++) {
			String barcode = sortedBarcodes[i];
			boolean barcodeFound = false;
			if(readSeq.startsWith(barcode)) barcodeFound = true;
			/*Pattern p = regExps.get(barcode);
			Matcher m = p.matcher(readSeq);
			if(m.find()) {
				if(m.start()==0) barcodeFound = true;
			}
			*/
			if(barcodeFound && (selectedBarcode==null || selectedBarcode.length()<barcode.length()) ) {
				selectedBarcode = barcode;
			}
		}
		return selectedBarcode;
	}
	
}
