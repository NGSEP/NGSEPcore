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
import java.util.TreeMap;
import java.util.logging.Logger;
import java.util.regex.Matcher;
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
	private ProgressNotifier progressNotifier = null;

	public static void main(String[] args) throws Exception {
		if (args.length == 0 || args[0].equals("-h") || args[0].equals("--help")){
			CommandsDescriptor.getInstance().printHelp(ReadsDemultiplex.class);
			return;
		}
		ReadsDemultiplex instance = new ReadsDemultiplex();
		int i=0;
		String flowcell = null;
		String lane = null;
		String laneFilesDescriptor = null;
		while(i<args.length && args[i].charAt(0)=='-') {
			if("-t".equals(args[i])) {
				i++;
				instance.trimSequence = args[i];
			} else if ("-p".equals(args[i])) {
				i++;
				instance.prefix = args[i];
			}  else if ("-o".equals(args[i])) {
				i++;
				instance.outDirectory = args[i];
			} else if ("-u".equals(args[i])) {
				instance.uncompressedOutput = true;
			} else if ("-f".equals(args[i])) {
				i++;
				flowcell = args[i]; 
			} else if ("-l".equals(args[i])) {
				i++;
				lane = args[i]; 
			} else if ("-d".equals(args[i])) {
				i++;
				laneFilesDescriptor = args[i]; 
			} else {
				System.err.println("Unrecognized option :"+args[i]);
				CommandsDescriptor.getInstance().printHelp(ReadsDemultiplex.class);
				System.exit(1);
			}
			i++;
		}
		String indexFile = args[i++];
		instance.loadIndex(indexFile);
		if(laneFilesDescriptor !=null) {
			instance.demultiplexGroup(laneFilesDescriptor);
			return;
		}
		if(flowcell== null || lane==null) {
			System.err.println("Either a lane files descriptor or a flow cell and a lane should be provided");
			System.exit(1);
		}
		instance.loadLaneInfo(flowcell, lane);
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

	public void demultiplexGroup(String laneFilesDescriptor) throws IOException {
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
				log.info("Barcode: "+items[2]+" sample: "+items[3]);
				laneSamplesPerBarcode.put(items[2], items[3]);
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
		sortedBarcodes = samplesPerBarcode.keySet().toArray(new String[0]);
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
			if(filename.endsWith(".gz")) {
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
		log.info("Deconvolution report");
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
		log.info("Deconvoluting paired end files "+filename1+" "+filename2);
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
			String barcode = findBarcode(read1.getSequenceString());
			if (barcode == null) {
				notFound++;
			} else {
				String sampleId = samplesPerBarcode.get(barcode);
				PrintStream out1 = outFiles1.get(sampleId);
				PrintStream out2 = outFiles2.get(sampleId);
				if(out1==null) {
					out1 = getOutputStream(sampleId+"_1");
					out2 = getOutputStream(sampleId+"_2");
					outFiles1.put(sampleId, out1);
					outFiles2.put(sampleId, out2);
				}
				int barcodeLength = barcode.length();
				//Trim barcode
				read1.trimFirstNucleotides(barcodeLength);
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
class RawRead extends QualifiedSequence {
	public RawRead(String id, CharSequence sequence, String qualityScores) {
		super(id,sequence,qualityScores);
	}
	public static RawRead load(BufferedReader in) throws IOException {
		String id = in.readLine();
		if(id==null) return null;
		String seq = in.readLine();
		if(seq==null) return null;
		String plus = in.readLine();
		if(plus==null) return null;
		String qs = in.readLine();
		if(qs==null) return null;
		return new RawRead(id, seq, qs);
	}
	public String getSequenceString() {
		return getCharacters().toString();
	}
	public void trimFirstNucleotides (int n) {
		setCharacters(getSequenceString().substring(n));
		setQualityScores(getQualityScores().substring(n));
	}
	public void trimFromSequence(Pattern regexp) {
		String sequence = getSequenceString();
		if(regexp==null) return;
		Matcher m = regexp.matcher(sequence);
		if(m.find()) {
			int start = m.start();
			trimToLength(start);
		}
	}
	public void trimToLength(int newLength) {
		String sequence = getSequenceString();
		setCharacters(sequence.substring(0,newLength));
		setQualityScores(getQualityScores().substring(0,newLength));
	}
	public void save (PrintStream out) {
		out.println(this.getName());
		out.println(this.getCharacters());
		out.println("+");
		out.println(this.getQualityScores());
	}
}
