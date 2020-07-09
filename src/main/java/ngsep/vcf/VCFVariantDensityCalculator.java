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
package ngsep.vcf;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

/**
 * Script to calculate density of variants in a VCF file across the genome
 * @author Jorge Duitama
 *
 */
public class VCFVariantDensityCalculator {

	// Constants for default values
	public static final int DEF_WINDOW_LENGTH=100000;
	
	// Logging and progress
	private Logger log = Logger.getLogger(VCFVariantDensityCalculator.class.getName());
	private ProgressNotifier progressNotifier=null;
	
	// Parameters
	private String inputFile = null;
	private ReferenceGenome genome;
	private String outputFile = null;
	private int windowLength = DEF_WINDOW_LENGTH;
	
	// Get and set methods
	
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}

	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}

	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	public void setGenome(String genomeFile) throws IOException {
		setGenome(OptionValuesDecoder.loadGenome(genomeFile,log));
	}
	
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}

	public int getWindowLength() {
		return windowLength;
	}
	public void setWindowLength(int windowLength) {
		this.windowLength = windowLength;
	}
	public void setWindowLength(String value) {
		setWindowLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public static void main(String[] args) throws Exception {
		VCFVariantDensityCalculator instance = new VCFVariantDensityCalculator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();

	}
	
	public void run() throws IOException {
		if (genome == null) throw new IOException("The file with the reference genome is a required parameter");
		log.info("Loaded reference genome from: "+genome.getFilename());
		log.info("Window length: "+getWindowLength());
		if(inputFile==null) {
			log.info("Reading from standard input");
			if(outputFile == null) run(System.in, System.out);
			else {
				try (PrintStream out = new PrintStream(outputFile)) {
					run(System.in, out);
				}
			}
		} else {
			log.info("Reading from file: "+inputFile);
			if(outputFile == null) run(inputFile,System.out);
			else {
				try (PrintStream out = new PrintStream(outputFile)) {
					run(inputFile, out);
				}
			}
		}
		log.info("Process finished");
	}

	public void run(String filename, PrintStream out) throws IOException {
		try (VCFFileReader in = new VCFFileReader(filename)) { 
			run(in, out);
		}
	}
	public void run(InputStream fis, PrintStream out) throws IOException {
		try (VCFFileReader in = new VCFFileReader(fis)) {
			run(in, out);
		}	
	}
	public void run(VCFFileReader in, PrintStream out) throws IOException {
		if(log!=null)in.setLog(log);
		in.setLoadMode(VCFFileReader.LOAD_MODE_MINIMAL);
		QualifiedSequenceList sequences = genome.getSequencesMetadata();
		in.setSequences(sequences);
		Iterator<VCFRecord> it = in.iterator();
		int n=0;
		
		String seqName = null;
		int endWindow = 0;
		int count = 0;
		int idxSeqNames = 0;
		while(it.hasNext()) {
			VCFRecord record = it.next();
			if(!record.getSequenceName().equals(seqName)) {
				if(seqName!=null) out.println(seqName+"\t"+(endWindow-windowLength+1)+"\t"+endWindow+"\t"+count);
				QualifiedSequence seq = sequences.get(idxSeqNames); 
				while (!seq.getName().equals(record.getSequenceName())) {
					processZeroes (seq,endWindow, out);
					idxSeqNames++;
					endWindow = 0;
					seq = sequences.get(idxSeqNames);
				}
				seqName = record.getSequenceName();
				endWindow = windowLength;
				count = 0;
			} else if (record.getFirst()>endWindow) {
				if(seqName!=null) out.println(seqName+"\t"+(endWindow-windowLength+1)+"\t"+endWindow+"\t"+count);
				endWindow+=windowLength;
				count = 0;
			}
			for(;endWindow<record.getFirst();endWindow+=windowLength) {
				out.println(seqName+"\t"+(endWindow-windowLength+1)+"\t"+endWindow+"\t0");
			}
			count++;
			n++;
			if (progressNotifier!=null && n%1000==0) {
				int progress = n/1000;
				if (!progressNotifier.keepRunning(progress)) {
					out.flush();
					return;
				}
			}
		}
		if(seqName!=null) out.println(seqName+"\t"+(endWindow-windowLength+1)+"\t"+endWindow+"\t"+count);
		while (idxSeqNames<sequences.size()) {
			QualifiedSequence seq = sequences.get(idxSeqNames);
			processZeroes (seq,endWindow, out);
			idxSeqNames++;
			endWindow = 0;
		}
		
	}

	private void processZeroes(QualifiedSequence seq, int endWindow, PrintStream out) {
		while(endWindow<seq.getLength()) {
			endWindow+=windowLength;
			out.println(seq.getName()+"\t"+(endWindow-windowLength+1)+"\t"+endWindow+"\t0");
			
		}
		
	}

}
