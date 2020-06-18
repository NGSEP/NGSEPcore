package ngsep.simulation;

import java.io.ByteArrayOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.RawRead;

public class SingleReadsSimulator {
	
	// Constants for default values
	public static final int DEF_NUM_READS = 30000;
	public static final int DEF_MEAN_READ_LENGTH = 20000;
	public static final int DEF_STDEV_READ_LENGTH = 10000;
	public static final int DEF_MIN_READ_LENGTH = 50;
	public static final double DEF_SUBSTITUTION_ERROR_RATE = 0.01;
	public static final double DEF_INDEL_ERROR_RATE = 0.01;
	public static final byte OUT_FORMAT_FASTQ = 0;
	public static final byte OUT_FORMAT_FASTA = 1;
	
	// Logging and progress
	private Logger log = Logger.getLogger(SingleReadsSimulator.class.getName());
	private ProgressNotifier progressNotifier = null;

	// Parameters
	private ReferenceGenome genome;
	private String outputFile = null;
	private int numberOfReads = DEF_NUM_READS;
	private int meanReadLength = DEF_MEAN_READ_LENGTH;
	private int stdevReadlength = DEF_STDEV_READ_LENGTH;
	private int minReadLength = DEF_MIN_READ_LENGTH;
	private double substitutionErrorRate = DEF_SUBSTITUTION_ERROR_RATE;
	private double indelErrorRate = DEF_INDEL_ERROR_RATE;
	private byte outFormat = OUT_FORMAT_FASTQ;

	private final static Random rnd = new Random();

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

	public int getNumberOfReads() {
		return numberOfReads;
	}
	public void setNumberOfReads(int numberOfReads) {
		this.numberOfReads = numberOfReads;
	}
	public void setNumberOfReads(String value) {
		this.setNumberOfReads((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	public int getMeanReadLength() {
		return meanReadLength;
	}
	public void setMeanReadLength(int meanReadLength) {
		this.meanReadLength = meanReadLength;
	}
	public void setMeanReadLength(String value) {
		this.setMeanReadLength((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	public int getStdevReadlength() {
		return stdevReadlength;
	}
	public void setStdevReadlength(int stdevReadlength) {
		this.stdevReadlength = stdevReadlength;
	}
	public void setStdevReadlength(String value) {
		this.setStdevReadlength((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	public double getSubstitutionErrorRate() {
		return substitutionErrorRate;
	}
	public void setSubstitutionErrorRate(double substitutionErrorRate) {
		this.substitutionErrorRate = substitutionErrorRate;
	}
	public void setSubstitutionErrorRate(String value) {
		this.setSubstitutionErrorRate((double) OptionValuesDecoder.decode(value, Double.class));
	}

	public double getIndelErrorRate() {
		return indelErrorRate;
	}
	public void setIndelErrorRate(double indelErrorRate) {
		this.indelErrorRate = indelErrorRate;
	}
	public void setIndelErrorRate(String value) {
		this.setIndelErrorRate((double) OptionValuesDecoder.decode(value, Double.class));
	}
	
	public byte getOutFormat() {
		return outFormat;
	}
	public void setOutFormat(byte outFormat) {
		this.outFormat = outFormat;
	}
	public void setOutFormat(String value) {
		this.setOutFormat((byte) OptionValuesDecoder.decode(value, Byte.class));
	}

	public static void main(String[] args) throws Exception {
		SingleReadsSimulator instance = new SingleReadsSimulator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}
	
	public void run() throws IOException {
		logParameters();
		if(genome==null) throw new IOException("A file with the reference genome is a required parameter");
		if(outputFile==null) throw new IOException("The output file is a required parameter");
		simulate(outputFile);
	}

	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		if (genome!=null) out.println("Genome for simulation loaded from file: "+genome.getFilename());
		out.println("Output file path: " + outputFile);
		out.println("Reference total length:" + genome.getTotalLength());
		out.println("Reads:" + numberOfReads + "   ~N(mean: " + meanReadLength + ", sdev: " + stdevReadlength + ")");
		out.println("Substitution error rate: " + substitutionErrorRate);
		out.println("Indel error rate: " + indelErrorRate);
	}
	
	public void simulate(String outPath) throws IOException {
		long totalLength = genome.getTotalLength();
		int nSeqs = genome.getNumSequences();
		long[] cumulativeStarts = new long[nSeqs];
		cumulativeStarts[0] = 0;
		for (int i = 1; i < nSeqs; i++) {
			cumulativeStarts[i] = cumulativeStarts[i - 1] + genome.getSequenceByIndex(i - 1).getLength();
		}
		try (OutputStream os = new GZIPOutputStream(new FileOutputStream(outPath));
			 PrintStream out = new PrintStream(os)) {
			for (int i = 0; i < numberOfReads; i++) {
				int readLength;
				long nextStart;
				QualifiedSequence seq = null;
				int relStart = 0;
				byte reverse = 0;
				String read = null;
				for (int j = 0; j < 100; j++) {
					readLength = (int) (rnd.nextGaussian() * stdevReadlength + meanReadLength);
					if(readLength<minReadLength) continue;
					long nextLong = rnd.nextLong();
					nextStart = nextLong % (totalLength - readLength);
					if (nextStart < 0) continue;

					int idx1 = Arrays.binarySearch(cumulativeStarts, nextStart);

					int sequenceIdx;
					if (idx1 >= 0)
						sequenceIdx = idx1;
					else {
						sequenceIdx = -idx1 - 2;
					}
					seq = genome.getSequenceByIndex(sequenceIdx);
					relStart = (int) (nextStart - cumulativeStarts[sequenceIdx]);
					int relEnd = relStart + readLength;
					if (relEnd <= seq.getLength()) {
						read = seq.getCharacters().subSequence(relStart, relEnd).toString();
						break;
					}
				}
				if (read == null) continue;
				if (rnd.nextBoolean()) {
					reverse = 1;
					read = DNAMaskedSequence.getReverseComplement(read).toString();
				}
				String finalRead = generateErrors(read);
				String readId = seq.getName() + "_" + (relStart+1) + "_" + reverse;
				if(outFormat == OUT_FORMAT_FASTA) {
					out.println(">" + readId);
					out.println(finalRead);
				} else {
					out.println("@" + readId);
					out.println(finalRead);
					out.println("+");
					out.println(RawRead.generateFixedQSString('5', finalRead.length()));
				}
				
			}
		}
	}

	private String generateErrors(String read) {
		String alphabet = DNASequence.BASES_STRING;
		int len = read.length();
		StringBuilder answer = new StringBuilder(len);

		for(int i=0;i<len;i++) {
			
			if(rnd.nextDouble()<indelErrorRate) {
				//Generate random indels
				int length = 0;
				for (int j=0; j<100 && length == 0;j++) {
					length = (int) Math.round(rnd.nextGaussian()*2.0);
				}
				length = Math.min(length, 10);
				length = Math.max(length, -10);
				if(length < 0) {
					//Deletion
					i-=(length+1);
				} else {
					//Insertion
					answer.append(read.charAt(i));
					for (int j=0; j<length;j++) answer.append(alphabet.charAt(0));
				}
			} else {
				char c = read.charAt(i);
				if(rnd.nextDouble()<substitutionErrorRate) {
					// Generate random substitution
					char c2 = c;
					for (int j=0; j<100 && c == c2;j++) {
						c2 = alphabet.charAt(rnd.nextInt(alphabet.length()));
					}
					c = c2;
				}
				answer.append(c);
			}
		}
		return answer.toString();
	}

}
