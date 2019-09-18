package ngsep.simulation;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Logger;

import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;

public class SingleReadsSimulator {
	private Logger log = Logger.getLogger(SingleReadsSimulator.class.getName());
	private ProgressNotifier progressNotifier = null;

	public final static int DEF_NUM_READS = 30000;
	public final static int DEF_MEAN_READ_LENGTH = 10000;
	public final static int DEF_STDEV_READ_LENGTH = 2000;
	public final static double DEF_SUBSTITUTION_ERROR_RATE = 0.02;
	public final static double DEF_INDEL_ERROR_RATE = 0.01;
	public final static byte OUT_FORMAT_FASTQ = 0;
	public final static byte OUT_FORMAT_FASTA = 1;

	private int numberOfReads = DEF_NUM_READS;
	private int meanReadLength = DEF_MEAN_READ_LENGTH;
	private int stdevReadlength = DEF_STDEV_READ_LENGTH;
	private double substitutionErrorRate = DEF_SUBSTITUTION_ERROR_RATE;
	private double indelErrorRate = DEF_INDEL_ERROR_RATE;
	private byte outFormat = OUT_FORMAT_FASTQ;

	private final static Random rnd = new Random();
	private ReferenceGenome genome;

	/**
	 * @return the progressNotifier
	 */
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}

	/**
	 * @param progressNotifier the progressNotifier to set
	 */
	public void setProgressNotifier(ProgressNotifier progressNotifier) {
		this.progressNotifier = progressNotifier;
	}

	/**
	 * @return the numberOfReads
	 */
	public int getNumberOfReads() {
		return numberOfReads;
	}

	/**
	 * @param numberOfReads the numberOfReads to set
	 */
	public void setNumberOfReads(int numberOfReads) {
		this.numberOfReads = numberOfReads;
	}

	public void setNumberOfReads(String value) {
		this.setNumberOfReads((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	/**
	 * @return the meanReadLength
	 */
	public int getMeanReadLength() {
		return meanReadLength;
	}

	/**
	 * @param meanReadLength the meanReadLength to set
	 */
	public void setMeanReadLength(int meanReadLength) {
		this.meanReadLength = meanReadLength;
	}

	public void setMeanReadLength(String value) {
		this.setMeanReadLength((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	/**
	 * @return the stdevReadlength
	 */
	public int getStdevReadlength() {
		return stdevReadlength;
	}

	/**
	 * @return the log
	 */
	public Logger getLog() {
		return log;
	}

	/**
	 * @param log the log to set
	 */
	public void setLog(Logger log) {
		this.log = log;
	}

	/**
	 * @param stdevReadlength the stdevReadlength to set
	 */
	public void setStdevReadlength(int stdevReadlength) {
		this.stdevReadlength = stdevReadlength;
	}

	public void setStdevReadlength(String value) {
		this.setStdevReadlength((int) OptionValuesDecoder.decode(value, Integer.class));
	}

	/**
	 * @return the substitutionErrorRate
	 */
	public double getSubstitutionErrorRate() {
		return substitutionErrorRate;
	}

	/**
	 * @param substitutionErrorRate the substitutionErrorRate to set
	 */
	public void setSubstitutionErrorRate(double substitutionErrorRate) {
		this.substitutionErrorRate = substitutionErrorRate;
	}

	public void setSubstitutionErrorRate(String value) {
		this.setSubstitutionErrorRate((double) OptionValuesDecoder.decode(value, Double.class));
	}

	/**
	 * @return the indelErrorRate
	 */
	public double getIndelErrorRate() {
		return indelErrorRate;
	}

	/**
	 * @param indelErrorRate the indelErrorRate to set
	 */
	public void setIndelErrorRate(double indelErrorRate) {
		this.indelErrorRate = indelErrorRate;
	}

	public void setIndelErrorRate(String value) {
		this.setIndelErrorRate((double) OptionValuesDecoder.decode(value, Double.class));
	}
	
	/**
	 * @return the outFormat
	 */
	public byte getOutFormat() {
		return outFormat;
	}

	/**
	 * @param outFormat the outFormat to set
	 */
	public void setOutFormat(byte outFormat) {
		this.outFormat = outFormat;
	}
	
	public void setOutFormat(String value) {
		this.setOutFormat((byte) OptionValuesDecoder.decode(value, Byte.class));
	}

	/**
	 * @return the genome
	 */
	public ReferenceGenome getGenome() {
		return genome;
	}

	/**
	 * @param genome the genome to set
	 */
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}

	public static void main(String[] args) throws Exception {
		SingleReadsSimulator instance = new SingleReadsSimulator();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String referenceFile = args[i++];
		String outReads = args[i++];

		instance.genome = new ReferenceGenome(referenceFile);
		instance.printInfo(referenceFile, outReads);
		instance.calculate(outReads);
	}

	private void printInfo(String referenceFile, String outReads) {
		System.out.println("Reference path: " + referenceFile);
		System.out.println("Out file path: " + outReads);
		System.out.println("Reference total length:" + genome.getTotalLength());
		System.out.println(
				"Reads:" + numberOfReads + "   ~N(mean: " + meanReadLength + ", sdev: " + stdevReadlength + ")");
		System.out.println("Substitution error rate: " + substitutionErrorRate);
		System.out.println("Indel error rate: " + indelErrorRate);
	}

	private void calculate(String outPath) throws IOException {
		long totalLength = genome.getTotalLength();
		int nSeqs = genome.getNumSequences();
		long[] cumulativeStarts = new long[nSeqs];
		cumulativeStarts[0] = 0;
		for (int i = 1; i < nSeqs; i++) {
			cumulativeStarts[i] = cumulativeStarts[i - 1] + genome.getSequenceByIndex(i - 1).getLength();
		}
		try (PrintStream out = new PrintStream(outPath)) {
			for (int i = 0; i < numberOfReads; i++) {
				int readLength;
				long nextStart;
				QualifiedSequence seq = null;
				int relStart = 0;
				byte reverse = 0;
				String read = null;
				for (int j = 0; j < 100; j++) {
					readLength = (int) (rnd.nextGaussian() * stdevReadlength + meanReadLength);
					long nextLong = rnd.nextLong();
					nextStart = nextLong % (totalLength - readLength);
					if (nextStart < 0)
						continue;

					int idx1 = Arrays.binarySearch(cumulativeStarts, nextStart);

					int sequenceIdx;
					if (idx1 >= 0)
						sequenceIdx = idx1;
					else {
						sequenceIdx = -idx1 - 2;
					}
					if (sequenceIdx < 0)
						System.out.println("Next start: " + nextStart + " idx: " + sequenceIdx);
					seq = genome.getSequenceByIndex(sequenceIdx);
					relStart = (int) (nextStart - cumulativeStarts[sequenceIdx]);
					if (relStart < 0)
						System.out.println("Next start: " + nextStart + " seq: " + seq.getName() + " idx: "
								+ sequenceIdx + " start seq: " + cumulativeStarts[sequenceIdx]);
					int relEnd = relStart + readLength;
					if (relEnd <= seq.getLength()) {
						read = seq.getCharacters().subSequence(relStart, relEnd).toString();
						break;
					}
				}
				if (read == null) {
					// TODO: Warning
					continue;
				}
				if (rnd.nextBoolean()) {
					reverse = 1;
					read = DNAMaskedSequence.getReverseComplement(read);
				}
				String finalRead = generateErrors(read);
				String readId = seq.getName() + "_" + relStart + "_" + reverse;
				if(outFormat == OUT_FORMAT_FASTA) {
					out.println(">" + readId);
					out.println(finalRead);
				} else {
					out.println("@" + readId);
					out.println(finalRead);
					out.println("+");
					out.println(simulateQualities(finalRead.length()));
				}
				
			}
		}
	}

	private String simulateQualities(int length) {
		char [] qualities = new char[length];
		Arrays.fill(qualities, '5');
		return new String(qualities);
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
