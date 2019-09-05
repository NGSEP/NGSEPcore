package ngsep.simulation;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Random;
import java.util.TreeSet;
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
	private ProgressNotifier progressNotifier=null;
	public final static int DEF_NUM_READS = 30000;
	public final static int DEF_MEAN_READ_LENGTH = 10000;
	public final static int DEF_STDEV_READ_LENGTH = 2000;
	public final static double DEF_SUBSTITUTION_ERROR_RATE = 0.02;
	public final static double DEF_INDEL_ERROR_RATE = 0.01;
	private final static Random rnd = new Random();
	private int numberOfReads = DEF_NUM_READS;
	private int meanReadLength = DEF_MEAN_READ_LENGTH;
	private int stdevReadlength = DEF_STDEV_READ_LENGTH;
	private double substitutionErrorRate = DEF_SUBSTITUTION_ERROR_RATE;
	private double indelErrorRate = DEF_INDEL_ERROR_RATE;
	
	private ReferenceGenome genome;
	
	
	
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

	/**
	 * @return the stdevReadlength
	 */
	public int getStdevReadlength() {
		return stdevReadlength;
	}

	/**
	 * @param stdevReadlength the stdevReadlength to set
	 */
	public void setStdevReadlength(int stdevReadlength) {
		this.stdevReadlength = stdevReadlength;
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
		System.out.println("Reads:" + numberOfReads + "   ~N(mean: " + meanReadLength + ", sdev: " + stdevReadlength + ")");
		System.out.println("Substitution error rate: " + substitutionErrorRate);
		System.out.println("Indel error rate: " + indelErrorRate);
	}

	private void calculate(String outPath) throws IOException {
		long totalLength = genome.getTotalLength();
		int nSeqs = genome.getNumSequences();
		long [] cumulativeStarts = new long [nSeqs];
		cumulativeStarts[0] = 0;
		for(int i=1;i<nSeqs;i++) {
			cumulativeStarts[i] = cumulativeStarts[i-1]+genome.getSequenceByIndex(i-1).getLength(); 
		}
		try (PrintStream out = new PrintStream(outPath)) {
			for (int i = 0; i < numberOfReads; i++) {
				int readLength;
				long nextStart;
				QualifiedSequence seq=null;
				int relStart=0;
				byte reverse=0;
				String read = null;
				for(int j=0;j<100;j++) {
					readLength = (int) (rnd.nextGaussian() * stdevReadlength + meanReadLength);
					nextStart = rnd.nextLong()%(totalLength-readLength);
					int idx1 = Arrays.binarySearch(cumulativeStarts, nextStart);
					
					int sequenceIdx;
					if(idx1>=0) sequenceIdx = idx1;
					else {
						//TODO: Choose actual chromosome
						sequenceIdx = 0;
					}
					
					seq = genome.getSequenceByIndex(sequenceIdx);
					relStart = (int) (nextStart-cumulativeStarts[sequenceIdx]);
					int relEnd = relStart+readLength; 
					if(relEnd<=seq.getLength()) {
						read = seq.getCharacters().subSequence(relStart, relEnd).toString();
						break;
					}
				}
				if(read==null) {
					//TODO: Warning
					continue;
				}
				if(rnd.nextBoolean()) {
					reverse = 1;
					read = DNAMaskedSequence.getReverseComplement(read);
				}
				String finalRead = generateErrors(read);
				String readId = seq.getName()+"_"+relStart+"_"+reverse;
				out.println(">"+readId);
				out.println(finalRead);
			}
		}
	}
	
	private String generateErrors(String read) {
		String alphabet = DNASequence.BASES_STRING;
		int len = read.length();
		Integer[] cuts = uniformSorted((int) (indelErrorRate * len), len);
		char[] ans = new char[len - cuts.length];

		// copy without cuts
		int i = 0, j = 0;
		for (int x : cuts) {
			while (i < x)
				ans[j++] = read.charAt((i++));
			i++;
		}
		while (j < ans.length)
			ans[j++] = read.charAt((i++));

		// change letters
		for (int x : uniformSorted((int) (ans.length * substitutionErrorRate), ans.length)) {
			char k = alphabet.charAt(rnd.nextInt(alphabet.length()));
			while (k == ans[x])
				k = alphabet.charAt(rnd.nextInt(alphabet.length()));
			ans[x] = k;
		}
		return new String(ans);
	}

	public static Integer[] uniformSorted(int samples, int N) {
		TreeSet<Integer> ans = new TreeSet<Integer>();
		for (int i = 0; i < samples; i++) {
			int j = rnd.nextInt(N);
			while (ans.contains(j))
				j = rnd.nextInt(N);
			ans.add(j);
		}
		return ans.toArray(new Integer[0]);
	}

}
