package ngsep.assembly;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import ngsep.alignments.ReadsAligner;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.ProgressNotifier;
import ngsep.sequences.QualifiedSequence;

public class AssemblyStatisticsCalculator {
	// Logging and progress
	private Logger log = Logger.getLogger(ReadsAligner.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	// Parameters
	private String inputFile = null;
	private String outputFile = null;
	
	
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
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public static void main(String[] args) throws Exception {
		AssemblyStatisticsCalculator instance = new AssemblyStatisticsCalculator();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
		
	}
	public void run () throws IOException {
		logParameters ();
		if(inputFile==null) throw new IOException("The input genome assembly is a required parameter");
		ReferenceGenome genome = new ReferenceGenome(inputFile);
		long[] stats = calculateNStatistics(genome);
		try(PrintStream out = new PrintStream(outputFile)) {
			printStatistics(genome, stats, out);
		}
		
		
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		log.info(os.toString());
		
	}
	
	public static long[] calculateNStatistics(ReferenceGenome genome) {
		List<Integer> lengths = new ArrayList<>(genome.getNumSequences());
		for(QualifiedSequence seq:genome.getSequencesMetadata()) lengths.add(seq.getLength());
		long [] stats = calculateNStatistics(lengths);
		return stats;
	}
	
	public static long[] calculateNStatistics(List<Integer> numbers) {
		long [] answer = new long [10];
		Arrays.fill(answer, 0);
		Collections.sort(numbers,(l1,l2)-> l2-l1);
		long total = 0;
		for(int i:numbers) total+=i;
		answer[0] = 0;
		double current = 0;
		int i=0;
		for(int n:numbers) {
			current += n;
			double p = 10.0*current/total;
			for(;i<answer.length && i<=p;i++) {
				answer[i] = n;
			}
		}
		return answer;
	}
	private static void printStatistics(ReferenceGenome genome, long[] stats, PrintStream out) {
		out.println("Contigs\t"+genome.getNumSequences());
		out.println("Total\t"+genome.getTotalLength());
		out.println("MaxLength\t"+ genome.getLongestSequenceLength());
		printNStatistics(stats, out);
	}
	public static void printNStatistics(long [] stats, PrintStream out) {
		for(int i=1;i<stats.length;i++) {
			out.println("N"+(10*i)+"\t"+stats[i]);
		}
	}
}
