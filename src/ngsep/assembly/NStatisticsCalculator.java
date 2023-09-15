package ngsep.assembly;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequence;

public class NStatisticsCalculator {
	public static void main(String[] args) throws Exception {
		ReferenceGenome genome = new ReferenceGenome(args[0]);
		List<Integer> lengths = new ArrayList<>(genome.getNumSequences());
		for(QualifiedSequence seq:genome.getSequencesMetadata()) lengths.add(seq.getLength());
		long [] stats = calculateNStatistics(lengths);
		printNStatistics(stats, System.out);
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
	public static void printNStatistics(long [] stats, PrintStream out) {
		for(int i=1;i<stats.length;i++) {
			out.println("N"+(10*i)+"\t"+stats[i]);
		}
	}
}
