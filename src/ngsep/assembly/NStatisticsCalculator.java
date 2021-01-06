package ngsep.assembly;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class NStatisticsCalculator {
	public static int[] calculateNStatistics(List<Integer> numbers) {
		int[] answer = new int [10];
		Arrays.fill(answer, 0);
		Collections.sort(numbers,(l1,l2)-> l2-l1);
		int total = 0;
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
	public static void printNStatistics(int [] stats, PrintStream out) {
		for(int i=1;i<stats.length;i++) {
			out.println("N"+(10*i)+"\t"+stats[i]);
		}
	}
}
