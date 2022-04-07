package ngsep.math;

import java.util.HashMap;
import java.util.Map;

public class ShannonEntropyCalculator {
	private static double [][] precalculatedTerms = new double[100][100];
	private static double LOG2_BASE10 = Math.log(2);
	static {
		initializeCache();
	}
	private static void initializeCache() {
		for(int i=0;i<100;i++) {
			int n = i+1;
			for(int j=0;j<=i;j++) {
				double count = j+1;
				precalculatedTerms[i][j] = calculateTerm(count,n);
				//if(n==5) System.out.println("count: "+count+" n: "+n+" term: "+precalculatedTerms[i][j] );
			}
		}
	}
	private static double calculateTerm(double count, int n) {
		double p = count / n;
		double invp = 1/p;
		return (p*Math.log(invp)/LOG2_BASE10);
	}
	public static double calculateEntropy(CharSequence sequence) {
		Map<Character,Integer> charCounts = new HashMap<>();
		int n = sequence.length();
		if(n==0) return 0;
		for(int i=0;i<n;i++) {
			charCounts.compute(sequence.charAt(i), (k,v)->(v==null)?1:v+1);
		}
		double answer = 0;
		for(int count:charCounts.values()) {
			if(n<100) answer+=precalculatedTerms[n-1][count-1];
			else answer+=calculateTerm(count, n);
			//System.out.println("Next count: "+count+" term: "+ precalculatedTerms[n-1][count-1]+" answer: "+answer);
		}
		return answer;
	}
	public static void main(String[] args) {
		System.out.println(""+calculateEntropy("KRGGR"));
		System.out.println(""+calculateEntropy("KRGRR"));
		System.out.println(""+calculateEntropy("KRTFC"));
	}
}
