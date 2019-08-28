package ngsep.simulation;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.TreeSet;
import java.util.AbstractMap.SimpleEntry;
import java.util.Map.Entry;

import ngsep.assembly.Assembler;
import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;

public class SingleReadsGenerator {
	private final static double DEF_SUBSTITUTION_ERROR_RATE = 0.02;
	private final static double DEF_INDEL_ERROR_RATE = 0.01;
	private final static Random rnd = new Random();

	public static void main(String[] args) throws Exception {
		String path = args[0];
		String pathlect = args[1];
		int numberOfSequences = Integer.parseInt(args[2]);
		int[] odist = new int[] { Integer.parseInt(args[3]), Integer.parseInt(args[4]) };
		int numberOfreads = Integer.parseInt(args[5]);
		int[] dist = new int[] { Integer.parseInt(args[6]), Integer.parseInt(args[7]) };
		double rateChanges = (args.length > 8) ? Double.parseDouble(args[8]) : DEF_SUBSTITUTION_ERROR_RATE;
		double rateIndels = (args.length > 9) ? Double.parseDouble(args[9]) : DEF_INDEL_ERROR_RATE;

		printInfo(path, pathlect, numberOfSequences, odist, numberOfreads, dist, rateChanges, rateIndels);
		calculate(path, pathlect, numberOfSequences, odist, numberOfreads, dist, rateChanges, rateIndels);
	}

	private static void printInfo(String path, String pathlect, int numberOfSequences, int[] odist, int numberOfreads,
			int[] dist, double rateChanges, double rateIndels) {
		System.out.println("Reference path: " + path);
		System.out.println("Samples path: " + pathlect);
		System.out.println("Reference:" + numberOfSequences + "   ~N(mean: " + odist[0] + ", sdev: " + odist[1] + ")");
		System.out.println("Samples:" + numberOfreads + "   ~N(mean: " + dist[0] + ", sdev: " + dist[1] + ")");
		System.out.println("Rate of changes: " + rateChanges);
		System.out.println("Rate of indels: " + rateIndels);
	}

	private static void calculate(String path, String pathlect, int numberOfSequences, int[] odist, int numberOfreads,
			int[] dist, double rateChanges, double rateIndels) throws IOException {

		String[] refs = new String[numberOfSequences];
		int[] lenRefs = normalValues(numberOfSequences, odist[0], odist[1]);
		for (int i = 0; i < numberOfSequences; i++)
			refs[i] = randomsequence(lenRefs[i], DNASequence.BASES_STRING);

		String[] lects = new String[numberOfreads];
		int[] lenLects = normalValues(numberOfreads, dist[0], dist[1]);
		int[] refLects = new int[numberOfreads];
		int[] posLects = new int[numberOfreads];
		boolean[] revLects = new boolean[numberOfreads];
		for (int i = 0; i < numberOfreads; i++) {
			refLects[i] = rnd.nextInt(refs.length);
			posLects[i] = rnd.nextInt(refs[refLects[i]].length() - lenLects[i]);
			revLects[i] = rnd.nextBoolean();
			String seq = subSequ(refs[refLects[i]], posLects[i], lenLects[i], rateChanges, rateIndels,
					DNASequence.BASES_STRING);
			lects[i] = revLects[i] ? DNASequence.getReverseComplement(seq) : seq;
		}

		Assembler.exportToFile(path, "ref", refs);
		exportToFileSorted(pathlect, lects, refLects, posLects, revLects);
	}

	private static void exportToFileSorted(String pathlect, String[] lects, int[] refLects, int[] posLects,
			boolean[] revLects) throws FileNotFoundException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> list = new ArrayList<QualifiedSequence>();
		PriorityQueue<Entry<Integer, Integer>> heap = new PriorityQueue<>(
				(Entry<Integer, Integer> l1, Entry<Integer, Integer> l2) -> l2.getKey() - l1.getKey());
		for (int i = 0; i < lects.length; i++)
			heap.add(new SimpleEntry<>(lects[i].length(), i));

		int j = 0;
		while (!heap.isEmpty()) {
			int ind = heap.poll().getValue();
			list.add(new QualifiedSequence(
					"lect" + (++j) + "-" + refLects[ind] + "-" + posLects[ind] + "-" + revLects[ind], lects[ind]));
		}

		try (PrintStream pr = new PrintStream(new FileOutputStream(pathlect))) {
			handler.saveSequences(list, pr, 1000);
		}
	}

	public static int[] normalValues(int samples, int mean, int strDesv) {
		int[] ans = new int[samples];
		for (int i = 0; i < ans.length; i++) {
			int aux = (int) (rnd.nextGaussian() * strDesv + mean);
			ans[i] = (aux > 0) ? aux : 1;
		}
		return ans;
	}

	public static String randomsequence(int len, String alphabet) {
		char[] str = new char[len];
		for (int i = 0; i < len; i++)
			str[i] = alphabet.charAt(rnd.nextInt(alphabet.length()));
		return new String(str);
	}

	private static String subSequ(String ref, int pos, int len, double rateChanges, double rateIndels,
			String alphabet) {
		Integer[] cuts = uniformSorted((int) (rateIndels * len), len);
		char[] ans = new char[len - cuts.length];

		// copy without cuts
		int i = 0, j = 0;
		for (int x : cuts) {
			while (i < x)
				ans[j++] = ref.charAt(pos + (i++));
			i++;
		}
		while (j < ans.length)
			ans[j++] = ref.charAt(pos + (i++));

		// change letters
		for (int x : uniformSorted((int) (ans.length * rateChanges), ans.length)) {
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
