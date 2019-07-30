package ngsep.assembly;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;
import ngsep.sequences.DNASequence;

public class TestReadGenerator {
	private final static double default_rate_of_changes = 0.07;
	private final static double default_rate_of_cuts = 0.03;
	private final static Random rnd = new Random();

	public static void main(String[] args) throws Exception {
		if (args.length < 8)
			throw new Exception("Faltan Parametros para generar los archivos");

		String path = args[0];
		String pathlect = args[1];
		int numberOfSequences = Integer.parseInt(args[2]);
		int[] odist = new int[] { Integer.parseInt(args[3]), Integer.parseInt(args[4]) };
		int numberOfreads = Integer.parseInt(args[5]);
		int[] dist = new int[] { Integer.parseInt(args[6]), Integer.parseInt(args[7]) };
		double rateChanges = (args.length > 8) ? Double.parseDouble(args[8]) : default_rate_of_changes;
		double rateIndels = (args.length > 9) ? Double.parseDouble(args[9]) : default_rate_of_cuts;
		String pathGraph = (args.length > 10) ? args[10] : null;
		int minOveralp = (args.length > 11) ? Integer.parseInt(args[11]) : 0;

		printInfo(path, pathlect, numberOfSequences, odist, numberOfreads, dist, rateChanges, rateIndels, pathGraph,
				minOveralp);
		calculate(path, pathlect, numberOfSequences, odist, numberOfreads, dist, rateChanges, rateIndels, pathGraph,
				minOveralp);
	}

	private static void printInfo(String path, String pathlect, int numberOfSequences, int[] odist, int numberOfreads,
			int[] dist, double rateChanges, double rateIndels, String pathGraph, int minOveralp) {
		System.out.println("Reference path: " + path);
		System.out.println("Samples path: " + pathlect);
		System.out.println("Reference:" + numberOfSequences + "   ~N(mean: " + odist[0] + ", sdev: " + odist[1] + ")");
		System.out.println("Samples:" + numberOfreads + "   ~N(mean: " + dist[0] + ", sdev: " + dist[1] + ")");
		System.out.println("Rate of changes: " + rateChanges);
		System.out.println("Rate of indels: " + rateIndels);
		if (pathGraph != null) {
			System.out.println();
			System.out.println("Graph Path: " + pathGraph);
			System.out.println("Minimun overlap in the graph: " + minOveralp);
		}
	}

	private static void calculate(String path, String pathlect, int numberOfSequences, int[] odist, int numberOfreads,
			int[] dist, double rateChanges, double rateIndels, String pathGraph, int minOveralp) throws IOException {

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
		Assembler.exportToFile(pathlect, "lect", lects);

		if (pathGraph != null) {
			SimplifiedAssemblyGraph sag = getGraph(lects, minOveralp, refLects, posLects, lenLects, revLects);
			sag.save(pathGraph);
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

	private static SimplifiedAssemblyGraph getGraph(String[] lects, int minOveralp, int[] refLects, int[] posLects,
			int[] lenLects, boolean[] revLects) throws FileNotFoundException {
		SimplifiedAssemblyGraph sag = new SimplifiedAssemblyGraph(Arrays.asList(lects));

		Map<Integer, List<Lect>> a = new HashMap<>();
		for (int i = 0; i < refLects.length; i++)
			a.computeIfAbsent(refLects[i], (x) -> new ArrayList<>())
					.add(new Lect(i, lenLects[i], revLects[i], posLects[i]));

		System.out.println("in");
		for (List<Lect> list : a.values())
			Collections.sort(list, (Lect x, Lect y) -> {
				int ans = x.pos - y.pos;
				if (ans != 0)
					return ans;
				return y.len - x.len;
			});

		System.out.println("in");
		try (PrintStream pr = new PrintStream(new FileOutputStream("kkkk"))) {
			for (List<Lect> list : a.values()) {
				for (int i = 0; i < list.size() - 1; i++) {
					Lect lect1 = list.get(i);
					pr.println(lect1.id + ") " + lect1.pos + " " + lect1.len);
					for (int j = i + 1; j < list.size() && list.get(j).pos < lect1.pos + lect1.len; j++) {
						Lect lect2 = list.get(j);

						int relativePos = lect2.pos - lect1.pos;
						if (relativePos + lect2.len > lect1.len)
							// lect1 -> lect2
							sag.addEdge((lect1.id << 1) + (lect1.rev ? 0 : 1), (lect2.id << 1) + (lect2.rev ? 1 : 0),
									lect1.len - relativePos, 1);
						else {
							// lect2 into lect1
							boolean reversed = lect1.rev ^ lect2.rev;
							relativePos = (lect1.rev) ? lect1.len - lect2.len - relativePos : relativePos;
							sag.addEmbedded(lect1.id, lect2.id, relativePos, reversed, 1);
						}
					}
				}
			}
		}
		System.out.println("out");

		sag.removeAllEmbeddedsIntoGraph();
		System.out.println("-------------Graph Properties----------");
		sag.printInfo();
		System.out.println("---------------------------------------");

		return sag;
	}

	private static class Lect {
		public int id, len, pos;
		public boolean rev;

		public Lect(int id, int len, boolean rev, int pos) {
			this.id = id;
			this.len = len;
			this.rev = rev;
			this.pos = pos;
		}
	}

}
