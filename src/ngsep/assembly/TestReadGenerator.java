package ngsep.assembly;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;

public class TestReadGenerator {
	private final static double Rate_of_changes = 0.00;
	private final static double Rate_of_cuts = 0.00;
	private static Random rnd = new Random();

	public TestReadGenerator(String path, String pathlect, int numberOfSequences, int[] sequencesdist,
			int numberOfreads, int[] distribution, double... tasas) throws FileNotFoundException {
		double rateChages = (tasas.length > 1) ? tasas[0] : Rate_of_changes;
		double rateCuts = (tasas.length > 2) ? tasas[1] : Rate_of_cuts;

		String[] ref = getRef(numberOfSequences, sequencesdist);
		String[] lects = new String[numberOfreads];

		for (int i = 0; i < lects.length; i++)
			lects[i] = getLect(ref, distribution, rateCuts, rateChages);

		print(path, ref, "ref");
		print(pathlect, lects, "lect");
	}

	private static String getLect(String[] ref, int[] dist, double Ncuts, double Nchng) {
		int i = rnd.nextInt(ref.length);
		int length = (int) (dist[0] + ((rnd.nextBoolean()) ? 1 : -1) * dist[1] * rnd.nextDouble());
		int pos = rnd.nextInt(ref[i].length() - length);
		String str = ref[i];

		StringBuilder strb = new StringBuilder();
		int j = 0;
		for (int x : indexTo(length, (int) (length * Ncuts))) {
			while (j < x) {
				strb.append(str.charAt(pos + j));
				j++;
			}
			j++;
		}
		while (j < length) {
			strb.append(str.charAt(pos + j));
			j++;
		}

		char[] arr = strb.toString().toCharArray();
		for (int x : indexTo(arr.length, (int) (Nchng * arr.length))) {
			char k = next();
			while (k == arr[x])
				k = next();
			arr[x] = k;
		}
		return new String(arr);
	}

	private static TreeSet<Integer> indexTo(int N, int cuts) {
		TreeSet<Integer> ans = new TreeSet<Integer>();
		for (int i = 0; i < cuts; i++) {
			int j = rnd.nextInt(N);
			while (ans.contains(j))
				j = rnd.nextInt(N);
			ans.add(j);
		}
		return ans;
	}

	private static String[] getRef(int N, int[] dist) {
		String[] ans = new String[N];
		for (int i = 0; i < ans.length; i++) {
			int length = (int) (dist[0] + ((rnd.nextBoolean()) ? 1 : -1) * dist[1] * rnd.nextDouble());
			ans[i] = createSequence(length);
		}
		return ans;
	}

	private static String createSequence(int length) {
		StringBuilder str = new StringBuilder();
		for (int i = 0; i < length; i++) {
			str.append(next());
		}
		return str.toString();
	}

	private static char next() {
		return DNASequence.BASES_ARRAY[rnd.nextInt(DNASequence.BASES_ARRAY.length)].charAt(0);
	}

	private static void print(String file, String[] sec, String name) throws FileNotFoundException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> list = new ArrayList<QualifiedSequence>();
		int i = 1;
		for (String str : sec)
			list.add(new QualifiedSequence(name + "_" + (i++), str));
		try (PrintStream pr = new PrintStream(new FileOutputStream(file))) {
			handler.saveSequences(list, pr, 1000);
		}
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 8)
			throw new Exception("faltan parametros");

		String path = args[0];
		String pathlect = args[1];
		int numberOfSequences = Integer.parseInt(args[2]);
		int[] odist = new int[] { Integer.parseInt(args[3]), Integer.parseInt(args[4]) };
		int numberOfreads = Integer.parseInt(args[5]);
		int[] dist = new int[] { Integer.parseInt(args[6]), Integer.parseInt(args[7]) };

		switch (args.length) {
		case 8:
			new TestReadGenerator(path, pathlect, numberOfSequences, odist, numberOfreads, dist);
			break;

		case 9:
			new TestReadGenerator(path, pathlect, numberOfSequences, odist, numberOfreads, dist,
					Double.parseDouble(args[5]));
			break;

		case 10:
			new TestReadGenerator(path, pathlect, numberOfSequences, odist, numberOfreads, dist,
					Double.parseDouble(args[5]), Double.parseDouble(args[6]));
			break;

		default:
			throw new Exception("to many parameters");
		}
	}
}
