package ngsep.assembly;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import ngsep.sequences.DNASequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;

public class TestReadGenerator {
	private final static double Rate_of_changes = 0.07;
	private final static double Rate_of_cuts = 0.03;
	private static Random rnd = new Random();

	public TestReadGenerator(String path, String pathlect, String pathLecRef, int numberOfSequences, int[] odist,
			int numberOfreads, int[] dist, double rateChanges, double rateIndels) throws FileNotFoundException {

		String[] ref = getRef(numberOfSequences, odist);
		String[] lects = new String[numberOfreads];
		String[] lectsRef = new String[numberOfreads];

		for (int j = 0; j < lects.length; j++) {
			int i = rnd.nextInt(ref.length);
			int length = (int) (dist[0] + ((rnd.nextBoolean()) ? 1 : -1) * dist[1] * rnd.nextDouble());
			int pos = rnd.nextInt(ref[i].length() - length);
			lects[j] = getLect(ref[i], pos, length, rateIndels, rateChanges);
			lectsRef[j] = ref[i].substring(pos, pos + length);
		}

		print(path, ref, "ref");
		print(pathlect, lects, "lect");
		print(pathLecRef, lectsRef, "lectRef");
	}

	private static String getLect(String str, int pos, int length, double Ncuts, double Nchng) {
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

		String ans = new String(arr);
		if (rnd.nextBoolean())
			ans = DNASequence.getReverseComplement(ans);
		return ans;
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
		char[] str = new char[length];
		for (int i = 0; i < length; i++)
			str[i] = next();
		return new String(str);
	}

	private static char next() {
		return DNASequence.BASES_ARRAY[rnd.nextInt(DNASequence.BASES_ARRAY.length)].charAt(0);
	}

	private static void print(String file, String[] sec, String name) throws FileNotFoundException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> list = new ArrayList<QualifiedSequence>();
		int i = 1;
		for (String str : sec) {
			if (str == null)
				System.out.println(i);
			list.add(new QualifiedSequence(name + "_" + (i++), str));
		}
		try (PrintStream pr = new PrintStream(new FileOutputStream(file))) {
			handler.saveSequences(list, pr, 1000);
		}
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 8)
			throw new Exception("Faltan Parametros para generar los archivos");

		String path = args[0];
		String pathlect = args[1];
		String pathLecRef = addToPathName(pathlect, "_ref");
		int numberOfSequences = Integer.parseInt(args[2]);
		int[] odist = new int[] { Integer.parseInt(args[3]), Integer.parseInt(args[4]) };
		int numberOfreads = Integer.parseInt(args[5]);
		int[] dist = new int[] { Integer.parseInt(args[6]), Integer.parseInt(args[7]) };
		double rateChanges = (args.length > 8) ? Double.parseDouble(args[8]) : Rate_of_changes;
		double rateIndels = (args.length > 9) ? Double.parseDouble(args[9]) : Rate_of_cuts;
		printInfo(path, pathlect, pathLecRef, numberOfSequences, odist, numberOfreads, dist, rateChanges, rateIndels);

		new TestReadGenerator(path, pathlect, pathLecRef, numberOfSequences, odist, numberOfreads, dist, rateChanges,
				rateIndels);
	}

	private static void printInfo(String path, String pathlect, String pathLecRef, int numberOfSequences, int[] odist,
			int numberOfreads, int[] dist, double rateChanges, double rateIndels) {
		System.out.println("Ruta archivo de referencia: " + path);
		System.out.println("Ruta lecturas de referencia: " + pathLecRef);
		System.out.println("Ruta levturas: " + pathlect);
		System.out.println("Referencia");
		System.out.println("Secuencias:" + numberOfSequences + "   ~Uniforme(" + (odist[0] - odist[1]) + ", "
				+ (odist[0] + odist[1]) + ")");
		System.out.println("Lecturas");
		System.out.println("Secuencias:" + numberOfreads + "   ~Uniforme(" + (dist[0] - dist[1]) + ", "
				+ (dist[0] + dist[1]) + ")");
		System.out.println("Rate of changes: " + rateChanges);
		System.out.println("Rate of indels: " + rateIndels);
	}

	private static String addToPathName(String foo, String extra) {
		String[] aux = foo.split("\\.");
		StringBuilder sb = new StringBuilder();
		sb.append(aux[0]);
		if (aux.length > 1) {
			for (int i = 1; i < aux.length - 1; i++)
				sb.append("." + aux[i]);
			sb.append(extra);
			sb.append("." + aux[aux.length - 1]);
			return sb.toString();
		}
		return foo + extra;
	}
}
