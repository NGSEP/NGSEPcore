package ngsep.assembly;

import java.text.DecimalFormat;
import java.util.concurrent.Callable;

public class TimeUtilities {
	private static DecimalFormat df = new DecimalFormat("00.00");

	public static void timeGroup(String str, Runnable runnable) {
		System.out.println(str);
		long ini = System.currentTimeMillis();
		runnable.run();
		System.out.println(str + " : " + (System.currentTimeMillis() - ini) + " ms");
	}

	public static <T extends Object> T timeGroup(String str, Callable<T> callable) throws RuntimeException {
		try {
			System.out.println(str);
			long ini = System.currentTimeMillis();
			T t = callable.call();
			System.out.println(str + (System.currentTimeMillis() - ini) + " ms");
			return t;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	public static void timeIt(String str, Runnable runnable) {
		System.out.print(str);
		long ini = System.currentTimeMillis();
		runnable.run();
		System.out
				.println("\r" + str + " : " + progressBar(1, 1, 50) + "	" + (System.currentTimeMillis() - ini) + " ms");
	}

	public static <T extends Object> T timeIt(String str, Callable<T> callable) throws RuntimeException {
		try {
			System.out.print(str);
			long ini = System.currentTimeMillis();
			T t = callable.call();
			System.out.println("\r" + str + " : " + progressBar(1, 1, 50).toString() + "	"
					+ (System.currentTimeMillis() - ini) + " ms");
			return t;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	public static void progress(String str, double act, double max) {
		System.out.print("\r" + str + " : " + progressBar(act, max, 50));
	}

	private static CharSequence progressBar(double act, double max, int N) {
		StringBuilder string = new StringBuilder();
		string.append("|");
		int i = 0;
		while (i < N * (act / max)) {
			string.append("█");
			i++;
		}
		while (i < N) {
			string.append(" ");
			i++;
		}
		string.append("| ");
		string.append(df.format(100 * act / max));
		string.append("%");
		return string;
	}
}
