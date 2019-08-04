package ngsep.assembly;

import java.text.DecimalFormat;
import java.util.concurrent.Callable;

public class TimeUtilities {
	private static DecimalFormat df = new DecimalFormat("00.00");

	/**
	 * used to print a group of timeIt elements
	 * 
	 * @param str      the name
	 * @param runnable the code
	 */
	public static void timeGroup(String str, Runnable runnable) {
		System.out.println(str);
		long ini = System.currentTimeMillis();
		runnable.run();
		System.out.println(str + " : " + (System.currentTimeMillis() - ini) + " ms");
	}

	/**
	 * used to print a group of timeIt elements
	 * 
	 * @param <T>      something
	 * @param str      the name
	 * @param callable the code
	 * @return whatever that returns the code inside
	 * @throws RuntimeException cast every exception inside the code
	 */
	public static <T extends Object> T timeGroup(String str, Callable<T> callable) throws RuntimeException {
		try {
			System.out.println(str);
			long ini = System.currentTimeMillis();
			T t = callable.call();
			System.out.println(str + " : " + (System.currentTimeMillis() - ini) + " ms");
			return t;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * time an specific event inside this code you can call progress to report
	 * 
	 * @param str      the name of the event
	 * @param runnable the code
	 */
	public static void timeIt(String str, Runnable runnable) {
		System.out.print(str);
		long ini = System.currentTimeMillis();
		runnable.run();
		System.out
				.println("\r" + str + " : " + progressBar(1, 1, 50) + "	" + (System.currentTimeMillis() - ini) + " ms");
	}

	/**
	 * time an specific event inside this code you can call progress to report
	 * 
	 * @param <T>      something
	 * @param str      the name of the event
	 * @param callable the code
	 * @return whatever that returns the code inside
	 * @throws RuntimeException cast every exception inside the code
	 */
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

	/**
	 * change the progress of timeIt
	 * 
	 * @param str the name (the same of the timeIt)
	 * @param act actual value
	 * @param max objective value
	 */
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
