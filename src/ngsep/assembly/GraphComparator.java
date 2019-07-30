package ngsep.assembly;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import ngsep.sequences.DNAMaskedSequence;

import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class GraphComparator implements Comparator<SimplifiedAssemblyGraph> {

	@Override
	public int compare(SimplifiedAssemblyGraph ref, SimplifiedAssemblyGraph lec) {
		ref.printInfo();
		lec.printInfo();

		int[] lecttoRef = new int[ref.getSequences().size()];
		int[] reftoLect = new int[ref.getSequences().size()];

		try (BufferedReader br = new BufferedReader(new FileReader("ind"))) {
			int j = 0;
			String line = br.readLine();
			while (line != null) {
				int i = Integer.valueOf(line);
				reftoLect[i] = j;
				lecttoRef[j] = i;
				j++;
				line = br.readLine();
			}
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		System.out.println(ref.getSequences().size());
		System.out.println(lec.getSequences().size());

		try (PrintStream pr = new PrintStream(new FileOutputStream("test.txt"))) {

			Set<Integer> refEmb = new HashSet<>();
			for (Map<Integer, Embedded> a : ref.getEmbbeded().values()) {
				refEmb.addAll(a.keySet());

			}
			Set<Integer> lectEmb = new HashSet<>();
			System.out.println("False Positives");
			for (Entry<Integer, Map<Integer, Embedded>> a : lec.getEmbbeded().entrySet()) {
				lectEmb.addAll(a.getValue().keySet());
				for (int i : a.getValue().keySet()) {
					
						pr.println(a.getKey() + "     " + i + "   (" + a.getValue().get(i).getPos() + ")  "
								+ a.getValue().get(i).getRate());
						pr.println(lecttoRef[a.getKey()] + "     " + lecttoRef[i] + "   ("
								+ a.getValue().get(i).getPos() + ")  " + a.getValue().get(i).getRate());
						pr.println(refEmb.contains(lecttoRef[i]));
						pr.println(lec.getSequences().get(a.getKey()));
						StringBuilder sb = new StringBuilder();
						for (int ll = 0; ll < a.getValue().get(i).getPos(); ll++) {
							sb.append("-");
						}
						String as = lec.getSequences().get(i).toString();
						sb.append(a.getValue().get(i).isReversed() ? DNAMaskedSequence.getReverseComplement(as) : as);
						pr.println(sb);
				}
			}

			System.out.println(refEmb.size());
			System.out.println(lectEmb.size());

			int trueP = 0, falseP = 0;
			for (int i : lectEmb) {
				if (refEmb.contains(lecttoRef[i]))
					trueP++;
				else
					falseP++;
			}

			int trueN = 0;
			for (int i : refEmb) {
				if (!lectEmb.contains(reftoLect[i]))
					trueN++;
			}
			int falseN = (ref.getSequences().size() - ref.amuontOfEmbeddedSequences()) - falseP;

			System.out.println("ref");
			System.out.println("          true|false");
			System.out.println("positive= " + trueP + "|" + falseP);
			System.out.println("negative= " + falseN + "|" + trueN);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return 0;
	}

	public static void main(String[] args) throws FileNotFoundException, ClassNotFoundException, IOException {
		SimplifiedAssemblyGraph ref = new SimplifiedAssemblyGraph(args[0]);
		SimplifiedAssemblyGraph sag = new SimplifiedAssemblyGraph(args[1]);
		(new GraphComparator()).compare(ref, sag);
	}

}
