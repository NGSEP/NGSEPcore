package ngsep.assembly;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.PriorityQueue;

import ngsep.sequences.FMIndex;

public class GraphBuilderFMIndex implements GraphBuilder {
	private final static int TALLY_DISTANCE = 100;
	private final static int SUFFIX_FRACTION = 25;

	@Override
	public SimplifiedAssemblyGraph buildAssemblyGraph(List<CharSequence> sequences) {
		System.out.println("	sorting sequences");
		long ini = System.currentTimeMillis();
		PriorityQueue<Entry<Integer, Integer>> heap = new PriorityQueue<>(
				(Entry<Integer, Integer> l1, Entry<Integer, Integer> l2) -> l2.getKey() - l1.getKey());
		for (int i = 0; i < sequences.size(); i++)
			heap.add(new SimpleEntry<>(sequences.get(i).length(), i));
		try (PrintStream pr = new PrintStream(new FileOutputStream("ind"))) {
			List<CharSequence> sorted = new ArrayList<>(sequences.size());
			while (!heap.isEmpty()) {
				Entry<Integer, Integer> an = heap.poll();
				sorted.add(sequences.get(an.getValue()));
				pr.println(an.getValue());
			}
			sequences = sorted;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		System.out.println("	sort sequeces: " + (System.currentTimeMillis() - ini) / (double) 1000 + " s");

		System.out.println("	building FMIndexes");
		ini = System.currentTimeMillis();
		FMIndex index = new FMIndex();
		index.loadUnnamedSequences(sequences, TALLY_DISTANCE, SUFFIX_FRACTION);
		System.out.println("	build FMIndexes: " + (System.currentTimeMillis() - ini) / (double) 1000 + " s");

		System.out.println("	indentifing overlaps");
		GraphBuilderOverlapFinder overlapFinder = new GraphBuilderOverlapFinderQueue2();
		overlapFinder.calculate(sequences, index);
		System.out.println("	indentify overlaps: " + (System.currentTimeMillis() - ini) / (double) 1000 + " s");
		return overlapFinder.getGrap();
	}

	public static void main(String[] args) throws Exception {
		System.out.println("¡¡¡ Esta funcion es solo para desarrollo !!!!");
		List<CharSequence> a = Assembler.load(args[0]);
		SimplifiedAssemblyGraph assemblyGraph = (new GraphBuilderFMIndex()).buildAssemblyGraph(a);
		assemblyGraph.save(args[1]);
	}

}
