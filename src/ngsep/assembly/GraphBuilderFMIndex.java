package ngsep.assembly;

import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import ngsep.alignments.ReadAlignment;
import ngsep.sequences.FMIndex;

import static ngsep.assembly.TimeUtilities.progress;
import static ngsep.assembly.TimeUtilities.timeGroup;
import static ngsep.assembly.TimeUtilities.timeIt;

public class GraphBuilderFMIndex implements GraphBuilder {
	private final static int TALLY_DISTANCE = 100;// *64
	private final static int SUFFIX_FRACTION = 16;

	private FMIndex index;
	private KmerIterator kmerIterator;
	private List<CharSequence> sequences;
	private AssemblyConfiguration config;

	private SimplifiedAssemblyGraph assemblyGraph;
	private Map<Integer, List<int[]>> hits;

	@Override
	public AssemblyGraph buildAssemblyGraph(List<CharSequence> sequences) {
		return buildSimplifiedAssemblyGraph(sequences).getAssemblyGraph();
	}

	public SimplifiedAssemblyGraph buildSimplifiedAssemblyGraph(List<CharSequence> sequences) {
		this.sequences = sequences;
		if (!isSorted(sequences)) {
			timeIt("Sort sequences", () -> Collections.sort(sequences, (l1, l2) -> l2.length() - l1.length()));
		}

		index = timeGroup("    Build FMindex", () -> {
			FMIndex ans = new FMIndex();
			ans.loadUnnamedSequences(sequences, TALLY_DISTANCE, SUFFIX_FRACTION);
			return ans;
		});

		timeGroup("    Indentify overlaps", () -> {
			kmerIterator = new KmerIterator(config);
			assemblyGraph = new SimplifiedAssemblyGraph(sequences);
			hits = new HashMap<Integer, List<int[]>>(sequences.size());

			printRates();
			timeIt("      Find overlaps ", () -> findOverlapsAndEmbedded());
			timeIt("      Clean Graph", () -> assemblyGraph.removeAllEmbeddedsIntoGraph());
			timeIt("      Extrapolate Aligns", () -> assemblyGraph.ExtrapolateAligns());
		});
		return assemblyGraph;
	}

	/**
	 * @return the config
	 */
	public AssemblyConfiguration getConfig() {
		return config;
	}

	/**
	 * @param config the config to set
	 */
	public void setConfig(AssemblyConfiguration config) {
		this.config = config;
	}

	/**
	 * @param sequences
	 */
	private boolean isSorted(List<CharSequence> sequences) {
		for (int i = 0; i < sequences.size() - 1; i++)
			if (sequences.get(i).length() < sequences.get(i + 1).length())
				return false;
		return true;
	}

	public void printRates() {
		System.out.println("      --------------------------------");
		System.out.println("      SEARCH_KMER_LENGTH: " + config.overlap().getKmerLength());
		System.out.println("      SEARCH_KMER_DISTANCE: " + config.overlap().getKmerDistance());
		System.out.println("      MAX_KMER_DES: " + config.overlap().getMaxKmerDiff());
		System.out.println("      MIN_COVER_RATE: " + config.overlap().getMinKmerCoverRate());
		System.out.println("      --------------------------------");
	}

	public void findOverlapsAndEmbedded() {
		HitsAligner aligner = new TreesHitAligner(assemblyGraph, config, sequences);
		for (int seqId = 0, excp = 0; seqId < sequences.size(); seqId++) {
			progress("      Find overlaps ", seqId + assemblyGraph.amuontOfEmbeddedSequences() - excp,
					sequences.size());
			if (assemblyGraph.isEmbedded(seqId)) {
				excp++;
				continue;
			}

			calculateHits(seqId, kmerIterator.positiveStrand(sequences.get(seqId)));
			for (Entry<Integer, List<int[]>> entry : hits.entrySet())
				if (entry.getValue().size() > 1)
					aligner.Aling(seqId, entry.getKey(), false, entry.getValue());

			calculateHits(seqId, kmerIterator.negativeStrand(sequences.get(seqId)));
			for (Entry<Integer, List<int[]>> entry : hits.entrySet())
				if (entry.getValue().size() > 1)
					aligner.Aling(seqId, entry.getKey(), true, entry.getValue());
		}
	}

	private void calculateHits(int id_Ref, Iterable<Entry<Integer, String>> kmerIters) {
		hits.clear();
		for (Entry<Integer, String> entry : kmerIters) {
			int pos_Ref = entry.getKey();
			for (ReadAlignment readAlignment : index.search(entry.getValue())) {
				int id_Lec = Integer.parseInt(readAlignment.getSequenceName());
				int pos_Lec = readAlignment.getFirst();

				if (id_Ref < id_Lec && !assemblyGraph.isEmbedded(id_Lec))
					hits.computeIfAbsent(id_Lec, x -> new LinkedList<int[]>()).add(new int[] { pos_Ref, pos_Lec });
			}
		}
	}

	public static void main(String[] args) throws Exception {
		List<CharSequence> sequences = timeIt("Load the sequences", () -> Assembler.load(args[0]));

		SimplifiedAssemblyGraph assemblyGraph = timeGroup("Build overlap Graph", () -> {
			AssemblyConfiguration config = (args.length > 3)
					? new AssemblyConfiguration(Double.valueOf(args[2]), Double.valueOf(args[3]))
					: new AssemblyConfiguration();
			GraphBuilderFMIndex builder = new GraphBuilderFMIndex();
			builder.setConfig(config);
			return builder.buildSimplifiedAssemblyGraph(sequences);
		});

		timeIt("Save the Graph", () -> {
			try {
				assemblyGraph.save(args[1]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		});

		System.out.println("---------Graph Properties--------------");
		assemblyGraph.printInfo();
		System.out.println("---------------------------------------");

	}
}
