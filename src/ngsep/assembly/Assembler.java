package ngsep.assembly;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.sequences.io.FastqFileReader;
import ngsep.assembly.SimplifiedAssemblyGraph;

import static ngsep.assembly.TimeUtilities.timeGroup;
import static ngsep.assembly.TimeUtilities.timeIt;

public class Assembler {
	private static final String[] fastq = { ".fastq", ".fastq.gz" };
	private static final String[] fasta = { ".fasta", ".fa" };

	private static enum Option {
		Normal, withGraph
	}

	private List<CharSequence> sequences;
	private AssemblyGraph graph;

	public Assembler(String fileIn, String fileOut) throws Exception {
		this(fileIn, fileOut, Option.Normal, new AssemblyConfiguration());
	}

	public Assembler(String fileIn, String fileOut, Option option, AssemblyConfiguration config) throws Exception {
		timeGroup("----------Assembly---------", () -> {
			switch (option) {
			case Normal:
				sequences = timeIt("  Load the sequences", () -> load(fileIn));
				graph = timeGroup("  Build overlap Graph",
						() -> (new GraphBuilderFMIndex()).buildAssemblyGraph(sequences, config).getAssemblyGraph());

				break;

			case withGraph:
				graph = timeIt("  Load the graph", () -> {
					SimplifiedAssemblyGraph sag = new SimplifiedAssemblyGraph(fileIn);
					sag.removeDuplicatedEmbeddes();
					return sag.getAssemblyGraph();
				});
				break;
			}

			timeGroup("  Build layouts", () -> {
				LayourBuilder pathsFinder = new LayoutBuilderGreedy();
				pathsFinder.findPaths(graph);
			});

			List<CharSequence> AssembleSequences = timeGroup("  Build consensus", () -> {
				ConsensusBuilder consensus = new ConsensusBuilderBidirectionalSimple(2, 20, 1, 8, 10, 0.07, 0.00, 1);
				return consensus.makeConsensus(graph);
			});

			timeIt("  Export", () -> {
				try {
					exportToFile(fileOut, "assembled", AssembleSequences);
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
			});
		});
	}

	/**
	 * Load the sequences of the file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	public static List<CharSequence> load(String filename) throws IOException {
		if (Stream.of(fastq)
				.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase()))) {
			return loadFastq(filename);
		} else if (Stream.of(fasta)
				.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase()))) {
			return loadFasta(filename);
		} else
			throw new IOException("the file not is a fasta or fastq file: " + filename);

	}

	/**
	 * Load the sequences of the Fasta file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private static List<CharSequence> loadFasta(String filename) throws IOException {
		List<CharSequence> sequences = new ArrayList<>();
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList seqsQl = handler.loadSequences(filename);
		for (QualifiedSequence seq : seqsQl) {
			DNAMaskedSequence characters = (DNAMaskedSequence) seq.getCharacters();
			sequences.add(characters);
		}
		return sequences;
	}

	/**
	 * Load the sequences of the Fastq file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private static List<CharSequence> loadFastq(String filename) throws IOException {
		List<CharSequence> sequences = new ArrayList<>();
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			reader.setLoadMode(FastqFileReader.LOAD_MODE_MINIMAL);
			reader.setSequenceType(DNAMaskedSequence.class);
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				DNAMaskedSequence characters = (DNAMaskedSequence) read.getCharacters();
				sequences.add(characters);
			}
		}
		return sequences;
	}

	public static void exportToFile(String fileName, String name, Iterable<? extends CharSequence> sequences)
			throws FileNotFoundException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> list = new ArrayList<QualifiedSequence>();
		int i = 1;
		for (CharSequence str : sequences)
			list.add(new QualifiedSequence(name + "_" + (i++), str));
		try (PrintStream pr = new PrintStream(new FileOutputStream(fileName))) {
			handler.saveSequences(list, pr, 1000);
		}
	}

	public static <T extends CharSequence> void exportToFile(String fileName, String name, T[] sequences)
			throws FileNotFoundException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> list = new ArrayList<QualifiedSequence>();
		int i = 1;
		for (CharSequence str : sequences)
			list.add(new QualifiedSequence(name + "_" + (i++), str));
		try (PrintStream pr = new PrintStream(new FileOutputStream(fileName))) {
			handler.saveSequences(list, pr, 1000);
		}
	}

	public static void main(String[] args) throws Exception {
		try {
			Option option = (args.length > 2) ? Option.valueOf(args[2].trim()) : Option.Normal;
			AssemblyConfiguration config = (args.length > 4)
					? new AssemblyConfiguration(Double.valueOf(args[3]), Double.valueOf(args[4]))
					: new AssemblyConfiguration();
			new Assembler(args[0], args[1], option, config);
		} catch (IllegalArgumentException e) {
			System.out.println(
					"Invalid option: '" + args[2] + "' the valid options are: " + Arrays.toString(Option.values()));
		}
	}
}
