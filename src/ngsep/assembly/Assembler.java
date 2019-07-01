package ngsep.assembly;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.sequences.io.FastqFileReader;

public class Assembler {
	private static final String[] fastq = { ".fastq", ".fastq.gz" };
	private static final String[] fasta = { ".fasta", ".fa" };

	private List<CharSequence> sequences;
	private AssemblyGraph graph;

	public Assembler(String fileIn, String fileOut) throws Exception {
		System.out.println("-----Assembler-----");
		sequences = load(fileIn);

		System.out.println("building overlap Graph");
		long ini = System.currentTimeMillis();
		GraphBuilder builder = new GraphBuilderFMIndex();
		graph = builder.buildAssemblyGraph(sequences);
		System.out.println("build overlap Graph: " + (System.currentTimeMillis() - ini) / (double) 1000 + " s");

		System.out.println("building layouts");
		ini = System.currentTimeMillis();
		LayourBuilder pathsFinder = new LayoutBuilderGreedy();
		pathsFinder.findPaths(graph);
		System.out.println("build layouts: " + (System.currentTimeMillis() - ini) / (double) 1000 + " s");
		
		System.out.println("building consensus");
		ini = System.currentTimeMillis();
		ConsensusBuilder consensus = new ConsensusBuilderBidirectionalFMIndex(2, 20, 1, 18, 10, 0.07, 0.03, 1);
		List<CharSequence> AssembleSequences = consensus.makeConsensus(graph);
		System.out.println("build consensus: " + (System.currentTimeMillis() - ini) / (double) 1000 + " s");

		exportToFile(fileOut,"assembled", AssembleSequences);
	}

	/**
	 * Load the sequences of the file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	public static List<CharSequence> load(String filename) throws IOException {
		System.out.println("loading file");
		long ini = System.currentTimeMillis();
		if (Stream.of(fastq)
				.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase()))) {
			System.out.println("load file" + ((System.currentTimeMillis() - ini) / (double) 1000) + " s");
			return loadFastq(filename);
		} else if (Stream.of(fasta)
				.anyMatch((String s) -> filename.endsWith(s.toLowerCase()) || filename.endsWith(s.toUpperCase()))) {
			System.out.println("load file: " + ((System.currentTimeMillis() - ini) / (double) 1000) + " s");
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

	private static void exportToFile(String fileName,String name, List<CharSequence> sequences) throws FileNotFoundException {
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
		new Assembler(args[0], args[1]);
	}
}
