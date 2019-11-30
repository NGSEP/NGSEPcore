/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.assembly;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.sequences.io.FastqFileReader;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;

/**
 * @author Jorge Duitama
 * @author Juan Camilo Bojaca
 * @author David Guevara
 */
public class Assembler {

	private Logger log = Logger.getLogger(Assembler.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	
	public static final String INPUT_FORMAT_FASTQ="fastq";
	public static final String INPUT_FORMAT_FASTA="fasta";
	public static final String INPUT_FORMAT_GRAPH="graph";
	public static final int DEF_KMER_SIZE = 15;
	public static final int DEF_KMER_OFFSET = 15;
	public static final int DEF_MIN_KMER_PCT = 40;
	
	private String inputFormat = INPUT_FORMAT_FASTQ;
	
	private String outFileGraph = null;
	private int kmerLength = DEF_KMER_SIZE;
	private int kmerOffset = DEF_KMER_OFFSET;
	private int minKmerPercentage = DEF_MIN_KMER_PCT;
	
	
	public static void main(String[] args) throws Exception {
		Assembler instance = new Assembler ();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		String inputFile = args[i++];
		String outputFile = args[i++];
		instance.run(inputFile, outputFile);
	}
	
	public void setProgressNotifier(ProgressNotifier progressNotifier) { 
		this.progressNotifier = progressNotifier;
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	
	/**
	 * @return the log
	 */
	public Logger getLog() {
		return log;
	}
	
	/**
	 * @param log the log to set
	 */
	public void setLog(Logger log) {
		this.log = log;
	}

	/**
	 * @return the inputFormat
	 */
	public String getInputFormat() {
		return inputFormat;
	}

	/**
	 * @param inputFormat the inputFormat to set
	 */
	public void setInputFormat(String inputFormat) {
		this.inputFormat = inputFormat;
	}
	
	/**
	 * @return the outFileGraph
	 */
	public String getOutFileGraph() {
		return outFileGraph;
	}

	/**
	 * @param outFileGraph the outFileGraph to set
	 */
	public void setOutFileGraph(String outFileGraph) {
		this.outFileGraph = outFileGraph;
	}

	/**
	 * @return the kmerLength
	 */
	public int getKmerLength() {
		return kmerLength;
	}

	/**
	 * @param kmerLength the kmerLength to set
	 */
	public void setKmerLength(int kmerLength) {
		this.kmerLength = kmerLength;
	}

	public void setKmerLength(String value) {
		setKmerLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	/**
	 * @return the kmerOffset
	 */
	public int getKmerOffset() {
		return kmerOffset;
	}

	/**
	 * @param kmerOffset the kmerOffset to set
	 */
	public void setKmerOffset(int kmerOffset) {
		this.kmerOffset = kmerOffset;
	}

	public void setKmerOffset(String value) {
		setKmerOffset((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	/**
	 * @return the minKmerPercentage
	 */
	public int getMinKmerPercentage() {
		return minKmerPercentage;
	}

	/**
	 * @param minKmerPercentage the minKmerPercentage to set
	 */
	public void setMinKmerPercentage(int minKmerPercentage) {
		this.minKmerPercentage = minKmerPercentage;
	}
	
	public void setMinKmerPercentage(String value) {
		setMinKmerPercentage((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	public void run(String inputFile, String outputFile) throws IOException {
		AssemblyGraph graph;
		if(INPUT_FORMAT_GRAPH.equals(inputFormat)) {
			graph = AssemblyGraph.load(inputFile);
			log.info("Loaded assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges");
		} else {
			List<CharSequence> sequences = load(inputFile);
			log.info("Loaded "+sequences.size()+" sequences");
			Collections.sort(sequences, (l1, l2) -> l2.length() - l1.length());
			log.info("Sorted "+sequences.size()+" sequences");
			List<CharSequence> finalSequences = Collections.unmodifiableList(sequences);
			GraphBuilderFMIndex gbIndex = new GraphBuilderFMIndex(kmerLength, kmerOffset, minKmerPercentage);
			gbIndex.setLog(log);
			graph =  gbIndex.buildAssemblyGraph(finalSequences);
			log.info("Built graph");
			if(outFileGraph!=null) {
				graph.serialize(outFileGraph);
				log.info("Saved graph in "+outFileGraph);
			}
		}

		LayourBuilder pathsFinder = new LayoutBuilderGreedyMinCost();
		pathsFinder.findPaths(graph);
		log.info("Layout complete. Paths: "+graph.getPaths().size());

		ConsensusBuilder consensus = new ConsensusBuilderBidirectionalSimple();
		List<CharSequence> assembledSequences =  consensus.makeConsensus(graph);
		log.info("Built consensus");
		saveAssembly(outputFile, "contig", assembledSequences);
	}

	/**
	 * Load the sequences of the file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	public List<CharSequence> load(String filename) throws IOException {
		if (INPUT_FORMAT_FASTQ.equals(inputFormat)) return loadFastq(filename);
		else if (INPUT_FORMAT_FASTA.equals(inputFormat)) return loadFasta(filename);
		else throw new IOException("the file not is a fasta or fastq file: " + filename);
	}

	/**
	 * Load the sequences of the Fasta file
	 * @param filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private List<CharSequence> loadFasta(String filename) throws IOException {
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
	private List<CharSequence> loadFastq(String filename) throws IOException {
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
	/**
	 * Saves the given sequences in fasta format
	 * @param filename name of the output file
	 * @param prefix of the sequence names
	 * @param sequences List of sequences corresponding to the final assembly
	 * @throws IOException If the file can not be generated
	 */
	public void saveAssembly(String filename, String prefix, List<CharSequence> sequences) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> list = new ArrayList<QualifiedSequence>();
		int i = 1;
		for (CharSequence str : sequences) {
			list.add(new QualifiedSequence(prefix + "_" + i, str));
			i++;
		}
			
		try (PrintStream out = new PrintStream(filename)) {
			handler.saveSequences(list, out, 100);
		}
	}

	
}
