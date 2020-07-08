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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmersExtractor;
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

	// Constants for default values
	public static final byte INPUT_FORMAT_FASTQ=KmersExtractor.INPUT_FORMAT_FASTQ;
	public static final byte INPUT_FORMAT_FASTA=KmersExtractor.INPUT_FORMAT_FASTA;
	public static final int DEF_KMER_LENGTH = KmersExtractor.DEF_KMER_LENGTH;
	public static final int DEF_KMER_OFFSET = 15;
	public static final int DEF_MIN_KMER_PCT = KmerHitsAssemblyEdgesFinder.DEF_MIN_KMER_PCT;
	public static final int DEF_NUM_THREADS = GraphBuilderMinimizers.DEF_NUM_THREADS;
	public static final String GRAPH_CONSTRUCTION_ALGORITHM_MINIMIZERS="Minimizers";
	public static final String GRAPH_CONSTRUCTION_ALGORITHM_FMINDEX="FMIndex";
	public static final String LAYOUT_ALGORITHM_MAX_OVERLAP="MaxOverlap";
	public static final String LAYOUT_ALGORITHM_KRUSKAL_PATH="KruskalPath";
	public static final String CONSENSUS_ALGORITHM_SIMPLE="Simple";
	public static final String CONSENSUS_ALGORITHM_POLISHING="Polishing";

	// Logging and progress
	private Logger log = Logger.getLogger(Assembler.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	// Parameters
	private String inputFile = null;
	private String outputPrefix = null;
	private int kmerLength = DEF_KMER_LENGTH;
	private int kmerOffset = DEF_KMER_OFFSET;
	private int minKmerPercentage = DEF_MIN_KMER_PCT;
	private byte inputFormat = INPUT_FORMAT_FASTQ;
	private String graphFile = null;
	private String graphConstructionAlgorithm=GRAPH_CONSTRUCTION_ALGORITHM_MINIMIZERS;
	private String layoutAlgorithm=LAYOUT_ALGORITHM_KRUSKAL_PATH;
	private String consensusAlgorithm=CONSENSUS_ALGORITHM_SIMPLE;
	private int numThreads = DEF_NUM_THREADS;
	
	// Get and set methods
	public Logger getLog() {
		return log;
	}
	public void setLog(Logger log) {
		this.log = log;
	}
	
	public ProgressNotifier getProgressNotifier() {
		return progressNotifier;
	}
	public void setProgressNotifier(ProgressNotifier progressNotifier) { 
		this.progressNotifier = progressNotifier;
	}
	
	public String getInputFile() {
		return inputFile;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	public String getOutputPrefix() {
		return outputPrefix;
	}
	public void setOutputPrefix(String outputPrefix) {
		this.outputPrefix = outputPrefix;
	}
	
	public int getKmerLength() {
		return kmerLength;
	}
	public void setKmerLength(int kmerLength) {
		if(kmerLength<=0) throw new IllegalArgumentException("Kmer length should be a positive number");
		this.kmerLength = kmerLength;
	}
	public void setKmerLength(String value) {
		setKmerLength((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getKmerOffset() {
		return kmerOffset;
	}
	public void setKmerOffset(int kmerOffset) {
		if(kmerOffset<=0) throw new IllegalArgumentException("Kmer offset should be a positive number");
		this.kmerOffset = kmerOffset;
	}
	public void setKmerOffset(String value) {
		setKmerOffset((int)OptionValuesDecoder.decode(value, Integer.class));
	}
	public int getMinKmerPercentage() {
		return minKmerPercentage;
	}
	public void setMinKmerPercentage(int minKmerPercentage) {
		if(minKmerPercentage<0) throw new IllegalArgumentException("Minimum kmer percentage should be a non-negative number");
		if(minKmerPercentage>100) throw new IllegalArgumentException("Minimum kmer percentage should be a number from 0 to 100");
		this.minKmerPercentage = minKmerPercentage;
	}
	public void setMinKmerPercentage(String value) {
		setMinKmerPercentage((int)OptionValuesDecoder.decode(value, Integer.class));
	}

	public byte getInputFormat() {
		return inputFormat;
	}
	public void setInputFormat(byte inputFormat) {
		if (inputFormat!=INPUT_FORMAT_FASTA && inputFormat != INPUT_FORMAT_FASTQ) {
			throw new IllegalArgumentException("Invalid input format "+inputFormat);
		}
		this.inputFormat = inputFormat;
	}
	public void setInputFormat(String value) {
		this.setInputFormat((byte) OptionValuesDecoder.decode(value, Byte.class));
	}
	
	public String getGraphFile() {
		return graphFile;
	}
	public void setGraphFile(String graphFile) {
		this.graphFile = graphFile;
	}
	
	public String getGraphConstructionAlgorithm() {
		return graphConstructionAlgorithm;
	}
	public void setGraphConstructionAlgorithm(String graphConstructionAlgorithm) {
		if(!GRAPH_CONSTRUCTION_ALGORITHM_FMINDEX.equals(graphConstructionAlgorithm) && !GRAPH_CONSTRUCTION_ALGORITHM_MINIMIZERS.equals(graphConstructionAlgorithm)) throw new IllegalArgumentException("Unrecognized graph construction algorithm "+graphConstructionAlgorithm);
		this.graphConstructionAlgorithm = graphConstructionAlgorithm;
	}
	
	public String getLayoutAlgorithm() {
		return layoutAlgorithm;
	}
	public void setLayoutAlgorithm(String layoutAlgorithm) {
		if(!LAYOUT_ALGORITHM_KRUSKAL_PATH.equals(layoutAlgorithm) && !LAYOUT_ALGORITHM_MAX_OVERLAP.equals(layoutAlgorithm)) throw new IllegalArgumentException("Unrecognized layout algorithm "+layoutAlgorithm);
		this.layoutAlgorithm = layoutAlgorithm;
	}
	public String getConsensusAlgorithm() {
		return consensusAlgorithm;
	}
	public void setConsensusAlgorithm(String consensusAlgorithm) {
		if(!CONSENSUS_ALGORITHM_SIMPLE.equals(consensusAlgorithm) && !CONSENSUS_ALGORITHM_POLISHING.equals(consensusAlgorithm)) throw new IllegalArgumentException("Unrecognized consensus algorithm "+consensusAlgorithm);
		this.consensusAlgorithm = consensusAlgorithm;
	}
	public int getNumThreads() {
		return numThreads;
	}
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	public void setNumThreads(String value) {
		this.setNumThreads((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public static void main(String[] args) throws Exception {
		Assembler instance = new Assembler ();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();
	}

	public void run() throws IOException {
		logParameters();
		if(inputFile==null) throw new IOException("The input file with raw reads is required");
		if(outputPrefix==null) throw new IOException("An output prefix is required");
		run (inputFile, outputPrefix);
		log.info("Process finished");
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Prefix for the output files:"+ outputPrefix);
		if (graphFile!=null) out.println("Load assembly graph from: "+graphFile);
		else out.println("Algorithm to build graph: "+graphConstructionAlgorithm);
		out.println("Algorithm to build layout: "+layoutAlgorithm);
		out.println("Algorithm to build consensus: "+consensusAlgorithm);
		out.println("K-mer length: "+ kmerLength);
		//out.println("K-mer offset for FM-index: "+ kmerOffset);
		if (inputFormat == INPUT_FORMAT_FASTQ)  out.println("Fastq format");
		if (inputFormat == INPUT_FORMAT_FASTA)  out.println("Fasta format");
		log.info(os.toString());
	}

	public void run(String inputFile, String outputPrefix) throws IOException {
		List<QualifiedSequence> sequences = load(inputFile,inputFormat);
		log.info("Loaded "+sequences.size()+" sequences");
		if(progressNotifier!=null && !progressNotifier.keepRunning(10)) return;
		AssemblyGraph graph;
		if(graphFile!=null) {
			graph = AssemblyGraph.load(sequences, graphFile);
		} else if (GRAPH_CONSTRUCTION_ALGORITHM_FMINDEX.equals(graphConstructionAlgorithm)) {
			GraphBuilderFMIndex gbIndex = new GraphBuilderFMIndex(kmerLength, kmerOffset, minKmerPercentage, numThreads);
			gbIndex.setLog(log);
			graph =  gbIndex.buildAssemblyGraph(sequences);
		} else {
			GraphBuilderMinimizers builder = new GraphBuilderMinimizers();
			builder.setKmerLength(kmerLength);
			//builder.setMinKmerPercentage(minKmerPercentage);
			builder.setNumThreads(numThreads);
			builder.setLog(log);
			graph = builder.buildAssemblyGraph(sequences);
		}
		log.info("Built assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges");
		graph.updateVertexDegrees();
		log.info("Built assembly graph");
		
		if(progressNotifier!=null && !progressNotifier.keepRunning(50)) return;
		if(graphFile==null) {
			String outFileGraph = outputPrefix+".graph.gz";
			graph.save(outFileGraph);
			log.info("Saved graph to "+outFileGraph);
		}
		graph.removeEdgesChimericReads();
		graph.filterEdgesAndEmbedded();
		log.info("Filtered graph. Vertices: "+graph.getVertices().size()+" edges: "+graph.getEdges().size());
		LayoutBuilder pathsFinder;
		if(LAYOUT_ALGORITHM_MAX_OVERLAP.equals(layoutAlgorithm)) {
			pathsFinder = new LayoutBuilderGreedyMaxOverlap();
			//LayoutBuilder pathsFinder = new LayoutBuilderGreedyMinCost();
		} else {
			pathsFinder= new LayoutBuilderSkeletonBestReciprocal();
			//LayourBuilder pathsFinder = new LayoutBuilderMetricMSTChristofides();
			//LayourBuilder pathsFinder = new LayoutBuilderModifiedKruskal();
		}
		pathsFinder.findPaths(graph);
		log.info("Layout complete. Paths: "+graph.getPaths().size());
		if(progressNotifier!=null && !progressNotifier.keepRunning(60)) return;
		ConsensusBuilder consensus;
		if(CONSENSUS_ALGORITHM_POLISHING.equals(consensusAlgorithm)) {
			consensus = new ConsensusBuilderBidirectionalWithPolishing();
		} else {
			consensus = new ConsensusBuilderBidirectionalSimple();
		}
		List<CharSequence> assembledSequences =  consensus.makeConsensus(graph);
		log.info("Built consensus");
		if(progressNotifier!=null && !progressNotifier.keepRunning(95)) return;
		saveAssembly(outputPrefix+".fa", "contig", assembledSequences);
	}

	/**
	 * Load the sequences of the file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	public static List<QualifiedSequence> load(String filename, byte inputFormat) throws IOException {
		List<QualifiedSequence> sequences;
		if (INPUT_FORMAT_FASTQ == inputFormat) sequences = loadFastq(filename);
		else if (INPUT_FORMAT_FASTA==inputFormat) sequences = loadFasta(filename);
		else throw new IOException("the file not is a fasta or fastq file: " + filename);
		Collections.sort(sequences, (l1, l2) -> l2.getLength() - l1.getLength());
		return sequences;
	}

	/**
	 * Load the sequences of the Fasta file
	 * @param filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private static List<QualifiedSequence> loadFasta(String filename) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList seqsQL = handler.loadSequences(filename);
		List<QualifiedSequence> answer = new ArrayList<>();
		for(QualifiedSequence seq:seqsQL) {
			QualifiedSequence upperCase = new QualifiedSequence(seq.getName(),new DNAMaskedSequence(seq.getCharacters().toString().toUpperCase()));
			answer.add(upperCase);
		}
		return answer;
	}

	/**
	 * Load the sequences of the Fastq file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private static List<QualifiedSequence> loadFastq(String filename) throws IOException {
		List<QualifiedSequence> sequences = new ArrayList<>();
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			reader.setSequenceType(DNAMaskedSequence.class);
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				DNAMaskedSequence characters = (DNAMaskedSequence) read.getCharacters();
				sequences.add(new QualifiedSequence(read.getName(), characters));
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
