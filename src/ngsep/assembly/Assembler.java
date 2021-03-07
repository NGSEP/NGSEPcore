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
import java.util.Set;
import java.util.logging.Logger;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastaSequencesHandler;
import ngsep.sequences.io.FastqFileReader;
import ngsep.assembly.io.AssemblyGraphFileHandler;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.math.Distribution;

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
	public static final int DEF_WINDOW_LENGTH = 30;
	public static final int DEF_MIN_READ_LENGTH = 5000;
	public static final int DEF_PLOIDY = AssemblyGraph.DEF_PLOIDY_ASSEMBLY;
	public static final int DEF_BP_HOMOPOLYMER_COMPRESSION = 0;
	public static final double DEF_MIN_SCORE_PROPORTION_EDGES = 0.3;
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
	private int windowLength = DEF_WINDOW_LENGTH;
	private int minReadLength = DEF_MIN_READ_LENGTH;
	private byte inputFormat = INPUT_FORMAT_FASTQ;
	private String graphFile = null;
	private String graphConstructionAlgorithm=GRAPH_CONSTRUCTION_ALGORITHM_MINIMIZERS;
	private String layoutAlgorithm=LAYOUT_ALGORITHM_KRUSKAL_PATH;
	private String consensusAlgorithm=CONSENSUS_ALGORITHM_POLISHING;
	private boolean correctReads = false;
	private int ploidy = DEF_PLOIDY;
	private int bpHomopolymerCompression = DEF_BP_HOMOPOLYMER_COMPRESSION;
	private double minScoreProportionEdges = DEF_MIN_SCORE_PROPORTION_EDGES;
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
	
	public int getWindowLength() {
		return windowLength;
	}
	public void setWindowLength(int windowLength) {
		if(windowLength<=0) throw new IllegalArgumentException("Window length should be a positive number");
		this.windowLength = windowLength;
	}
	public void setWindowLength(String value) {
		setWindowLength((int)OptionValuesDecoder.decode(value, Integer.class));
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
	
	public int getPloidy() {
		return ploidy;
	}
	public void setPloidy(int ploidy) {
		this.ploidy = ploidy;
	}
	public void setPloidy(String value) {
		this.setPloidy((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public int getBpHomopolymerCompression() {
		return bpHomopolymerCompression;
	}
	public void setBpHomopolymerCompression(int bpHomopolymerCompression) {
		this.bpHomopolymerCompression = bpHomopolymerCompression;
	}
	public void setBpHomopolymerCompression(String value) {
		this.setBpHomopolymerCompression((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public double getMinScoreProportionEdges() {
		return minScoreProportionEdges;
	}
	public void setMinScoreProportionEdges(double minScoreProportionEdges) {
		this.minScoreProportionEdges = minScoreProportionEdges;
	}
	public void setMinScoreProportionEdges(String value) {
		this.setMinScoreProportionEdges((double) OptionValuesDecoder.decode(value, Double.class));
	}
	
	public boolean isCorrectReads() {
		return correctReads;
	}
	public void setCorrectReads(boolean correctReads) {
		this.correctReads = correctReads;
	}
	public void setCorrectReads(Boolean correctReads) {
		this.setCorrectReads(correctReads.booleanValue());
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
		out.println("Window length for minimizers: "+windowLength);
		if(bpHomopolymerCompression>0) out.println("Run homopolymer compression keeping at most "+bpHomopolymerCompression+" consecutive base pairs");
		out.println("Minimum score proportion (from the maximum score) to keep edges of a sequence: "+ minScoreProportionEdges);
		//out.println("K-mer offset for FM-index: "+ kmerOffset);
		if (inputFormat == INPUT_FORMAT_FASTQ)  out.println("Fastq format");
		if (inputFormat == INPUT_FORMAT_FASTA)  out.println("Fasta format");
		out.println("Number of threads "+numThreads);
		log.info(os.toString());
	}

	public void run(String inputFile, String outputPrefix) throws IOException {
		Runtime runtime = Runtime.getRuntime();
		long startTime = System.currentTimeMillis();
		List<QualifiedSequence> sequences = load(inputFile,inputFormat, minReadLength);
		long totalBp = 0;
		Distribution distReadLength = new Distribution(0, 100000, 1000);
		for(QualifiedSequence seq:sequences) {
			distReadLength.processDatapoint(seq.getLength());
			totalBp+=seq.getLength();
		}
		log.info("Loaded "+sequences.size()+" sequences. Total basepairs: "+totalBp);
		distReadLength.printDistributionInt(System.out);
		if(progressNotifier!=null && !progressNotifier.keepRunning(10)) return;
		AssemblyGraph graph;
		if(graphFile!=null) {
			graph = AssemblyGraphFileHandler.load(sequences, graphFile);
			log.info("Loaded assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges");
		} else {
			double [] compressionFactors =null;
			if (bpHomopolymerCompression>0) {
				compressionFactors = runHomopolymerCompression (sequences);
				log.info("Performed homopolymer compression");
			}
			
			GraphBuilderMinimizers builder = new GraphBuilderMinimizers();
			builder.setKmerLength(kmerLength);
			builder.setWindowLength(windowLength);
			builder.setPloidy(ploidy);
			builder.setNumThreads(numThreads);
			builder.setLog(log);
			graph = builder.buildAssemblyGraph(sequences,compressionFactors);
			if(bpHomopolymerCompression>0) {
				List<QualifiedSequence> originalSeqs = load(inputFile,inputFormat, minReadLength);
				log.info("Loaded original sequences to restore. Compressed sequences: "+sequences.size()+". Loaded: "+originalSeqs.size());
				for(int i=0;i<sequences.size();i++) {
					QualifiedSequence seq = sequences.get(i);
					CharSequence original = originalSeqs.get(i).getCharacters();
					seq.setCharacters(original);
				}
			}
			log.info("Built assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges");
		}
		
		if(progressNotifier!=null && !progressNotifier.keepRunning(50)) return;
		if(graphFile==null) {
			String outFileGraph = outputPrefix+".graph.gz";
			AssemblyGraphFileHandler.save(graph, outFileGraph);
			log.info("Saved graph to "+outFileGraph);
		}
		
		LayoutBuilder pathsFinder;
		if(LAYOUT_ALGORITHM_MAX_OVERLAP.equals(layoutAlgorithm)) {
			pathsFinder = new LayoutBuilderGreedyMaxOverlap();
			//LayoutBuilder pathsFinder = new LayoutBuilderGreedyMinCost();
		} else {
			pathsFinder= new LayoutBuilderKruskalPath();
			//LayourBuilder pathsFinder = new LayoutBuilderMetricMSTChristofides();
			//LayourBuilder pathsFinder = new LayoutBuilderModifiedKruskal();
		}
		ConsensusBuilder consensus;
		if(CONSENSUS_ALGORITHM_POLISHING.equals(consensusAlgorithm)) {
			ConsensusBuilderBidirectionalWithPolishing consensusP = new ConsensusBuilderBidirectionalWithPolishing();
			consensusP.setNumThreads(numThreads);
			if(correctReads) consensusP.setCorrectedReadsFile(outputPrefix+"_correctedReads.fa.gz");
			consensus = consensusP;
		} else {
			consensus = new ConsensusBuilderBidirectionalSimple();
		}
		
		graph.removeVerticesChimericReads();
		log.info("Filtered chimeric reads. Vertices: "+graph.getVertices().size()+" edges: "+graph.getEdges().size());
		//graph.updateScores(true);
		long time2 = System.currentTimeMillis();
		AssemblySequencesRelationshipFilter filter = new AssemblySequencesRelationshipFilter();
		List<QualifiedSequence> assembledSequences = new ArrayList<QualifiedSequence>();
		if(ploidy > 1) {
			AssemblyGraph diploidGraph = graph.buildSubgraph(null);
			log.info("Copied graph. New graph has "+diploidGraph.getVertices().size()+" vertices and "+diploidGraph.getEdges().size()+" edges");
			filter.filterEdgesAndEmbedded(diploidGraph, minScoreProportionEdges);
			diploidGraph.updateScores(true);
			log.info("Filtered graph. New graph has now "+diploidGraph.getVertices().size()+" vertices and "+diploidGraph.getEdges().size()+" edges");
			if (pathsFinder instanceof LayoutBuilderKruskalPath) ((LayoutBuilderKruskalPath)pathsFinder).setMinPathLength(0);
			pathsFinder.findPaths(diploidGraph);
			log.info("Building haplotype subgraphs");
			HaplotypeReadsClusterCalculator hapsCalculator = new HaplotypeReadsClusterCalculator();
			hapsCalculator.setLog(log);
			hapsCalculator.setNumThreads(numThreads);
			List<Set<Integer>> readIdsClusters =  hapsCalculator.clusterReads(diploidGraph, ploidy);
			log.info("Separated reads in "+readIdsClusters.size()+" clusters");
			if (pathsFinder instanceof LayoutBuilderKruskalPath) ((LayoutBuilderKruskalPath)pathsFinder).setMinPathLength(5);
			int haplotypeNumber= 0;
			for(Set<Integer> readIdsCluster: readIdsClusters) {
				AssemblyGraph haplotypeGraph = graph.buildSubgraph(readIdsCluster);
				log.info("Built haplotype subgraph with "+haplotypeGraph.getVertices().size()+" vertices and "+haplotypeGraph.getNumEdges()+ " edges from "+readIdsCluster.size()+" reads");
				String outFileGraph = outputPrefix+"_hap"+haplotypeNumber+".graph.gz";
				AssemblyGraphFileHandler.save(haplotypeGraph, outFileGraph);
				log.info("Saved graph to "+outFileGraph);
				filter.filterEdgesAndEmbedded(haplotypeGraph, minScoreProportionEdges);
				haplotypeGraph.updateScores(true);
				pathsFinder.findPaths(haplotypeGraph);
				log.info("Built "+haplotypeGraph.getPaths().size()+" paths for next haplotype cluster with "+readIdsCluster.size()+" reads");
				consensus.setSequenceNamePrefix("ContigHap"+haplotypeNumber);
				List<QualifiedSequence> sequencesCluster = consensus.makeConsensus(haplotypeGraph);
				log.info("Assembled "+sequencesCluster.size()+" sequences for next haplotype cluster");
				assembledSequences.addAll(sequencesCluster);
				haplotypeNumber++;
			}
		} else {
			filter.filterEdgesAndEmbedded(graph, minScoreProportionEdges);
			graph.updateScores(true);
			
			pathsFinder.findPaths(graph);
			if(progressNotifier!=null && !progressNotifier.keepRunning(60)) return;
			assembledSequences.addAll(consensus.makeConsensus(graph));
			
		}
		
		FastaSequencesHandler handler = new FastaSequencesHandler();
		try (PrintStream out = new PrintStream(outputPrefix+"_initial.fa")) {
			handler.saveSequences(assembledSequences, out, 100);
		}
		long usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time3 = System.currentTimeMillis();
		long diff1 = (time3-time2)/1000;
		long diff2 = (time3-startTime)/1000;
		log.info("Layout and initial consensus complete. Memory: "+usedMemory+" Time process (s): "+diff1+" total time (s): "+diff2);
		List<Integer> lengths = new ArrayList<Integer>();
		for(QualifiedSequence seq:assembledSequences) lengths.add(seq.getLength());
		long [] nStats = NStatisticsCalculator.calculateNStatistics(lengths);
		System.out.println("Initial consensus N statistics");
		NStatisticsCalculator.printNStatistics(nStats, System.out);
		if(progressNotifier!=null && !progressNotifier.keepRunning(95)) return;
		log.info("Built initial consensus. Merging contig ends");
		ContigEndsMerger merger = new ContigEndsMerger();
		assembledSequences = merger.mergeContigs(assembledSequences);
		try (PrintStream out = new PrintStream(outputPrefix+".fa")) {
			handler.saveSequences(assembledSequences, out, 100);
		}
		
		lengths.clear();
		for(QualifiedSequence seq:assembledSequences) lengths.add(seq.getLength());
		nStats = NStatisticsCalculator.calculateNStatistics(lengths);
		System.out.println("Final assembly N statistics");
		NStatisticsCalculator.printNStatistics(nStats, System.out);
		
		
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time4 = System.currentTimeMillis();
		diff1 = (time4-time3)/1000;
		diff2 = (time4-startTime)/1000;
		log.info("Finished consensus. Memory: "+usedMemory+" Time consensus (s): "+diff1+" total time (s): "+diff2);
	}

	private double [] runHomopolymerCompression(List<QualifiedSequence> sequences) {
		double [] compressionFactors = new double[sequences.size()];
		for(int i=0;i<sequences.size();i++) {
			QualifiedSequence seq = sequences.get(i);
			compressionFactors[i] = compressHomopolymers(seq);
		}
		return compressionFactors;
	}
	private double compressHomopolymers(QualifiedSequence seq) {
		String seqStr = seq.getCharacters().toString();
		int n = seqStr.length();
		StringBuilder compressed = new StringBuilder(n);
		char c2 = 0;
		int homopolymerCount = 0;
		for (int i=0;i<n;i++) {
			char c = seqStr.charAt(i);
			if (c==c2) homopolymerCount++;
			else homopolymerCount = 1;
			if(homopolymerCount<=bpHomopolymerCompression) compressed.append(c);
			c2=c;
		}
		double answer = compressed.length();
		if(n>0) answer /=n;
		seq.setCharacters(new DNAMaskedSequence(compressed));
		return answer;
	}
	/**
	 * Load the sequences of the file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	public static List<QualifiedSequence> load(String filename, byte inputFormat, int minReadLength) throws IOException {
		List<QualifiedSequence> sequences;
		if (INPUT_FORMAT_FASTQ == inputFormat) sequences = loadFastq(filename,minReadLength);
		else if (INPUT_FORMAT_FASTA==inputFormat) sequences = loadFasta(filename, minReadLength);
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
	private static List<QualifiedSequence> loadFasta(String filename, int minReadLength) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList seqsQL = handler.loadSequences(filename);
		List<QualifiedSequence> answer = new ArrayList<QualifiedSequence>();
		for(QualifiedSequence seq:seqsQL) {
			if(seq.getLength()>=minReadLength) answer.add(seq);
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
	private static List<QualifiedSequence> loadFastq(String filename, int minReadLength) throws IOException {
		List<QualifiedSequence> sequences = new ArrayList<>();
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			reader.setSequenceType(DNAMaskedSequence.class);
			//TODO: Option to load quality scores
			reader.setLoadMode(FastqFileReader.LOAD_MODE_WITH_NAME);
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				DNAMaskedSequence characters = (DNAMaskedSequence) read.getCharacters();
				if(characters.length()>=minReadLength) sequences.add(new QualifiedSequence(read.getName(), characters));
			}
		}
		return sequences;
	}
}
