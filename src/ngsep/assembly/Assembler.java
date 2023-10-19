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
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.DNASequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMap;
import ngsep.sequences.QualifiedSequence;
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
	public static final int DEF_KMER_LENGTH = 25;
	public static final int DEF_WINDOW_LENGTH = 40;
	public static final int DEF_MIN_READ_LENGTH = 5000;
	public static final int DEF_PLOIDY = AssemblyGraph.DEF_PLOIDY_ASSEMBLY;
	public static final int DEF_BP_HOMOPOLYMER_COMPRESSION = 0;
	public static final int DEF_ERROR_CORRCTION_ROUNDS = 0;
	public static final double DEF_MIN_SCORE_PROPORTION_EDGES = 0.5;
	public static final int DEF_NUM_THREADS = GraphBuilderMinimizers.DEF_NUM_THREADS;
	public static final int DEF_CIRCULAR_MAX_LENGTH = 0;
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
	private int ploidy = DEF_PLOIDY;
	private int errorCorrectionRounds = DEF_ERROR_CORRCTION_ROUNDS;
	private int circularMoleculesMaxLength = DEF_CIRCULAR_MAX_LENGTH;
	private String circularMoleculesStartsFile;
	private int bpHomopolymerCompression = DEF_BP_HOMOPOLYMER_COMPRESSION;
	private double minScoreProportionEdges = DEF_MIN_SCORE_PROPORTION_EDGES;
	private boolean saveCorrected = false;
	private int numThreads = DEF_NUM_THREADS;
	
	//Model objects
	private AssemblySequencesRelationshipFilter relationshipsFilter = new AssemblySequencesRelationshipFilter();
	private LayoutBuilder pathsFinder = createLayoutBuilder();
	private ConsensusBuilder consensus = createConsensusBuilder();
	
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
	
	public int getMinReadLength() {
		return minReadLength;
	}
	public void setMinReadLength(int minReadLength) {
		this.minReadLength = minReadLength;
	}
	public void setMinReadLength(String value) {
		setMinReadLength((int)OptionValuesDecoder.decode(value, Integer.class));
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
		pathsFinder = createLayoutBuilder();	
	}
	public String getConsensusAlgorithm() {
		return consensusAlgorithm;
	}
	public void setConsensusAlgorithm(String consensusAlgorithm) {
		if(!CONSENSUS_ALGORITHM_SIMPLE.equals(consensusAlgorithm) && !CONSENSUS_ALGORITHM_POLISHING.equals(consensusAlgorithm)) throw new IllegalArgumentException("Unrecognized consensus algorithm "+consensusAlgorithm);
		this.consensusAlgorithm = consensusAlgorithm;
		consensus = createConsensusBuilder();
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
	
	public int getCircularMoleculesMaxLength() {
		return circularMoleculesMaxLength;
	}
	public void setCircularMoleculesMaxLength(int circularMoleculesMaxLength) {
		this.circularMoleculesMaxLength = circularMoleculesMaxLength;
	}
	public void setCircularMoleculesMaxLength(String value) {
		this.setCircularMoleculesMaxLength((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public String getCircularMoleculesStartsFile() {
		return circularMoleculesStartsFile;
	}
	public void setCircularMoleculesStartsFile(String circularMoleculesStartsFile) {
		this.circularMoleculesStartsFile = circularMoleculesStartsFile;
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
	
	public int getErrorCorrectionRounds() {
		return errorCorrectionRounds;
	}
	public void setErrorCorrectionRounds(int errorCorrectionRounds) {
		this.errorCorrectionRounds = errorCorrectionRounds;
		if(errorCorrectionRounds>0) saveCorrected = true;
	}
	public void setErrorCorrectionRounds(String value) {
		this.setErrorCorrectionRounds((int) OptionValuesDecoder.decode(value, Integer.class));
	}
	
	public boolean isSaveCorrected() {
		return saveCorrected;
	}
	public void setSaveCorrected(boolean saveCorrected) {
		this.saveCorrected = saveCorrected;
	}
	public void setSaveCorrected(Boolean saveCorrected) {
		this.setSaveCorrected(saveCorrected.booleanValue());
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

	public void run() throws IOException, InterruptedException {
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
		//else out.println("Algorithm to build graph: "+graphConstructionAlgorithm);
		//out.println("Algorithm to build layout: "+layoutAlgorithm);
		out.println("Algorithm to build consensus: "+consensusAlgorithm);
		out.println("Kmer length: "+kmerLength);
		out.println("Window length for minimizers: "+windowLength);
		out.println("Minimum read length: "+minReadLength);
		if(bpHomopolymerCompression>0) out.println("Run homopolymer compression keeping at most "+bpHomopolymerCompression+" consecutive base pairs");
		out.println("Minimum score proportion (from the maximum score) to keep edges of a sequence: "+ minScoreProportionEdges);
		out.println("Sample ploidy: "+ploidy);
		out.println("Maximum length of circular molecules: "+circularMoleculesMaxLength);
		out.println("Fasta file with known start sequences of circular molecules: "+circularMoleculesStartsFile);
		if (inputFormat == INPUT_FORMAT_FASTQ)  out.println("Fastq format");
		if (inputFormat == INPUT_FORMAT_FASTA)  out.println("Fasta format");
		out.println("Number of threads "+numThreads);
		log.info(os.toString());
	}

	public void run(String inputFile, String outputPrefix) throws IOException, InterruptedException {
		Runtime runtime = Runtime.getRuntime();
		long startTime = System.currentTimeMillis();
		KmersExtractor extractor = new KmersExtractor();
		extractor.setLog(log);
		extractor.setNumThreads(numThreads);
		extractor.setReadNCharacters(false);
		
		
		
		List<QualifiedSequence> sequences;
		//correctReads(sequences,map);
		AssemblyGraph graph;
		KmersMap map = null;
		long usedMemory;
		List<QualifiedSequence> compressedSeqs;
		if(graphFile!=null) {
			sequences = load(inputFile, inputFormat, minReadLength);
			compressedSeqs = sequences;
			graph = AssemblyGraphFileHandler.load(sequences, graphFile);
			usedMemory = runtime.totalMemory()-runtime.freeMemory();
			usedMemory/=1000000000;
			log.info("Loaded assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges. Memory: "+usedMemory);
		} else {
			log.info("Calculating kmers distribution");
			extractor.setLoadSequences(true);
			extractor.setInputFormat(inputFormat);
			extractor.setMinReadLength(minReadLength);
			//The conditional avoids creating twice the large array in ShortArrayKmersMapImpl
			//if(extractor.getKmerLength()!=kmerLength) extractor.setKmerLength(kmerLength);
			extractor.processFile(inputFile);
			
			sequences = extractor.getLoadedSequences();
			Collections.sort(sequences, (l1, l2) -> l2.getLength() - l1.getLength());
			map = extractor.getKmersMap();
			long totalBp = 0;
			Distribution distReadLength = new Distribution(0, 100000, 1000);
			for(QualifiedSequence seq:sequences) {
				distReadLength.processDatapoint(seq.getLength());
				totalBp+=seq.getLength();
			}
			log.info("Loaded "+sequences.size()+" sequences. Total basepairs: "+totalBp);
			distReadLength.printDistributionInt(System.out);
			usedMemory = runtime.totalMemory()-runtime.freeMemory();
			usedMemory/=1000000000;
			long time1 = System.currentTimeMillis();
			long diff1 = (time1-startTime)/1000;
			log.info("Reads loaded. Time(s): "+diff1+" Memory (Gbp): "+usedMemory);
			if(progressNotifier!=null && !progressNotifier.keepRunning(10)) return;
			compressedSeqs = sequences;
			if (bpHomopolymerCompression>0) {
				compressedSeqs = runHomopolymerCompression (sequences);
				log.info("Performed homopolymer compression");
				extractor.dispose();
				extractor.processQualifiedSequences(compressedSeqs);
				map = extractor.getKmersMap();
			}
			graph = buildGraph(compressedSeqs, map);
			usedMemory = runtime.totalMemory()-runtime.freeMemory();
			usedMemory/=1000000000;
			long time2 = System.currentTimeMillis();
			long diff = (time2-startTime)/1000;
			log.info("Built assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges. Total Time: "+diff+" Memory: "+usedMemory);
		}
		
		if(progressNotifier!=null && !progressNotifier.keepRunning(50)) return;
		if(graphFile==null && errorCorrectionRounds>0) {
			String outFileGraph = outputPrefix+"_uncorrected.graph.gz";
			AssemblyGraphFileHandler.save(graph, outFileGraph);
			log.info("Saved uncorrected graph to "+outFileGraph);
		}
		//Current error correction. By now not phase sensitive
		AlignmentBasedIndelErrorsCorrector indelCorrector = new AlignmentBasedIndelErrorsCorrector();
		indelCorrector.setNumThreads(numThreads);
		indelCorrector.setLog(log);
		for(int i=0;i<errorCorrectionRounds;i++) {
			long startRound = System.currentTimeMillis();
			usedMemory = runtime.totalMemory()-runtime.freeMemory();
			usedMemory/=1000000000;
			log.info("Started round "+(i+1)+" of error correction. Memory: "+usedMemory);
			indelCorrector.correctErrors(graph,(sequences!=compressedSeqs?sequences:null));
			sequences = new ArrayList<>();
			sequences.addAll(graph.getSequences());
			Collections.sort(sequences,(s1,s2)->s2.getLength()-s1.getLength());
			compressedSeqs = sequences;
			long timeRound = System.currentTimeMillis()-startRound;
			usedMemory = runtime.totalMemory()-runtime.freeMemory();
			usedMemory/=1000000000;
			log.info("Finished error correction process "+(i+1)+". Time: "+(timeRound/1000)+" Memory: "+usedMemory);
			extractor.setMinReadLength(0);
			extractor.dispose();
			extractor.processQualifiedSequences(sequences);
			map = extractor.getKmersMap();
			graph = buildGraph(sequences, map);
		}
		long time2 = System.currentTimeMillis();
		List<QualifiedSequence> assembledSequences = new ArrayList<QualifiedSequence>();
		int value = ploidy; 
		while(value > 1) {
		//while(value > 0) {
			AssemblyGraph copyGraph = graph.buildSubgraph(null);
			log.info("Copied graph. New graph has "+copyGraph.getVertices().size()+" vertices and "+copyGraph.getEdges().size()+" edges");
			copyGraph.removeVerticesChimericReads();
			copyGraph.updateScores(0);
			relationshipsFilter.filterEdgesAndEmbedded(copyGraph, minScoreProportionEdges);
			//diploidGraph.updateScores();
			log.info("Filtered copy graph. New graph has now "+copyGraph.getVertices().size()+" vertices and "+copyGraph.getEdges().size()+" edges");
			//if (pathsFinder instanceof LayoutBuilderKruskalPath) ((LayoutBuilderKruskalPath)pathsFinder).setMinPathLength(0);
			pathsFinder.findPaths(copyGraph);
			log.info("Filtering graph by phasing");
			if (compressedSeqs!=sequences) {
				copyGraph.replaceSequences(sequences);
			}
			HaplotypeReadsClusterCalculator hapsCalculator = new HaplotypeReadsClusterCalculator();
			hapsCalculator.setLog(log);
			hapsCalculator.setNumThreads(numThreads);
			hapsCalculator.setGlobalPloidy(value);
			Map<Integer,ReadPathPhasingData> readsData = hapsCalculator.calculatePathReadsPhasingData(copyGraph, ploidy);
			saveReadsPhasingData (graph, readsData, outputPrefix+"_FV"+value+"_phasedReadsData.txt");
			filterGraphWithPhasingData(graph,readsData);
			String outFileGraph = outputPrefix+"_FV"+value+".graph.gz";
			AssemblyGraphFileHandler.save(graph, outFileGraph);
			log.info("Saved graph with phase filtering to "+outFileGraph);
			value /=2;
		}
		//Save final graph and corrected reads
		if(graphFile==null || errorCorrectionRounds>0) {
			saveGraphAndCorrectedReads(outputPrefix, sequences, graph);
		}
		graph.removeVerticesChimericReads();
		log.info("Filtered chimeric reads. Vertices: "+graph.getVertices().size()+" edges: "+graph.getEdges().size());
		
		graph.updateScores(0);
		relationshipsFilter.filterEdgesAndEmbedded(graph, minScoreProportionEdges);
		//graph.updateScores();
		//if (pathsFinder instanceof LayoutBuilderKruskalPath) ((LayoutBuilderKruskalPath)pathsFinder).setMinPathLength(6);
		pathsFinder.findPaths(graph);
		if(progressNotifier!=null && !progressNotifier.keepRunning(60)) return;
		if (compressedSeqs!=sequences) {
			graph.replaceSequences(sequences);
		}
		assembledSequences.addAll(consensus.makeConsensus(graph));
		
		FastaSequencesHandler handler = new FastaSequencesHandler();
		try (PrintStream out = new PrintStream(outputPrefix+"_initial.fa")) {
			handler.saveSequences(assembledSequences, out, 100);
		}
		usedMemory = runtime.totalMemory()-runtime.freeMemory();
		usedMemory/=1000000000;
		long time3 = System.currentTimeMillis();
		long diff1 = (time3-time2)/1000;
		long diff2 = (time3-startTime)/1000;
		log.info("Layout and initial consensus complete. Time process (s): "+diff1+" total time (s): "+diff2+" Memory: "+usedMemory);
		List<Integer> lengths = new ArrayList<Integer>();
		for(QualifiedSequence seq:assembledSequences) lengths.add(seq.getLength());
		long [] nStats = NStatisticsCalculator.calculateNStatistics(lengths);
		System.out.println("Initial consensus N statistics");
		NStatisticsCalculator.printNStatistics(nStats, System.out);
		if(progressNotifier!=null && !progressNotifier.keepRunning(95)) return;
		//Final merginig
		log.info("Built initial consensus. Merging contig ends");
		ContigEndsMerger merger = new ContigEndsMerger();
		assembledSequences = merger.mergeContigs(assembledSequences);
		//Circularization
		CircularSequencesProcessor circularizator = new CircularSequencesProcessor();
		circularizator.setLog(log);
		circularizator.setMaxLength(circularMoleculesMaxLength);
		circularizator.setStarts(circularMoleculesStartsFile);
		circularizator.processContigs(assembledSequences);
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
	private void saveGraphAndCorrectedReads(String outputPrefix, List<QualifiedSequence> sequences, AssemblyGraph graph)
			throws IOException, FileNotFoundException {
		String outFileGraph = outputPrefix+".graph.gz";
		if(errorCorrectionRounds>0 && saveCorrected) {
			String outFileCorrectedReads = outputPrefix+"_correctedReads.fa.gz";
			try (OutputStream os = new GZIPOutputStream(new FileOutputStream(outFileCorrectedReads));
				 PrintStream outReads = new PrintStream(os)) {
				for(QualifiedSequence seq:sequences) {
					outReads.println(">"+seq.getName());
					outReads.println(seq.getCharacters());
				}
			}
			log.info("Saved corrected reads to "+outFileCorrectedReads);
			outFileGraph = outputPrefix+"_corrected.graph.gz";
		}
		if(errorCorrectionRounds==0 || saveCorrected) {
			AssemblyGraphFileHandler.save(graph, outFileGraph);
			log.info("Saved graph to "+outFileGraph);
		}
	}
	private ConsensusBuilder createConsensusBuilder() {
		ConsensusBuilder consensus;
		if(CONSENSUS_ALGORITHM_POLISHING.equals(consensusAlgorithm)) {
			ConsensusBuilderBidirectionalWithPolishing consensusP = new ConsensusBuilderBidirectionalWithPolishing();
			consensusP.setNumThreads(numThreads);
			consensus = consensusP;
		} else {
			consensus = new ConsensusBuilderBidirectionalSimple();
			consensus.setNumThreads(numThreads);
		}
		return consensus;
	}
	private LayoutBuilder createLayoutBuilder() {
		LayoutBuilder pathsFinder;
		if(LAYOUT_ALGORITHM_MAX_OVERLAP.equals(layoutAlgorithm)) {
			pathsFinder = new LayoutBuilderGreedyMaxOverlap();
			//LayoutBuilder pathsFinder = new LayoutBuilderGreedyMinCost();
		} else {
			pathsFinder= new LayoutBuilderKruskalPath();
			//((LayoutBuilderKruskalPath)pathsFinder).setRunImprovementAlgorithms(false);
		}
		return pathsFinder;
	}
	private void saveReadsPhasingData(AssemblyGraph graph, Map<Integer, ReadPathPhasingData> readsData, String outFilename) throws IOException {
		try (PrintStream out = new PrintStream(outFilename)) {
			for(int i=0;i<graph.getNumSequences();i++) {
				QualifiedSequence seq = graph.getSequence(i);
				out.print(i+"\t"+seq.getName());
				if(graph.isChimeric(i)) {
					out.println("\tChimeric");
					continue;
				}
				ReadPathPhasingData data = readsData.get(i);
				if(data==null) {
					out.println("\tNoData");
					continue;
				}
				out.println(" \t"+data.getPathId()+"\t"+data.getBlockNumber()+"\t"+data.getPhaseWithinBlock()+"\t"+data.isInHomozygousRegion()+"\t"+data.getReadDepth());
			}
		}
		
		
	}
	private void filterGraphWithPhasingData(AssemblyGraph graph, Map<Integer, ReadPathPhasingData> readsData) {
		for(AssemblyEdge edge: graph.getEdges()) {
			if(edge.isSameSequenceEdge()) continue;
			ReadPathPhasingData d1 = readsData.get(edge.getVertex1().getSequenceIndex());
			ReadPathPhasingData d2 = readsData.get(edge.getVertex2().getSequenceIndex());
			if(d1==null || d2==null) continue;
			if(d1.isOppositePhase(d2)) {
				if(ploidy==1) System.out.println("Filtering by phasing edge: "+edge);
				graph.removeEdge(edge);
			}
		}
		for (AssemblyEmbedded embedded:graph.getAllEmbedded()) {
			ReadPathPhasingData d1 = readsData.get(embedded.getSequenceId());
			ReadPathPhasingData d2 = readsData.get(embedded.getHostId());
			if(d1==null || d2==null) continue;
			if(d1.isOppositePhase(d2)) {
				if(ploidy==1) System.out.println("Filtering by phasing embedded: "+embedded);
				graph.removeEmbedded(embedded);
			}
		}
		
	}
	private AssemblyGraph buildGraph(List<QualifiedSequence> sequences, KmersMap map) {
		AssemblyGraph graph;
		GraphBuilderMinimizers builder = new GraphBuilderMinimizers();
		builder.setKmerLength(kmerLength);
		builder.setWindowLength(windowLength);
		builder.setPloidy(ploidy);
		builder.setNumThreads(numThreads);
		builder.setKmersMap(map);
		builder.setLog(log);
		graph = builder.buildAssemblyGraph(sequences);
		return graph;
	}

	private List<QualifiedSequence> runHomopolymerCompression(List<QualifiedSequence> sequences) {
		List<QualifiedSequence> answer = new ArrayList<>(sequences.size());
		for(int i=0;i<sequences.size();i++) {
			QualifiedSequence seq = sequences.get(i);
			answer.add(compressHomopolymers(seq));
		}
		return answer;
	}
	private QualifiedSequence compressHomopolymers(QualifiedSequence seq) {
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
		return new QualifiedSequence(seq.getName(), new DNAMaskedSequence(compressed));
	}
	/**
	 * Load the sequences of the file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot be opened
	 */
	private List<QualifiedSequence> load(String filename, byte inputFormat, int minReadLength) throws IOException {
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
	private List<QualifiedSequence> loadFasta(String filename, int minReadLength) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		handler.setSequenceType(DNASequence.class);
		List<QualifiedSequence> seqsQL = handler.loadSequences(filename);
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
	private List<QualifiedSequence> loadFastq(String filename, int minReadLength) throws IOException {
		List<QualifiedSequence> sequences = new ArrayList<>();
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			reader.setSequenceType(DNASequence.class);
			//TODO: Option to load quality scores
			reader.setLoadMode(FastqFileReader.LOAD_MODE_WITH_NAME);
			Iterator<RawRead> it = reader.iterator();
			while (it.hasNext()) {
				RawRead read = it.next();
				CharSequence characters = read.getCharacters();
				if(characters.length()>=minReadLength) sequences.add(new QualifiedSequence(read.getName(), characters));
			}
		}
		return sequences;
	}
}
