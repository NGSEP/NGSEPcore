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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import ngsep.alignments.ReadAlignment;
import ngsep.alignments.io.ReadAlignmentFileReader;
import ngsep.assembly.io.AssemblyGraphFileHandler;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.main.io.ParseUtils;
import ngsep.math.Distribution;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class AssemblyGraphStatistics {

	// Constants for default values
	public static final int DEF_MIN_READ_LENGTH = Assembler.DEF_MIN_READ_LENGTH;
	public static final byte READS_FORMAT_FASTQ=KmersExtractor.INPUT_FORMAT_FASTQ;
	public static final byte READS_FORMAT_FASTA=KmersExtractor.INPUT_FORMAT_FASTA;
	public static final String LAYOUT_ALGORITHM_MAX_OVERLAP=Assembler.LAYOUT_ALGORITHM_MAX_OVERLAP;
	public static final String LAYOUT_ALGORITHM_KRUSKAL_PATH=Assembler.LAYOUT_ALGORITHM_KRUSKAL_PATH;
	public static final double DEF_MIN_SCORE_PROPORTION_EDGES = Assembler.DEF_MIN_SCORE_PROPORTION_EDGES;
	
	// Logging and progress
	private Logger log = Logger.getLogger(AssemblyGraphStatistics.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	//Parameters
	private String inputFile = null;
	private String outputFile = null;
	private ReferenceGenome genome = null;
	private String alignmentsFile = null;
	private String layoutAlgorithm=LAYOUT_ALGORITHM_KRUSKAL_PATH;
	private double minScoreProportionEdges = DEF_MIN_SCORE_PROPORTION_EDGES;
	private boolean simulated = false;
	
	//Statistics
	private Distribution distLengthsLayoutReads = new Distribution(0,100000,2000);
	private Distribution distSumLengthsIKbpLayout = new Distribution(0,100000,2000);
	private Distribution distSumLengthsLayout = new Distribution(0,100000,2000);
	private Distribution distLengthsEmbeddedReads = new Distribution(0,100000,2000);
	private Distribution distSumLengthsIKbpEmbedded = new Distribution(0,100000,2000);
	private Distribution distSumLengthsEmbedded = new Distribution(0,100000,2000);
	
	private Distribution distScoresTPEmbedded = new Distribution(0,200000,2000);
	private Distribution distCostsTPEmbedded = new Distribution(0,200000,2000);
	private Distribution distCSKTPEmbedded = new Distribution(0, 30000, 500);
	private Distribution distWCSKTPEmbedded = new Distribution(0, 30000, 500);
	private Distribution distEvidencePropLengthTPEmbedded = new Distribution(0, 1.1, 0.01);
	private Distribution distWCSKPropSelfTPEmbedded = new Distribution(0, 1.1, 0.01);
	private Distribution distIndelsKbpTPEmbedded = new Distribution(0, 50, 1);
	
	private Distribution distScoresFPEmbedded = new Distribution(0,200000,2000);
	private Distribution distCostsFPEmbedded = new Distribution(0,200000,2000);
	private Distribution distCSKFPEmbedded = new Distribution(0, 30000, 500);
	private Distribution distWCSKFPEmbedded = new Distribution(0, 30000, 500);
	private Distribution distEvidencePropLengthFPEmbedded = new Distribution(0, 1.1, 0.01);
	private Distribution distWCSKPropSelfFPEmbedded = new Distribution(0, 1.1, 0.01);
	private Distribution distIndelsKbpFPEmbedded = new Distribution(0, 50, 1);
	
	private Distribution distOverlapsTPPathEdges = new Distribution(0,100000,2000);
	private Distribution distScoresTPPathEdges = new Distribution(0,200000,2000);
	private Distribution distCostsTPPathEdges = new Distribution(0,200000,2000);
	private Distribution distSharedKmersTPPathEdges = new Distribution(0,25000,500);
	private Distribution distCSKTPPathEdges = new Distribution(0, 30000, 500);
	private Distribution distWCSKTPPathEdges = new Distribution(0, 30000, 500);
	private Distribution distEvidencePropLengthTPPathEdges = new Distribution(0, 1.1, 0.01);
	private Distribution distSharedKmersProportionTPPathEdges = new Distribution(0, 1.1, 0.01);
	private Distribution distWCSKProportionTPPathEdges = new Distribution(0, 1.1, 0.01);
	private Distribution distOverlapSDTPPathEdges = new Distribution(0, 500, 10);
	private Distribution distNumIndelsTPPathEdges = new Distribution(0,500,10);
	private Distribution distIndelsKbpTPPathEdges = new Distribution(0,50,1);
	
	
	private Distribution distOverlapsFPEdges = new Distribution(0,100000,2000);
	private Distribution distScoresFPEdges = new Distribution(0,200000,2000);
	private Distribution distCostsFPEdges = new Distribution(0,200000,2000);
	private Distribution distSharedKmersFPEdges = new Distribution(0,25000,500);
	private Distribution distCSKFPEdges = new Distribution(0, 30000, 500);
	private Distribution distWCSKFPEdges = new Distribution(0, 30000, 500);
	private Distribution distEvidencePropLengthFPEdges = new Distribution(0, 1.1, 0.01);
	private Distribution distSharedKmersProportionFPEdges = new Distribution(0, 1.1, 0.01);
	private Distribution distWCSKProportionFPEdges = new Distribution(0, 1.1, 0.01);
	private Distribution distOverlapSDFPEdges = new Distribution(0, 500, 10);
	private Distribution distNumIndelsFPEdges = new Distribution(0,500,10);
	private Distribution distIndelsKbpFPEdges = new Distribution(0,50,1);
	
	private Distribution distOverlapsFPPathEdges = new Distribution(0,100000,2000);
	private Distribution distScoresFPPathEdges = new Distribution(0,200000,2000);
	private Distribution distCostsFPPathEdges = new Distribution(0,200000,2000);
	private Distribution distCSKFPPathEdges = new Distribution(0, 30000, 500);
	private Distribution distWCSKFPPathEdges = new Distribution(0, 30000, 500);
	private Distribution distEvidencePropLengthFPPathEdges = new Distribution(0, 1.1, 0.01);
	private Distribution distWCSKProportionFPPathEdges = new Distribution(0, 1.1, 0.01);
	private Distribution distNumIndelsFPPathEdges = new Distribution(0,500,10);
	private Distribution distIndelsKbpFPPathEdges = new Distribution(0,50,1);
	
	private Distribution distOverlapsFNPathEdges = new Distribution(0,100000,2000);
	
	private int tpEmbSeqs = 0;
	private int fpEmbSeqs = 0;
	private int fnEmbSeqs = 0;
	private int tpEmbRel = 0;
	private int fpEmbRel = 0;
	private int fnEmbRel = 0;
	private int tpEdgesNotEmbedded = 0;
	private int tpEdgesEmbedded = 0;
	private int fpEdges = 0;
	private int fnEdgesNotEmbedded = 0;
	private int fnEdgesEmbedded = 0;
	private int tpPathEdges = 0;
	private int totalPathEdges = 0;
	
	private int tpLayoutEdges = 0;
	private int errorsTPEdgeNoLayout = 0;
	private int errorsEdgeEmbeddedNoLayout = 0;
	private int errorsFPEdge = 0;
	private int errorsFNLayoutEdge = 0;
	private int totalTestLayoutPaths = 0;
	private int totalTestLayoutEdges = 0;
	private int totalGSLayoutEdges = 0;
	
	private double rmsePredictedOverlap=0;
	private int countPredictedOverlap = 0;
	private Distribution distOverlapError = new Distribution(-500,500,50);
	private double rmseAveragePredictedOverlap=0;
	private int countAveragePredictedOverlap = 0;
	private Distribution distAverageOverlapError = new Distribution(-500,500,50);
	private double rmseMedianPredictedOverlap=0;
	private int countMedianPredictedOverlap = 0;
	private Distribution distMedianOverlapError = new Distribution(-500,500,50);
	private double rmseFromLimitsPredictedOverlap=0;
	private int countFromLimitsPredictedOverlap = 0;
	private Distribution distFromLimitsOverlapError = new Distribution(-500,500,50);
	
	private boolean filteredGraph = false;
	private boolean logErrors = false;
	
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
	public String getOutputFile() {
		return outputFile;
	}
	public void setOutputFile(String outputFile) {
		this.outputFile = outputFile;
	}
	
	public ReferenceGenome getGenome() {
		return genome;
	}
	public void setGenome(ReferenceGenome genome) {
		this.genome = genome;
	}
	public void setGenome(String genomeFile) throws IOException {
		setGenome(OptionValuesDecoder.loadGenome(genomeFile,log));
	}
	public String getAlignmentsFile() {
		return alignmentsFile;
	}
	public void setAlignmentsFile(String alignmentsFile) {
		this.alignmentsFile = alignmentsFile;
	}
	
	
	public String getLayoutAlgorithm() {
		return layoutAlgorithm;
	}
	public void setLayoutAlgorithm(String layoutAlgorithm) {
		this.layoutAlgorithm = layoutAlgorithm;
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
	
	public boolean isSimulated() {
		return simulated;
	}
	public void setSimulated(boolean simulated) {
		this.simulated = simulated;
	}
	public void setSimulated(Boolean simulated) {
		setSimulated(simulated.booleanValue());
	}
	
	public static void main(String[] args) throws Exception {
		AssemblyGraphStatistics instance = new AssemblyGraphStatistics();
		CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.run();

	}
	public void run() throws IOException {
		logParameters();
		if(inputFile==null) throw new IOException("The input graph is required");
		if(outputFile==null) throw new IOException("An output file path is required");
		//if(!simulated && readsFile==null && alignmentsFile==null) throw new IOException("For non simulated reads either the original reads or the alignments are required");
		run (inputFile, outputFile);
		log.info("Process finished");
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		if(simulated) out.println("Reads were simulated using SingleReadsSimulator");
		out.println("Layout algorithm:"+ layoutAlgorithm);
		if (genome!=null) out.println("Target genome for benchmark loaded from file: "+genome.getFilename());
		if (alignmentsFile!=null) out.println("Alignments file for benchmark "+alignmentsFile);
		out.println("Minimum score proportion (from the maximum score) to keep edges of a sequence: "+ minScoreProportionEdges);
		log.info(os.toString());
	}
	private void run(String inputFile, String outputFile) throws IOException {
		List<ReadAlignment> alignments = null;
		List<QualifiedSequence> sequences = AssemblyGraphFileHandler.loadSequenceNamesFromGraphFile(inputFile);
		AssemblyGraph graph = AssemblyGraphFileHandler.load(sequences, inputFile);
		Map<String,Integer> seqIds = new HashMap<String, Integer>();
		for(int i=0;i<sequences.size();i++) {
			QualifiedSequence seq = sequences.get(i);
			if(seqIds.containsKey(seq.getName())) {
				log.warning("Duplicated read name: "+seq.getName());
			} else {
				seqIds.put(seq.getName(), i);
			}
		}
		if (alignmentsFile!=null) {
			alignments = new ArrayList<ReadAlignment>();
			try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(alignmentsFile)) {
				log.info("Loading alignments from: "+alignmentsFile);
				reader.setLoadMode(ReadAlignmentFileReader.LOAD_MODE_ALIGNMENT_NAME);
				reader.setFilterFlags(ReadAlignment.FLAG_SECONDARY);
				Iterator<ReadAlignment> it = reader.iterator();
				while (it.hasNext()) {
					ReadAlignment aln = it.next();
					if(aln.getReadLength()<Assembler.DEF_MIN_READ_LENGTH) continue;
					if(aln.isReadUnmapped()) continue;
					Integer idx = seqIds.get(aln.getReadName());
					if(idx==null) {
						log.warning("Aligned read: "+aln.getReadName()+" not found in graph");
						continue;
					}
					else aln.setReadNumber(idx);
					QualifiedSequence refSeq = genome.getSequenceByName(aln.getSequenceName());
					if(refSeq == null) {
						log.warning("Reference sequence not found for alignment: "+aln);
						continue;
					}
					int refSeqLength = refSeq.getLength();
					
					if(aln.getFirst()-aln.getSoftClipStart()>-1000 && aln.getLast()+aln.getSoftClipEnd()<refSeqLength+1000 && aln.getSoftClipStart()>3000 || aln.getSoftClipEnd()>3000) {
						//log.warning("Alignment of read idx "+idx+ " name "+aln.getReadName()+" has a large soft clip: "+aln.getSoftClipStart()+" "+aln.getSoftClipEnd()+" aligned to "+aln.getSequenceName()+":"+aln.getFirst()+"-"+aln.getLast());
						continue;
					}
					
					//sequences.add(new QualifiedSequence(aln.getReadName(), characters));
					
					String newReadName = aln.getSequenceName()+"_"+aln.getFirst()+"_"+(aln.isNegativeStrand()?1:0);
					aln.setReadName(newReadName);
					sequences.get(idx).setName(newReadName);
					alignments.add(aln);
				}
				
			}
		} else if (simulated) {
			log.info("Building simulated alignments");
			alignments = buildAlignmentsFromSimulatedReads(sequences);
		}
		log.info("Loaded "+alignments.size()+" alignments");
		AssemblyGraph goldStandardGraph = null;
		if(alignments!=null) goldStandardGraph = buildGoldStandardGraph(alignments, sequences);
		try (PrintStream out=new PrintStream(outputFile)) {
			Distribution degreeDist = graph.getVertexDegreeDistribution ();
			out.println("Vertex degree distribution");
			degreeDist.printDistributionInt(out);
			//Infer distributions to calculate costs
			log.info("Updating scores");
			graph.updateScores(0);
			log.info("Comparing initial graph");
			if(goldStandardGraph!=null) compareGraphs(goldStandardGraph, graph, out);
			out.println("Initial graph statistics. Vertices: "+graph.getVertices().size()+" edges: "+graph.getNumEdges());
			printStatistics(out);
			resetStatistics();
			log.info("Removing vertices chimeric reads");
			graph.removeVerticesChimericReads();
			log.info("Filtered chimeric reads. Vertices: "+graph.getVertices().size()+" edges: "+graph.getEdges().size());
			//findProblematicVertices(graph);
			log.info("Filtering edges and embedded");
			(new AssemblySequencesRelationshipFilter()).filterEdgesAndEmbedded(graph, minScoreProportionEdges);
			//log.info("Updating scores after filtering");
			//graph.updateScores();
			filteredGraph = true;
			log.info("Building layout");
			LayoutBuilder pathsFinder;
			if(LAYOUT_ALGORITHM_MAX_OVERLAP.equals(layoutAlgorithm)) {
				//LayoutBuilder pathsFinder = new LayoutBuilderGreedyMinCost();
				//LayoutBuilder pathsFinder = new LayoutBuilderGreedyMaxCoverageSharedKmers();
				pathsFinder = new LayoutBuilderGreedyMaxOverlap();
			} else {
				LayoutBuilderKruskalPath kruskal = new LayoutBuilderKruskalPath();
				//kruskal.setRunImprovementAlgorithms(false);
				pathsFinder = kruskal;
				
			}
			//System.out.println("Searched vertex: "+graph.getVertex(1222, true));
			//System.out.println("Searched vertex: "+graph.getVertex(963, false));
			//System.out.println("Searched edge. "+graph.getEdge(graph.getVertex(963, false), graph.getVertex(1222, true)));
			pathsFinder.findPaths(graph);
			if(goldStandardGraph!=null) {
				//logErrors=true;
				log.info("Comparing graph and layout");
				compareGraphs(goldStandardGraph, graph, out);
				//logErrors = true;
				compareLayouts(goldStandardGraph, graph, out);
			}
			out.println("Filtered graph statistics");
			printStatistics(out);
		}
	}
	
	public List<ReadAlignment> buildAlignmentsFromSimulatedReads(List<QualifiedSequence> sequences) {
		//Create true alignments of simulated reads to the target genome
		QualifiedSequenceList seqNames = genome.getSequencesMetadata();
		List<ReadAlignment> alignments = new ArrayList<>();
		List<String> readNames = new ArrayList<String>(sequences.size());
		for(int i=0;i<sequences.size();i++) {
			QualifiedSequence seq = sequences.get(i);
			String readName = seq.getName();
			readNames.add(readName);
			int i1 = -1;
			int i2 = -1;
			int i3 = -1;
			for(int j=readName.length()-1;j>=0;j--) {
				if(readName.charAt(j)!='_') continue;
				if(i3==-1) i3 = j;
				else if(i2==-1) i2 = j;
				else if (i1==-1) {
					i1 = j;
					break;
				}
			}
			QualifiedSequence seqName = seqNames.get(readName.substring(0,i1));
			int first = Integer.parseInt(readName.substring(i1+1, i2));
			boolean reverse = readName.charAt(i2+1)=='1';
			int flags = 0;
			if (reverse) flags = ReadAlignment.FLAG_READ_REVERSE_STRAND;
			//System.out.println("Next sequence: "+readName+" first: "+first+" reverse: "+reverse+" flags: "+flags+" seq: "+seqName);
			ReadAlignment aln = new ReadAlignment(seqName.getName(), first, first+seq.getLength()-1, seq.getLength(), flags);
			aln.setReadNumber(i);
			aln.setReadName(readName);
			aln.setCigarString(seq.getLength()+"M");
			alignments.add(aln);
		}
		return alignments;
	}
	private AssemblyGraph buildGoldStandardGraph(List<ReadAlignment> alignments, List<QualifiedSequence> sequences) {
		log.info("Building gold standard graph from "+alignments.size()+" alignments on "+sequences.size()+" sequences");
		AssemblyGraph graph = new AssemblyGraph(sequences);
		//Sort by target genome location to calculate edges efficiently
		GenomicRegionComparator comparator = new GenomicRegionComparator(genome.getSequencesMetadata());
		Collections.sort(alignments, comparator);
		for(int i=0;i<alignments.size();i++) {
			ReadAlignment left = alignments.get(i);
			boolean debug = left.getReadNumber()==-1;
			QualifiedSequence leftSeq = new QualifiedSequence(left.getReadName());
			leftSeq.setLength(left.getReadLength());
			int leftStart = left.getFirst()-left.getSoftClipStart();
			int leftEnd = left.getLast()+left.getSoftClipEnd();
			if(debug) System.out.println("GoldStandardGraph. Left aln: "+left+" limits: "+leftStart+" "+leftEnd);
			AssemblyVertex vertexLeft = graph.getVertex(left.getReadNumber(), left.isNegativeStrand());
			for(int j=i+1;j<alignments.size();j++) {
				ReadAlignment right = alignments.get(j);
				int cmp = comparator.compare(right, left);
				if(debug) System.out.println("GoldStandardGraph. Next aln to compare: "+right+" cmp: "+cmp);
				if(cmp>1) break;
				QualifiedSequence rightSeq = new QualifiedSequence(right.getReadName());
				rightSeq.setLength(right.getReadLength());
				int rightStart = right.getFirst()-right.getSoftClipStart();
				int rightEnd = right.getLast()+right.getSoftClipEnd();
				if(rightStart>=leftStart && rightEnd > leftEnd-50) {
					AssemblyVertex vertexRight = graph.getVertex(right.getReadNumber(), !right.isNegativeStrand());
					int overlap = leftEnd - rightStart + 1;
					if(debug) System.out.println("GoldStandardGraph. Left right edge Next aln: "+right+" overlap: "+overlap+" limits: "+rightStart+" "+rightEnd);
					AssemblyEdge edge = new AssemblyEdge(vertexLeft, vertexRight, overlap);
					edge.setAverageOverlap(overlap);
					edge.setMedianOverlap(overlap);
					edge.setFromLimitsOverlap(overlap);
					graph.addEdge(edge);
				} else if (rightStart<leftStart && rightEnd < leftEnd) {
					AssemblyVertex vertexRight = graph.getVertex(right.getReadNumber(), !right.isNegativeStrand());
					int overlap = rightEnd - leftStart + 1;
					if(debug) System.out.println("GoldStandardGraph. Right left edge. Next aln: "+right+" overlap: "+overlap+" limits: "+rightStart+" "+rightEnd);
					AssemblyEdge edge = new AssemblyEdge(vertexRight, vertexLeft, overlap);
					edge.setAverageOverlap(overlap);
					edge.setMedianOverlap(overlap);
					edge.setFromLimitsOverlap(overlap);
					graph.addEdge(edge);
				}
				
				boolean relativeNegative = left.isNegativeStrand()!=right.isNegativeStrand();
				int relativeStart;
				
				if(leftStart >= rightStart && leftEnd<=rightEnd) {
					//left is embedded in right
					int relativeEnd;
					if(right.isNegativeStrand()) {
						relativeStart = left.getLast() - right.getFirst();
						relativeEnd = 0;
					} else {
						relativeStart = 0;
						relativeEnd = left.getLast() - right.getFirst();
					}
					AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(left.getReadNumber(), leftSeq, relativeNegative, right.getReadNumber(), rightSeq, relativeStart, relativeEnd);
					graph.addEmbedded(embeddedEvent);
					if(debug) System.out.println("GoldStandardGraph. Added embedded relationship right - left");
					if(leftStart == rightStart && rightEnd<=leftEnd) {
						//right is also embedded in left
						if(left.isNegativeStrand()) {
							relativeStart = right.getLast() - left.getFirst();
							relativeEnd = 0;
						} else {
							relativeStart = 0;
							relativeEnd = right.getLast() - left.getFirst();
						}
						embeddedEvent = new AssemblyEmbedded(right.getReadNumber(), rightSeq, relativeNegative, left.getReadNumber(), leftSeq, relativeStart, relativeEnd );
						graph.addEmbedded(embeddedEvent);
						if(debug) System.out.println("GoldStandardGraph. Added reciprocal embedded relationship left-right");
					}
				} else if(leftStart <= rightStart && rightEnd<=leftEnd) {
					//Right is embedded in left
					int relativeEnd;
					if(left.isNegativeStrand()) {
						relativeStart = right.getLast() - left.getFirst();
						relativeEnd = right.getFirst()-left.getFirst();
					} else {
						relativeStart = right.getFirst()-left.getFirst();
						relativeEnd = right.getLast() - left.getFirst();
					}
					AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(right.getReadNumber(), rightSeq, relativeNegative, left.getReadNumber(), leftSeq, relativeStart, relativeEnd );
					graph.addEmbedded(embeddedEvent);
					if(debug) System.out.println("GoldStandardGraph. Added embedded relationship left-right");
				}
			}
		}
		log.info("Created gold standard assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges. Embedded: "+graph.getEmbeddedCount());
		//Build gold standard layouts it must be done after knowing which sequences are embedded
		String lastSeqName = null;
		AssemblyPath nextPath = null;
		AssemblyVertex lastVertex = null;
		for(int i=0;i<alignments.size();i++) {
			ReadAlignment aln = alignments.get(i);
			if(graph.isEmbedded(aln.getReadNumber())) continue;
			if(!aln.getSequenceName().equals(lastSeqName)) {
				if(nextPath!=null) graph.addPath(nextPath);
				lastSeqName = aln.getSequenceName();
				nextPath = null;
				lastVertex = null;
			}
			AssemblyVertex leftSequenceVertex = graph.getVertex(aln.getReadNumber(), !aln.isNegativeStrand());
			AssemblyVertex rightSequenceVertex = graph.getVertex(aln.getReadNumber(), aln.isNegativeStrand());
			if(lastVertex!=null) {
				AssemblyEdge connectingEdge = graph.getEdge(lastVertex, leftSequenceVertex);
				if(connectingEdge==null) {
					log.info("Discontiguity in gold standard layout. Last vertex: "+lastVertex+" next vertex: "+leftSequenceVertex);
					if(nextPath!=null) graph.addPath(nextPath);
					nextPath = null;
					lastVertex = null;
				} else {
					//System.out.println("Next layout edge between "+lastVertex.getSequenceIndex()+" start " +lastVertex.isStart()+" and "+leftSequenceVertex.getSequenceIndex()+" start "+leftSequenceVertex.isStart()+" read id: "+aln.getSequenceName()+" first: "+aln.getFirst());
					AssemblyEdge edgeSequence = graph.getSameSequenceEdge(lastVertex);
					if(nextPath==null) nextPath = new AssemblyPath(edgeSequence);
					nextPath.connectEdgeRight(graph, connectingEdge);
					connectingEdge.setLayoutEdge(true);
				}
			}
			lastVertex = rightSequenceVertex;
		}
		if(nextPath!=null) graph.addPath(nextPath);
		return graph;
	}

	private void compareGraphs(AssemblyGraph goldStandardGraph, AssemblyGraph testGraph, PrintStream out) {
		
		int n = goldStandardGraph.getNumSequences();
		for(int i=0;i<n;i++) {
			QualifiedSequence sequence = goldStandardGraph.getSequence(i);
			AssemblyEdge selfEdge = testGraph.getSameSequenceEdge(i);
			//Check embedded status
			boolean gsE = goldStandardGraph.isEmbedded(i);
			if(gsE) distLengthsEmbeddedReads.processDatapoint(goldStandardGraph.getSequenceLength(i));
			boolean testE = testGraph.isEmbedded(i);
			if(gsE && testE) {
				tpEmbSeqs++;
				List<AssemblyEmbedded> hosts = testGraph.getEmbeddedBySequenceId(i);
				int maxScore = 0;
				int minCost = 10000000;
				double maxEvidenceProp = 0;
				double maxCSK = 0;
				double maxWCSK = 0;
				double minIndelsKbp = 1000000;
				AssemblyEmbedded minCostE = null;
				for(AssemblyEmbedded embedded: hosts) {
					maxScore = Math.max(maxScore, embedded.getScore());
					if(minCostE==null || minCostE.getCost()>embedded.getCost()) minCostE = embedded;
					minCost = Math.min(minCost, embedded.getCost());
					
					maxEvidenceProp = Math.max(maxEvidenceProp, embedded.getEvidenceProportion());
					maxCSK = Math.max(maxCSK, embedded.getCoverageSharedKmers());
					maxWCSK = Math.max(maxWCSK, embedded.getWeightedCoverageSharedKmers());
					minIndelsKbp = Math.min(minIndelsKbp, embedded.getIndelsPerKbp());
				}
				distScoresTPEmbedded.processDatapoint(maxScore);
				distCostsTPEmbedded.processDatapoint(minCost);
				distCSKTPEmbedded.processDatapoint(maxCSK);
				distWCSKTPEmbedded.processDatapoint(maxWCSK);
				distEvidencePropLengthTPEmbedded.processDatapoint(maxEvidenceProp);
				distIndelsKbpTPEmbedded.processDatapoint(minIndelsKbp);
				if(minIndelsKbp>0) {
					int totalLength = testGraph.getSequenceLength(minCostE.getSequenceId())+testGraph.getSequenceLength(minCostE.getHostId());
					distSumLengthsEmbedded.processDatapoint(totalLength);
					distSumLengthsIKbpEmbedded.processDatapoint(minIndelsKbp,totalLength);
					/*if(minIndelsKbp>10) {
						System.out.println("Large indels per KBP for true embedded. Relationships: ");
						for(AssemblyEmbedded embedded: hosts) System.out.println(""+embedded);
					}*/
				}
				if(selfEdge!=null) {
					double p2 = maxWCSK/selfEdge.getCoverageSharedKmers();
					distWCSKPropSelfTPEmbedded.processDatapoint(p2);
				}
				//if(minCostE!=null && minCost>200000) System.out.println("TPEmbedded with high cost: "+minCostE+" score: "+minCostE.getScore()+" cost: "+minCostE.getCost()+" host: "+logSequence(minCostE.getHostId(), goldStandardGraph.getSequence(minCostE.getHostId())));
			}
			else if (gsE) {
				fnEmbSeqs++;
				if (logErrors) System.err.println("Embedded sequence not called: "+logSequence(i, sequence));
				double maxEvidenceProp = 0;
				for(AssemblyEmbedded embedded:goldStandardGraph.getEmbeddedBySequenceId(i)) {
					maxEvidenceProp = Math.max(maxEvidenceProp, embedded.getEvidenceProportion());
					QualifiedSequence seqHost = goldStandardGraph.getSequence(embedded.getHostId()) ;
					if (logErrors) System.err.println("Next true host "+logSequence(embedded.getHostId(),seqHost));
				}
			}
			else if (testE) {
				fpEmbSeqs++;
				List<AssemblyEmbedded> falseHosts = testGraph.getEmbeddedBySequenceId(i);
				int maxScore = 0;
				AssemblyEmbedded relMaxScore = null;
 				int minCost = 10000000;
				double maxEvidenceProp = 0;
				double maxCSK = 0;
				double maxWCSK = 0;
				double minIndelsKbp = 1000000;
				if (logErrors) System.err.println("False embedded sequence "+logSequence(i, sequence)+" false hosts: "+falseHosts.size());
				for(AssemblyEmbedded embedded:falseHosts) {
					QualifiedSequence seqHost = testGraph.getSequence(embedded.getHostId()) ;
					if (logErrors) System.err.println("Next false relation with host: "+seqHost.getName()+" length: "+seqHost.getLength()+" "+embedded);
					if(relMaxScore==null || maxScore<embedded.getScore()) {
						relMaxScore = embedded;
						maxScore = embedded.getScore();
					}
					minCost = Math.min(minCost, embedded.getCost());
					maxEvidenceProp = Math.max(maxEvidenceProp, embedded.getEvidenceProportion());
					maxCSK = Math.max(maxCSK, embedded.getCoverageSharedKmers());
					maxWCSK = Math.max(maxWCSK, embedded.getWeightedCoverageSharedKmers());
					minIndelsKbp = Math.min(minIndelsKbp, embedded.getIndelsPerKbp());
					
				}
				distScoresFPEmbedded.processDatapoint(maxScore);
				distCostsFPEmbedded.processDatapoint(minCost);
				distCSKFPEmbedded.processDatapoint(maxCSK);
				distWCSKFPEmbedded.processDatapoint(maxWCSK);
				distEvidencePropLengthFPEmbedded.processDatapoint(maxEvidenceProp);
				distIndelsKbpFPEmbedded.processDatapoint(minIndelsKbp);
				if(selfEdge!=null) {
					double p2 = maxCSK/selfEdge.getCoverageSharedKmers();
					distWCSKPropSelfFPEmbedded.processDatapoint(p2);
				}
				int totalEdgesFP = 0;
				AssemblyVertex v1 = testGraph.getVertex(i, true);
				if(v1!=null) totalEdgesFP += testGraph.getEdges(v1).size();
				AssemblyVertex v2 = testGraph.getVertex(i, false);
				if(v2!=null) totalEdgesFP += testGraph.getEdges(v2).size();
				if(!filteredGraph && totalEdgesFP==0) System.err.println("Zero edges for false positive embedded "+relMaxScore);
			}
			//Check embedded relationships
			List<AssemblyEmbedded> embeddedGS = goldStandardGraph.getEmbeddedByHostId(i);
			List<AssemblyEmbedded> embeddedTest = testGraph.getEmbeddedByHostId(i);
			int tpM = calculateIntersection (embeddedGS,embeddedTest);
			tpEmbRel+=tpM;
			fpEmbRel+=(embeddedTest.size()-tpM);
			fnEmbRel+=(embeddedGS.size()-tpM);
			
			AssemblyVertex testVertex = testGraph.getVertex(i, true);
			AssemblyVertex gsVertex = goldStandardGraph.getVertex(i, true);
			if(!validateEqualVertices(i, sequence, testVertex, gsVertex)) continue;
			
			calculateComparisonStats (goldStandardGraph, gsVertex, testGraph, testVertex, out);
			testVertex = testGraph.getVertex(i, false);
			gsVertex = goldStandardGraph.getVertex(i, false);
			if(!validateEqualVertices(i, sequence, testVertex, gsVertex)) continue;
			
			calculateComparisonStats (goldStandardGraph, gsVertex, testGraph, testVertex, out);
		}
	}
	private int calculateIntersection(List<AssemblyEmbedded> embeddedGS, List<AssemblyEmbedded> embeddedTest) {
		int count = 0;
		for(AssemblyEmbedded eGS:embeddedGS) {
			for(AssemblyEmbedded eTest:embeddedTest) {
				if(eGS.getSequenceId()== eTest.getSequenceId() && eGS.isReverse()== eTest.isReverse()) {
					count++;
					break;
				}
			}
		}
		return count;
	}

	private boolean validateEqualVertices(int seqIndex, QualifiedSequence sequence, AssemblyVertex testVertex, AssemblyVertex gsVertex) {
		if(testVertex==null) {
			return false;
		}
		if(!testVertex.getRead().getName().equals(sequence.getName())) {
			log.warning("Inconsistent sequence for test vertex. test name: "+logSequence(seqIndex, testVertex.getRead())+" expected: "+logSequence(seqIndex, sequence));
			return false;
		}
		if(gsVertex==null) {
			log.warning("Gold standard vertex not found for start of sequence "+logSequence(seqIndex, sequence));
			return false;
		}
		if(!gsVertex.getRead().getName().equals(sequence.getName())) {
			log.warning("Inconsistent sequence for gold standard vertex. test name: "+logSequence(seqIndex, gsVertex.getRead())+" expected: "+logSequence(seqIndex, sequence));
			return false;
		}
		if(testVertex.getUniqueNumber()!=gsVertex.getUniqueNumber()) {
			log.warning("Inconsistent number ids for test and gold standard vertices. test number: "+testVertex.getUniqueNumber()+" gold standard number: "+gsVertex.getUniqueNumber());
			return false;
		}
		return true;
	}
	private String logSequence(int idx, QualifiedSequence sequence) {
		return ""+idx+" "+sequence.getName()+" "+sequence.getLength();
	}
	private void calculateComparisonStats(AssemblyGraph goldStandardGraph, AssemblyVertex gsVertex, AssemblyGraph testGraph,  AssemblyVertex testVertex, PrintStream out) {
		//Find path edge of this vertex
		List<AssemblyEdge> gsEdges = goldStandardGraph.getEdges(gsVertex);
		List<AssemblyEdge> testEdges = testGraph.getEdges(testVertex);
		boolean debug = gsVertex.getSequenceIndex()==-1;
		//boolean debug = gsVertex.getSequenceIndex()==513290 || gsVertex.getSequenceIndex()== 213638 || gsVertex.getSequenceIndex()==1267;
		//boolean debug = gsVertex.getSequenceIndex()==2123 || gsVertex.getSequenceIndex()==4681 || gsVertex.getSequenceIndex()==239; 
		if(debug) {
			printEdgeList("Gold standard", gsVertex, gsEdges, goldStandardGraph, false, false, out);
			printEdgeList("Test", testVertex, testEdges, testGraph, true, true, out);
		}
		Map<Integer,Boolean> testEdgesMatched = new HashMap<Integer, Boolean>();
		Map<Integer,AssemblyEdge> testEdgesByConnectingVertex = new HashMap<Integer, AssemblyEdge>();
		for(AssemblyEdge edge:testEdges) {
			testEdgesMatched.put(edge.getConnectingVertex(testVertex).getUniqueNumber(), false);
			testEdgesByConnectingVertex.put(edge.getConnectingVertex(testVertex).getUniqueNumber(),edge);
		}
		boolean vertexEmbedded = testGraph.isEmbedded(gsVertex.getSequenceIndex());
		for(AssemblyEdge gsEdge:gsEdges) {
			AssemblyVertex gsConnectingVertex = gsEdge.getConnectingVertex(gsVertex);
			boolean edgeEmbeddedInTest = vertexEmbedded || testGraph.isEmbedded(gsConnectingVertex.getSequenceIndex());
			int number = gsConnectingVertex.getUniqueNumber();
			AssemblyEdge matchedTestEdge = testEdgesByConnectingVertex.get(number);
			boolean match = testEdgesMatched.containsKey(number);
			if(match) {
				//True positive
				if(edgeEmbeddedInTest) tpEdgesEmbedded++;
				else {
					tpEdgesNotEmbedded++;
					if(!gsEdge.isSameSequenceEdge()) updateOverlapStats(gsEdge, matchedTestEdge);
				}
				testEdgesMatched.put(number, true);
			}
			else {
				//False negative
				if(edgeEmbeddedInTest) fnEdgesEmbedded++;
				else {
					fnEdgesNotEmbedded++;
					//out.println("False negative edge between vertex: "+logVertex(gsVertex)+" and "+logVertex(gsConnectingVertex)+" overlap: "+gsEdge.getOverlap());
				}
			}
			/*if(debug) {
				out.println("Next gs: "+number+" match: "+match+" TP: "+answer[0]+" test entries: "+testEdgesMatched);
			}*/
			if(gsEdge.isLayoutEdge()) {
				if(match) {
					tpPathEdges++;
							//distSumLengthsIKbpLayout
					if (!gsEdge.isSameSequenceEdge()) {
						distOverlapsTPPathEdges.processDatapoint(matchedTestEdge.getOverlap());
						distScoresTPPathEdges.processDatapoint(matchedTestEdge.getScore());
						distCostsTPPathEdges.processDatapoint(matchedTestEdge.getCost());
						distSharedKmersTPPathEdges.processDatapoint(matchedTestEdge.getNumSharedKmers());
						distCSKTPPathEdges.processDatapoint(matchedTestEdge.getCoverageSharedKmers());
						distWCSKTPPathEdges.processDatapoint(matchedTestEdge.getWeightedCoverageSharedKmers());
						distEvidencePropLengthTPPathEdges.processDatapoint(matchedTestEdge.getEvidenceProportion());
						distSharedKmersProportionTPPathEdges.processDatapoint((double)matchedTestEdge.getNumSharedKmers()/(matchedTestEdge.getOverlap()+1));
						distWCSKProportionTPPathEdges.processDatapoint((double)matchedTestEdge.getWeightedCoverageSharedKmers()/(matchedTestEdge.getOverlap()+1));
						distOverlapSDTPPathEdges.processDatapoint(matchedTestEdge.getRawKmerHitsSubjectStartSD());
						distNumIndelsTPPathEdges.processDatapoint(matchedTestEdge.getNumIndels());
						distIndelsKbpTPPathEdges.processDatapoint(matchedTestEdge.getIndelsPerKbp());
						if(matchedTestEdge.getIndelsPerKbp()>7) log.info("Large indels per kbp for path edge: "+matchedTestEdge);
						if(matchedTestEdge.getEvidenceProportion()<0.9) log.info("Low evidence proportion for path edge: "+matchedTestEdge);
						int lengthSum = testGraph.getSequenceLength(matchedTestEdge.getVertex1().getSequenceIndex())+testGraph.getSequenceLength(matchedTestEdge.getVertex2().getSequenceIndex());
						distSumLengthsIKbpLayout.processDatapoint(matchedTestEdge.getIndelsPerKbp(), lengthSum);
						distSumLengthsLayout.processDatapoint(lengthSum);
					}
				} else {
					distOverlapsFNPathEdges.processDatapoint(gsEdge.getOverlap());
					//log.info("Path edge not matched. Embedded vertices in test. "+edgeEmbeddedInTest+" edge: "+gsEdge);
				}
				totalPathEdges++;
				if (gsEdge.isSameSequenceEdge()) distLengthsLayoutReads.processDatapoint(goldStandardGraph.getSequenceLength(gsVertex.getSequenceIndex()));
			}
		}
		for(AssemblyEdge edge:testEdges) {
			AssemblyVertex vertex = edge.getConnectingVertex(testVertex); 
			int number = vertex.getUniqueNumber();
			if(!testEdgesMatched.get(number)) {
				//False positive
				fpEdges++;
				distOverlapsFPEdges.processDatapoint(edge.getOverlap());
				distScoresFPEdges.processDatapoint(edge.getScore());
				distCostsFPEdges.processDatapoint(edge.getCost());
				distSharedKmersFPEdges.processDatapoint(edge.getNumSharedKmers());
				distCSKFPEdges.processDatapoint(edge.getCoverageSharedKmers());
				distWCSKFPEdges.processDatapoint(edge.getWeightedCoverageSharedKmers());
				distEvidencePropLengthFPEdges.processDatapoint(edge.getEvidenceProportion());
				distSharedKmersProportionFPEdges.processDatapoint((double)edge.getNumSharedKmers()/(edge.getOverlap()+1));
				distWCSKProportionFPEdges.processDatapoint((double)edge.getWeightedCoverageSharedKmers()/(edge.getOverlap()+1));
				distOverlapSDFPEdges.processDatapoint(edge.getRawKmerHitsSubjectStartSD());
				distNumIndelsFPEdges.processDatapoint(edge.getNumIndels());
				distIndelsKbpFPEdges.processDatapoint(edge.getIndelsPerKbp());
				if (logErrors) System.err.println("False positive edge "+edge);
			}
		}
	}
	private void updateOverlapStats(AssemblyEdge gsEdge, AssemblyEdge testEdge) {
		if(testEdge==null || testEdge.isSameSequenceEdge()) return;
		double error = gsEdge.getOverlap()-testEdge.getOverlap();
		rmsePredictedOverlap+=error*error;
		countPredictedOverlap++;
		distOverlapError.processDatapoint(error);
		if(logErrors && error < -200 && gsEdge.isLayoutEdge()) System.out.println("Large overlap error in edge. GS edge: "+gsEdge+"\ntest edge: "+testEdge );
		error = gsEdge.getOverlap()-testEdge.getAverageOverlap();
		rmseAveragePredictedOverlap+=error*error;
		countAveragePredictedOverlap++;
		distAverageOverlapError.processDatapoint(error);
		error = gsEdge.getOverlap()-testEdge.getMedianOverlap();
		rmseMedianPredictedOverlap+=error*error;
		countMedianPredictedOverlap++;
		distMedianOverlapError.processDatapoint(error);
		error = gsEdge.getOverlap()-testEdge.getFromLimitsOverlap();
		rmseFromLimitsPredictedOverlap+=error*error;
		countFromLimitsPredictedOverlap++;
		distFromLimitsOverlapError.processDatapoint(error);
	}
	private boolean sameEdges(AssemblyEdge nextGSEdge, AssemblyEdge nextTestEdge) {
		AssemblyVertex testV1 = nextTestEdge.getVertex1();
		AssemblyVertex testV2 = nextTestEdge.getVertex2();
		AssemblyVertex gsV1 = nextGSEdge.getVertex1();
		AssemblyVertex gsV2 = nextGSEdge.getVertex2();
		if (testV1.getUniqueNumber()==gsV1.getUniqueNumber()) {
			if(testV2.getUniqueNumber()==gsV2.getUniqueNumber()) return true;
		} else if (testV1.getUniqueNumber()==gsV2.getUniqueNumber()) {
			if(testV2.getUniqueNumber()==gsV1.getUniqueNumber()) return true;
		}
		return false;
	}

	public void printEdgeList(String text, AssemblyVertex v, List<AssemblyEdge> edges, AssemblyGraph graph, boolean includeEmbedded, boolean sortByCost, PrintStream out) {
		out.println(text+" vertex "+v);
		List<AssemblyEdge> copy = new ArrayList<AssemblyEdge>();
		copy.addAll(edges);
		if(sortByCost) Collections.sort(copy,(e1,e2)->e1.getCost()-e2.getCost());
		else Collections.sort(copy,(e1,e2)->e2.getOverlap()-e1.getOverlap());
		for(AssemblyEdge edge:copy) {
			if(includeEmbedded || !graph.isEmbedded(edge.getConnectingVertex(v).getSequenceIndex())) out.println(edge);
		}
		out.println();
		
	}
	
	private void compareLayouts(AssemblyGraph goldStandardGraph, AssemblyGraph testGraph, PrintStream out) {
		List<AssemblyPath> gsPaths = goldStandardGraph.getPaths();
		List<AssemblyPath> testPaths = testGraph.getPaths();
		errorsTPEdgeNoLayout = 0;
		errorsEdgeEmbeddedNoLayout = 0;
		errorsFPEdge = 0;
		errorsFNLayoutEdge = 0;
		totalGSLayoutEdges = 0;
		totalTestLayoutPaths = 0;
		totalTestLayoutEdges = 0;
		for(AssemblyPath gsPath:gsPaths) {
			totalGSLayoutEdges+=gsPath.getPathLength();
		}
		System.out.println();
		for(int i=0;i<testPaths.size();i++) {
			AssemblyPath nextPath = testPaths.get(i);
			if(nextPath.getPathLength()<=1) continue;
			Map<String,Integer> sequencesPathCounts = new HashMap<String,Integer>();
			long estimatedLength=0;
			int lastOverlap = 0;
			log.info("Compare layouts. Next path: "+(i+1)+" Limits "+nextPath.getVertexLeft()+" to "+nextPath.getVertexRight()+" Edges: "+nextPath.getPathLength());
			totalTestLayoutPaths++;
			totalTestLayoutEdges+=nextPath.getPathLength();
			AssemblyPath nextGSPath = null;
			List<AssemblyEdge> nextGSPathEdges = null;
			int nextGSEdgeIdx = -1;
			int direction = 0;
			List<AssemblyEdge> edges = nextPath.getEdges();
			for(AssemblyEdge nextTestEdge:edges) {
				if(nextTestEdge.isSameSequenceEdge()) {
					String readId = nextTestEdge.getVertex1().getRead().getName();
					int idx = readId.lastIndexOf("_");
					if(idx>0) {
						String s1 = readId.substring(0,idx);
						idx = s1.lastIndexOf("_");
						//if(nextTestEdge.getVertex1().getRead().getLength()>20000) System.out.println("Parsing read id: Read id"+readId+" s1: "+s1+" idx: "+idx);
						if(idx>0) sequencesPathCounts.compute(s1.substring(0, idx),(k,v)->v==null?1:v+1);
					}
					estimatedLength+= nextTestEdge.getVertex1().getRead().getLength()-lastOverlap;
				} else lastOverlap = nextTestEdge.getOverlap();
				boolean searchGSEdge = false;
				if(nextGSPath==null) {
					searchGSEdge = true;
				} else {
					if (direction == 0 && nextGSEdgeIdx<nextGSPathEdges.size()-1 && sameEdges(nextGSPathEdges.get(nextGSEdgeIdx+1), nextTestEdge)) direction = 1;
					else if (direction == 0) direction = -1;
					nextGSEdgeIdx+=direction;
					if(nextGSEdgeIdx>=0 && nextGSEdgeIdx<nextGSPathEdges.size()) {
						AssemblyEdge nextGSEdge = nextGSPathEdges.get(nextGSEdgeIdx);
						if(sameEdges(nextGSEdge, nextTestEdge)) {
							tpLayoutEdges++;
						} else {
							searchGSEdge = true;
						}
					} else {
						System.err.println("Compare layouts. Test path went over gs path. Edge leaving: "+nextTestEdge);
						searchGSEdge = true;
					}
				}
				if(searchGSEdge) {
					direction = 0;
					int [] gsEdgeLocation = findGSEdgeLocation(gsPaths, nextTestEdge );
					if(gsEdgeLocation == null) {
						AssemblyEdge gsEdge = findEdge(goldStandardGraph, nextTestEdge);
						AssemblyEmbedded gsEmbedded = findEmbeddedRelationship(goldStandardGraph, nextTestEdge);
						if(gsEdge==null && gsEmbedded==null) {
							errorsFPEdge++;
							System.err.println("Compare layouts. False positive edge: "+nextTestEdge);
							distOverlapsFPPathEdges.processDatapoint(nextTestEdge.getOverlap());
							distScoresFPPathEdges.processDatapoint(nextTestEdge.getScore());
							distCostsFPPathEdges.processDatapoint(nextTestEdge.getCost());
							distCSKFPPathEdges.processDatapoint(nextTestEdge.getCoverageSharedKmers());
							distWCSKFPPathEdges.processDatapoint(nextTestEdge.getWeightedCoverageSharedKmers());
							distEvidencePropLengthFPPathEdges.processDatapoint(nextTestEdge.getEvidenceProportion());
							distWCSKProportionFPPathEdges.processDatapoint((double)nextTestEdge.getWeightedCoverageSharedKmers()/(nextTestEdge.getOverlap()+1));
							distNumIndelsFPPathEdges.processDatapoint(nextTestEdge.getNumIndels());
							distIndelsKbpFPPathEdges.processDatapoint(nextTestEdge.getIndelsPerKbp());
						} else if (gsEdge==null) { 
							errorsEdgeEmbeddedNoLayout++;
							System.err.println("Compare layouts. Edge from false negative embedded relationship "+nextTestEdge);
						} else  if (goldStandardGraph.isEmbedded(gsEdge.getVertex1().getSequenceIndex()) || goldStandardGraph.isEmbedded(gsEdge.getVertex2().getSequenceIndex())) {
							errorsEdgeEmbeddedNoLayout++;
							System.err.println("Compare layouts. True edge between embedded sequences "+nextTestEdge);
						} else {
							errorsTPEdgeNoLayout++;
							if(!nextTestEdge.isSameSequenceEdge()) System.err.println("Compare layouts. True positive no GS layout "+nextTestEdge);
						}
						if(nextGSPath!=null && nextGSEdgeIdx>=0 && nextGSEdgeIdx<nextGSPathEdges.size()) System.err.println("Last GS layout edge "+nextGSPathEdges.get(nextGSEdgeIdx));
						nextGSPath = null;
						nextGSEdgeIdx = -1;
					} else  {
						if (nextGSPath!=null) {
							log.warning("Found path edge connecting to different gold standard paths. Edge "+nextTestEdge);
						}
						tpLayoutEdges++;
						nextGSPath = gsPaths.get(gsEdgeLocation[0]);
						nextGSPathEdges = new ArrayList<AssemblyEdge>(nextGSPath.getEdges());
						nextGSEdgeIdx = gsEdgeLocation[1];
					}
				}
			}
			if(nextGSEdgeIdx>0 && nextGSEdgeIdx<nextGSPathEdges.size()-1) {
				nextGSEdgeIdx+=direction;
				errorsFNLayoutEdge++;
				AssemblyEdge nextGSEdge = nextGSPathEdges.get(nextGSEdgeIdx);
				log.info("Compare layouts. Finished test path before GS path. last test edge: "+edges.get(edges.size()-1)+"\nNext GS edge after end of test path: "+nextGSEdge);	
			} else if (nextGSEdgeIdx==-1) {
				log.info("Compare layouts. Finished test path without concordance with GS path. last test edge: "+edges.get(edges.size()-1));
			}
			log.info("Compare layouts. Finished path: "+(i+1)+" edges: "+nextPath.getPathLength()+" estimated length: "+estimatedLength+" Sequences "+sequencesPathCounts);
			System.out.println();
		}
	}
	
	private void findProblematicVertices(AssemblyGraph testGraph) {
		System.out.println("Finding problematic vertices.");
		int n = testGraph.getNumSequences();
		for(int i=0;i<n;i++) {
			//Check embedded status
			if(!testGraph.isEmbedded(i)) {
				checkProblematicVertex(testGraph, testGraph.getVertex(i, true));
				checkProblematicVertex(testGraph, testGraph.getVertex(i, false));
			}
		}
	}
	private void checkProblematicVertex(AssemblyGraph testGraph, AssemblyVertex vTest) {
		if(vTest==null) return;
		if(testGraph.isEmbedded(vTest.getSequenceIndex())) return;
		int debugIdx = -1;
		List<AssemblyEdge> edgesTest = new ArrayList<>(testGraph.getEdges(vTest));
		if(edgesTest.size()<3) return;
		Collections.sort(edgesTest, (e1,e2)->e1.getCost()-e2.getCost());
		int d1 = edgesTest.size();
		AssemblyEdge e1 = edgesTest.get(0);
		if(e1.isSameSequenceEdge()) e1 = edgesTest.get(1);
		if(vTest.getSequenceIndex()==debugIdx) System.out.println("E1 with embedded: "+testGraph.isEmbedded(e1.getConnectingVertex(vTest).getSequenceIndex()));
		int compared = 0;
		for(int j=1;j<edgesTest.size() && compared<5;j++) {
			AssemblyEdge e2 = edgesTest.get(j);
			if(vTest.getSequenceIndex()==debugIdx) System.out.println("Number of comparisons: "+compared+" j:"+j+" Comparing edge: "+e1+" with "+e2);
			if(e2==e1 || e2.isSameSequenceEdge()) continue;
			if(vTest.getSequenceIndex()==debugIdx) System.out.println("Pass 2. e2emb: "+testGraph.isEmbedded(e2.getConnectingVertex(vTest).getSequenceIndex()));
			if(testGraph.isEmbedded(e2.getConnectingVertex(vTest).getSequenceIndex())) {
				compared++;
				continue;
			}
			if(testGraph.isEmbedded(e1.getConnectingVertex(vTest).getSequenceIndex())) {
				e1=e2;
				continue;
			}
			//if(vTest.getSequenceIndex()==debugIdx) System.out.println("Pass 3");
			if(e2.getOverlap()<0.5*e1.getOverlap()) break;
			else if (e2.getOverlap()>=e1.getOverlap()) {
				AssemblyEdge eT = e1;
				e1 = e2;
				e2 = eT;
			}
			AssemblyVertex v1 = e1.getConnectingVertex(vTest);
			AssemblyVertex v1C = testGraph.getSameSequenceEdge(v1).getConnectingVertex(v1);
			int remainingLength1 = testGraph.getSequenceLength(v1.getSequenceIndex())-e1.getOverlap();
			AssemblyVertex v2 = e2.getConnectingVertex(vTest);
			AssemblyVertex v2C = testGraph.getSameSequenceEdge(v2).getConnectingVertex(v2);
			int remainingLength2 = testGraph.getSequenceLength(v2.getSequenceIndex())-e2.getOverlap();
			if(remainingLength1+100<remainingLength2) {
				//There should be an edge between v1V and v2
				AssemblyEdge edge = testGraph.getEdge(v1C, v2);
				if(edge==null) System.out.println("Possible edge within repetitive sequence. Vertex: "+vTest+" disconnected vertices: "+v1C+" and "+v2+" remaining lengths: "+remainingLength1+" "+remainingLength2+" degrees "+d1+" "+testGraph.getEdges(v1).size()+" "+testGraph.getEdges(v2).size()+" embedded: "+testGraph.isEmbedded(vTest.getSequenceIndex())+" "+testGraph.isEmbedded(v1.getSequenceIndex())+" "+" "+testGraph.isEmbedded(v2.getSequenceIndex()));
			} else {
				//TODO: v2 should be embedded within v1
				
			}
			compared++;
		}
	}

	private AssemblyEdge findEdge(AssemblyGraph goldStandardGraph, AssemblyEdge nextTestEdge) {
		AssemblyVertex v1 = goldStandardGraph.getVertexByUniqueId(nextTestEdge.getVertex1().getUniqueNumber());
		AssemblyVertex v2 = goldStandardGraph.getVertexByUniqueId(nextTestEdge.getVertex2().getUniqueNumber());
		return goldStandardGraph.getEdge(v1, v2);
	}
	private AssemblyEmbedded findEmbeddedRelationship(AssemblyGraph goldStandardGraph, AssemblyEdge nextTestEdge) {
		int seqId1 = nextTestEdge.getVertex1().getSequenceIndex();
		int seqId2 = nextTestEdge.getVertex2().getSequenceIndex();
		List<AssemblyEmbedded> embeddedList1 = goldStandardGraph.getEmbeddedBySequenceId(seqId1);
		for(AssemblyEmbedded embedded:embeddedList1) {
			if(embedded.getHostId()==seqId2) return embedded;
		}
		List<AssemblyEmbedded> embeddedList2 = goldStandardGraph.getEmbeddedBySequenceId(seqId2);
		for(AssemblyEmbedded embedded:embeddedList2) {
			if(embedded.getHostId()==seqId1) return embedded;
		}
		return null;
	}
	private int[] findGSEdgeLocation(List<AssemblyPath> gsPaths, AssemblyEdge nextTestEdge) {
		for(int i=0;i<gsPaths.size();i++) {
			AssemblyPath nextPath = gsPaths.get(i);
			List<AssemblyEdge> nextPathEdges = new ArrayList<AssemblyEdge>(nextPath.getEdges());
			for(int j=0;j<nextPathEdges.size();j++) {
				AssemblyEdge nextGSEdge = nextPathEdges.get(j);
				if(sameEdges(nextGSEdge, nextTestEdge)) {
					int [] answer = {i,j};
					return answer;
				}
			}
			
		}
		return null;
	}
	private void printStatistics(PrintStream out) {
		double precision = (double)tpEmbSeqs/(tpEmbSeqs+fpEmbSeqs);
		double recall = (double)tpEmbSeqs/(tpEmbSeqs+fnEmbSeqs);
		out.println("EMBEDDED_SEQUENCES\t"+tpEmbSeqs+"\t"+fpEmbSeqs+"\t"+fnEmbSeqs+"\t"+precision+"\t"+recall);
		precision = (double)tpEmbRel/(tpEmbRel+fpEmbRel);
		recall = (double)tpEmbRel/(tpEmbRel+fnEmbRel);
		out.println("EMBEDDED_RELATIONS\t"+tpEmbRel+"\t"+fpEmbRel+"\t"+fnEmbRel+"\t"+precision+"\t"+recall);
		tpEdgesEmbedded/=2;
		tpEdgesNotEmbedded/=2;
		fpEdges/=2;
		fnEdgesEmbedded/=2;
		fnEdgesNotEmbedded/=2;
		tpPathEdges/=2;
		totalPathEdges/=2;
		double precisionEdges = (double)tpEdgesNotEmbedded/(tpEdgesNotEmbedded+fpEdges);
		double recallEdges = (double)tpEdgesNotEmbedded/(tpEdgesNotEmbedded+fnEdgesNotEmbedded);
		double recallPathEdges = (double)tpPathEdges/totalPathEdges;
		out.println("EDGES\t"+tpEdgesNotEmbedded+"\t"+fpEdges+"\t"+fnEdgesNotEmbedded+"\t"+precisionEdges+"\t"+recallEdges+"\t"+tpPathEdges+"\t"+totalPathEdges+"\t"+recallPathEdges+"\t"+tpEdgesEmbedded+"\t"+fnEdgesEmbedded);
		precision = 0;
		if(totalTestLayoutEdges>0) precision = (double)tpLayoutEdges/totalTestLayoutEdges;
		recall = (double)tpLayoutEdges/totalGSLayoutEdges;
		out.println("PATHS\t"+totalTestLayoutPaths+"\t"+tpLayoutEdges+"\t"+errorsFPEdge+"\t"+errorsEdgeEmbeddedNoLayout+"\t"+errorsTPEdgeNoLayout+"\t"+errorsFNLayoutEdge+"\t"+totalTestLayoutEdges+"\t"+totalGSLayoutEdges+"\t"+precision+"\t"+recall);
		out.println();
		
		out.println("Predicted overlap estimation errors");
		out.println("Current: "+Math.sqrt(rmsePredictedOverlap/countPredictedOverlap)+" count "+countPredictedOverlap);
		out.println("Average: "+Math.sqrt(rmseAveragePredictedOverlap/countAveragePredictedOverlap)+" count "+countAveragePredictedOverlap);
		out.println("Median: "+Math.sqrt(rmseMedianPredictedOverlap/countMedianPredictedOverlap)+" count "+countMedianPredictedOverlap);
		out.println("FromLimits: "+Math.sqrt(rmseFromLimitsPredictedOverlap/countFromLimitsPredictedOverlap)+" count "+countFromLimitsPredictedOverlap);
		out.println();
		double [] d1 = distOverlapError.getDistribution();
		double [] d2 = distAverageOverlapError.getDistribution();
		double [] d3 = distMedianOverlapError.getDistribution();
		double [] d4 = distFromLimitsOverlapError.getDistribution();
		out.println("Number\tCurrent\tAverage\tMedian\tFromLimits");
		out.println("Less+\t"+distOverlapError.countOutliersLess()+"\t"+distAverageOverlapError.countOutliersLess()+"\t"+distMedianOverlapError.countOutliersLess()+"\t"+distFromLimitsOverlapError.countOutliersLess());
		for(int i=0;i<d1.length;i++) {
			int min = -500+i*(int)distOverlapError.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]+"\t"+d4[i]);
		}
		out.println("More+\t"+distOverlapError.countOutliersMore()+"\t"+distAverageOverlapError.countOutliersMore()+"\t"+distMedianOverlapError.countOutliersMore()+"\t"+distFromLimitsOverlapError.countOutliersMore());
		
		d1 = distLengthsLayoutReads.getDistribution();
		d2 = distSumLengthsIKbpLayout.getDistribution();
		d3 = distSumLengthsLayout.getDistribution();
		d4 = distLengthsEmbeddedReads.getDistribution();
		double [] d5 = distSumLengthsIKbpEmbedded.getDistribution();
		double [] d6 = distSumLengthsEmbedded.getDistribution();
		out.println("Lengths reads");
		out.println("Number\tLayout\tLayoutIKbp\tEmbedded\tEmbeddedIKbp");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distLengthsLayoutReads.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]/(d3[i]+0.01)+"\t"+d4[i]+"\t"+d5[i]/(d6[i]+0.01));
		}
		
		d1 = distScoresTPEmbedded.getDistribution();
		d2 = distScoresFPEmbedded.getDistribution();
		d3 = distCostsTPEmbedded.getDistribution();
		d4 = distCostsFPEmbedded.getDistribution();
		out.println("Score/Cost distributions embedded");
		out.println("Number\tScoreTP\tScoreFP\tCostTP\tCostFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distCostsTPEmbedded.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]+"\t"+d4[i]);
		}
		out.println("More+\t"+distScoresTPEmbedded.countOutliersMore()+"\t"+distScoresFPEmbedded.countOutliersMore()+"\t"+distCostsTPEmbedded.countOutliersMore()+"\t"+distCostsFPEmbedded.countOutliersMore());
		
		
		d1 = distCSKTPEmbedded.getDistribution();
		d2 = distCSKFPEmbedded.getDistribution();
		d3 = distWCSKTPEmbedded.getDistribution();
		d4 = distWCSKFPEmbedded.getDistribution();
		out.println("CSK embedded");
		out.println("Number\tCSK_TP\tCSK_FP\tWCSK_TP\tWCSK_FP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distCSKTPEmbedded.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]+"\t"+d4[i]);
		}
		
		d1 = distEvidencePropLengthTPEmbedded.getDistribution();
		d2 = distEvidencePropLengthFPEmbedded.getDistribution();
		d3 = distWCSKPropSelfTPEmbedded.getDistribution();
		d4 = distWCSKPropSelfFPEmbedded.getDistribution();
		out.println("Proportions embedded");
		out.println("Number\tPropEvidenceTP\tPropEvidenceFP\tPropWCSKSelfTP\tPropWCSKSelfFP");
		for(int i=0;i<d1.length;i++) {
			double min = distEvidencePropLengthTPEmbedded.getBinLength()*i;
			out.println(ParseUtils.ENGLISHFMT_PROBABILITIES.format(min)+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]+"\t"+d4[i]);
		}
		
		d1 = distIndelsKbpTPEmbedded.getDistribution();
		d2 = distIndelsKbpFPEmbedded.getDistribution();
		out.println("Indels per kbp distributions embedded");
		out.println("Number\tTP\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distIndelsKbpTPEmbedded.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]);
		}
		out.println("More\t"+distIndelsKbpTPEmbedded.countOutliersMore()+"\t"+distIndelsKbpFPEmbedded.countOutliersMore());
		
		d1 = distOverlapsTPPathEdges.getDistribution();
		d2 = distOverlapsFPPathEdges.getDistribution();
		d3 = distOverlapsFNPathEdges.getDistribution();
		d4 = distOverlapsFPEdges.getDistribution();
		out.println("Overlap distributions");
		out.println("Number\tTPpath\tFPpath\tFNPath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distOverlapsTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]+"\t"+d4[i]);
		}
		
		d1 = distScoresTPPathEdges.getDistribution();
		d2 = distScoresFPPathEdges.getDistribution();
		d3 = distScoresFPEdges.getDistribution();
		d4 = distCostsTPPathEdges.getDistribution();
		d5 = distCostsFPPathEdges.getDistribution();
		d6 = distCostsFPEdges.getDistribution();
		out.println("Score/Cost distributions edges");
		out.println("Number\tScoreTPpath\tScoreFPPath\tScoreFP\tCostTPPath\tCostFPPath\tCostFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distCostsTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]+"\t"+d4[i]+"\t"+d5[i]+"\t"+d6[i]);
		}
		out.println("More\t"+distScoresTPPathEdges.countOutliersMore()+"\t"+distScoresFPPathEdges.countOutliersMore()+"\t"+distScoresFPEdges.countOutliersMore()+"\t"+distCostsTPPathEdges.countOutliersMore()+"\t"+distCostsFPPathEdges.countOutliersMore()+"\t"+distCostsFPEdges.countOutliersMore());
		
		d1 = distSharedKmersTPPathEdges.getDistribution();
		d2 = distSharedKmersFPEdges.getDistribution();
		out.println("Shared kmers distributions edges");
		out.println("Number\tTPpath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distSharedKmersTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]);
		}
		
		d1 = distCSKTPPathEdges.getDistribution();
		d2 = distCSKFPPathEdges.getDistribution();
		d3 = distCSKFPEdges.getDistribution();
		out.println("Coverage shared kmers distributions edges");
		out.println("Number\tTPpath\tFPpath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distCSKTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]);
		}
		
		d1 = distWCSKTPPathEdges.getDistribution();
		d2 = distWCSKFPPathEdges.getDistribution();
		d3 = distWCSKFPEdges.getDistribution();
		out.println("Weighted Coverage shared kmers distributions edges");
		out.println("Number\tTPpath\tFPpath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distWCSKTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]);
		}
		
		d1 = distNumIndelsTPPathEdges.getDistribution();
		d2 = distNumIndelsFPPathEdges.getDistribution();
		d3 = distNumIndelsFPEdges.getDistribution();
		out.println("Number of indels distributions edges");
		out.println("Number\tTPpath\tFPpath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distNumIndelsTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]);
		}
		
		d1 = distIndelsKbpTPPathEdges.getDistribution();
		d2 = distIndelsKbpFPPathEdges.getDistribution();
		d3 = distIndelsKbpFPEdges.getDistribution();
		out.println("Indels per kbp distributions edges");
		out.println("Number\tTPpath\tFPpath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distIndelsKbpTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]);
		}
		out.println("More\t"+distIndelsKbpTPPathEdges.countOutliersMore()+"\t"+distIndelsKbpFPPathEdges.countOutliersMore()+"\t"+distIndelsKbpFPEdges.countOutliersMore());
		
		d1 = distOverlapSDTPPathEdges.getDistribution();
		d2 = distOverlapSDFPEdges.getDistribution();
		out.println("Overlap standard deviation distributions edges");
		out.println("Number\tTPpath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distOverlapSDTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]);
		}
		out.println("More\t"+distOverlapSDTPPathEdges.getOutliers().size()+"\t"+distOverlapSDFPEdges.getOutliers().size());
		
		d1 = distEvidencePropLengthTPPathEdges.getDistribution();
		d2 = distEvidencePropLengthFPPathEdges.getDistribution();
		d3 = distEvidencePropLengthFPEdges.getDistribution();
		d4 = distWCSKProportionTPPathEdges.getDistribution();
		d5 = distWCSKProportionFPPathEdges.getDistribution();
		d6 = distWCSKProportionFPEdges.getDistribution();
		out.println("Proportions distributions edges");
		out.println("Number\tPropEvidenceTP\tPropEvidenceFPpath\tPropEvidenceFP\tPropWCSKOverlapTP\tPropWCSKOverlapFPpath\tPropWCSKOverlapFP");
		for(int i=0;i<d1.length;i++) {
			double min = distEvidencePropLengthTPPathEdges.getBinLength()*i;
			out.println(ParseUtils.ENGLISHFMT_PROBABILITIES.format(min)+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]+"\t"+d4[i]+"\t"+d5[i]+"\t"+d6[i]);
		}
		out.print("\t"+distEvidencePropLengthTPPathEdges.getOutliers().size());
		out.print("\t"+distEvidencePropLengthFPPathEdges.getOutliers().size());
		out.print("\t"+distEvidencePropLengthFPEdges.getOutliers().size());
		out.print("\t"+distWCSKProportionTPPathEdges.getOutliers().size());
		out.print("\t"+distWCSKProportionFPPathEdges.getOutliers().size());
		out.println("\t"+distWCSKProportionFPEdges.getOutliers().size());
		
	}
	private void resetStatistics() {
		
		tpEmbSeqs = 0;
		fpEmbSeqs = 0;
		fnEmbSeqs = 0;
		tpEmbRel = 0;
		fpEmbRel = 0;
		fnEmbRel = 0;
		tpEdgesNotEmbedded = 0;
		tpEdgesEmbedded = 0;
		fpEdges = 0;
		fnEdgesNotEmbedded = 0;
		fnEdgesEmbedded = 0;
		tpPathEdges = 0;
		totalPathEdges = 0;
		
		tpLayoutEdges = 0;
		errorsFPEdge = 0;
		errorsEdgeEmbeddedNoLayout = 0;
		errorsTPEdgeNoLayout = 0;
		errorsFNLayoutEdge = 0;
		totalTestLayoutPaths = 0;
		totalTestLayoutEdges = 0;
		totalGSLayoutEdges = 0;
		
		distLengthsLayoutReads.reset();
		distLengthsEmbeddedReads.reset();
		
		rmsePredictedOverlap=0;
		countPredictedOverlap = 0;
		distOverlapError.reset();
		rmseAveragePredictedOverlap=0;
		countAveragePredictedOverlap = 0;
		distAverageOverlapError.reset();
		rmseMedianPredictedOverlap=0;
		countMedianPredictedOverlap = 0;
		distMedianOverlapError.reset();
		rmseFromLimitsPredictedOverlap=0;
		countFromLimitsPredictedOverlap = 0;
		distFromLimitsOverlapError.reset();
		
		distScoresTPEmbedded.reset();
		distScoresFPEmbedded.reset();
		distCostsTPEmbedded.reset();
		distCostsFPEmbedded.reset();
		distCSKTPEmbedded.reset();
		distCSKFPEmbedded.reset();
		distWCSKTPEmbedded.reset();
		distWCSKFPEmbedded.reset();
		distEvidencePropLengthTPEmbedded.reset();
		distEvidencePropLengthFPEmbedded.reset();
		distWCSKPropSelfTPEmbedded.reset();
		distWCSKPropSelfFPEmbedded.reset();
		distIndelsKbpTPEmbedded.reset();
		distIndelsKbpFPEmbedded.reset();
		
		
		distOverlapsTPPathEdges.reset();
		distOverlapSDTPPathEdges.reset();
		distScoresTPPathEdges.reset();
		distCostsTPPathEdges.reset();
		distSharedKmersTPPathEdges.reset();
		distCSKTPPathEdges.reset();
		distWCSKTPPathEdges.reset();
		distEvidencePropLengthTPPathEdges.reset();
		distSharedKmersProportionTPPathEdges.reset();
		distWCSKProportionTPPathEdges.reset();
		distNumIndelsTPPathEdges.reset();
		distIndelsKbpTPPathEdges.reset();
		
		
		distOverlapsFPEdges.reset();
		distOverlapSDFPEdges.reset();
		distScoresFPEdges.reset();
		distCostsFPEdges.reset();
		distSharedKmersFPEdges.reset();
		distCSKFPEdges.reset();
		distWCSKFPEdges.reset();
		distEvidencePropLengthFPEdges.reset();
		distSharedKmersProportionFPEdges.reset();
		distWCSKProportionFPEdges.reset();
		distNumIndelsFPEdges.reset();
		distIndelsKbpFPEdges.reset();
		
		distOverlapsFNPathEdges.reset();
		
		//FP path edges do not need to be reset
		
		
		
	}
}
