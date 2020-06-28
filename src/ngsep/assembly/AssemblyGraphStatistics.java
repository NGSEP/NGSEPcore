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
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
import ngsep.main.CommandsDescriptor;
import ngsep.main.OptionValuesDecoder;
import ngsep.main.ProgressNotifier;
import ngsep.main.io.ParseUtils;
import ngsep.math.Distribution;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class AssemblyGraphStatistics {

	// Constants for default values
	public static final byte READS_FORMAT_FASTQ=KmersExtractor.INPUT_FORMAT_FASTQ;
	public static final byte READS_FORMAT_FASTA=KmersExtractor.INPUT_FORMAT_FASTA;
	
	// Logging and progress
	private Logger log = Logger.getLogger(AssemblyGraphStatistics.class.getName());
	private ProgressNotifier progressNotifier = null;
	
	//Parameters
	private String inputFile = null;
	private String outputFile = null;
	private String readsFile = null;
	private byte readsFormat = READS_FORMAT_FASTQ;
	private ReferenceGenome genome = null;
	private String alignmentsFile = null;
	private boolean simulated = false;
	
	//Statistics
	private Distribution distOverlapsTPPathEdges = new Distribution(0,100000,2000);
	private Distribution distCostsTPPathEdges = new Distribution(0,100000,2000);
	private Distribution distMismatchesTPPathEdges = new Distribution(0,10000,200);
	private Distribution distSharedKmersTPPathEdges = new Distribution(0,25000,500);
	private Distribution distCoverageSharedKmersTPPathEdges = new Distribution(0, 50000, 1000);
	private Distribution distMismatchesProportionTPPathEdges = new Distribution(0, 1, 0.01);
	private Distribution distSharedKmersProportionTPPathEdges = new Distribution(0, 1, 0.01);
	private Distribution distCoverageSharedKmersProportionTPPathEdges = new Distribution(0, 1, 0.01);
	
	
	private Distribution distOverlapsFPEdges = new Distribution(0,100000,2000);
	private Distribution distCostsFPEdges = new Distribution(0,100000,2000);
	private Distribution distMismatchesFPEdges = new Distribution(0,10000,200);
	private Distribution distSharedKmersFPEdges = new Distribution(0,25000,500);
	private Distribution distCoverageSharedKmersFPEdges = new Distribution(0, 50000, 1000);
	private Distribution distMismatchesProportionFPEdges = new Distribution(0, 1, 0.01);
	private Distribution distSharedKmersProportionFPEdges = new Distribution(0, 1, 0.01);
	private Distribution distCoverageSharedKmersProportionFPEdges = new Distribution(0, 1, 0.01);
	
	private Distribution distOverlapsFNPathEdges = new Distribution(0,100000,2000);
	private Distribution distCostsFNPathEdges = new Distribution(0,100000,2000);
	
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
	
	public String getReadsFile() {
		return readsFile;
	}
	public void setReadsFile(String readsFile) {
		this.readsFile = readsFile;
	}
	
	
	public byte getReadsFormat() {
		return readsFormat;
	}
	public void setReadsFormat(byte readsFormat) {
		this.readsFormat = readsFormat;
	}
	public void setReadsFormat(String value) {
		this.setReadsFormat((byte) OptionValuesDecoder.decode(value, Byte.class));
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
		run (inputFile, outputFile);
		log.info("Process finished");
	}
	private void logParameters() {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		out.println("Input file:"+ inputFile);
		out.println("Output file:"+ outputFile);
		if (readsFile!=null) out.println("File with original reads: "+readsFile);
		if(simulated) out.println("Reads were simulated using SingleReadsSimulator");
		if (genome!=null) out.println("Target genome for benchmark loaded from file: "+genome.getFilename());
		if (alignmentsFile!=null) out.println("Alignments file for benchmark "+alignmentsFile);
		log.info(os.toString());
	}
	private void run(String inputFile, String outputFile) throws IOException {
		AssemblyGraph graph = AssemblyGraph.load(inputFile);
		try (PrintStream out=new PrintStream(outputFile)) {
			Distribution degreeDist = graph.getVertexDegreeDistribution ();
			out.println("Vertex degree distribution");
			degreeDist.printDistributionInt(out);
			if(genome==null)  return;
			List<ReadAlignment> alignments = null;
			List<QualifiedSequence> sequences = graph.getSequences();
			if (alignmentsFile!=null) {
				alignments = new ArrayList<ReadAlignment>();
				sequences = new ArrayList<QualifiedSequence>();
				try (ReadAlignmentFileReader reader = new ReadAlignmentFileReader(alignmentsFile)) {
					reader.setFilterFlags(ReadAlignment.FLAG_SECONDARY);
					Iterator<ReadAlignment> it = reader.iterator();
					while (it.hasNext()) {
						ReadAlignment aln = it.next();
						CharSequence characters = aln.getReadCharacters();
						if(aln.isNegativeStrand()) {
							characters = DNAMaskedSequence.getReverseComplement(characters);
						}
						sequences.add(new QualifiedSequence(aln.getReadName(), characters));
						if(!aln.isReadUnmapped()) alignments.add(aln);
					}
				}
				Collections.sort(sequences, (l1, l2) -> l2.getLength() - l1.getLength());
			} else if (simulated) {
				alignments = buildAlignmentsFromSimulatedReads(sequences);
			}
			if(alignments==null) return;
			AssemblyGraph goldStandardGraph = buildGoldStandardGraph(alignments, sequences);
			compareGraphs(goldStandardGraph, graph, out);
			out.println("Initial graph statistics. Vertices: "+graph.getVertices().size()+" edges: "+graph.getNumEdges());
			printStatistics(out);
			resetStatistics();
			graph.removeVerticesChimericReads();
			graph.filterEdgesAndEmbedded();
			log.info("Filtered graph. Vertices: "+graph.getVertices().size()+" edges: "+graph.getEdges().size());
			graph.filterEdgesCloseRelationships();
			log.info("Filtered inconsistent transitive. Vertices: "+graph.getVertices().size()+" edges: "+graph.getEdges().size());
			//LayoutBuilder pathsFinder = new LayoutBuilderGreedyMinCost();
			//LayoutBuilder pathsFinder = new LayoutBuilderGreedyMaxCoverageSharedKmers();
			//LayoutBuilder pathsFinder = new LayoutBuilderGreedyMaxOverlap();
			LayoutBuilder pathsFinder = new LayoutBuilderSkeletonBestReciprocal();
			//LayourBuilder pathsFinder = new LayoutBuilderMetricMSTChristofides();
			//LayourBuilder pathsFinder = new LayoutBuilderModifiedKruskal();
			pathsFinder.findPaths(graph);
			//logErrors=true;
			compareGraphs(goldStandardGraph, graph, out);
			compareLayouts(goldStandardGraph, graph, out);
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
			String [] items = readName.split("_");
			QualifiedSequence seqName = seqNames.get(items[0]);
			int first = Integer.parseInt(items[1]);
			boolean reverse = items[2].charAt(0)=='1';
			int flags = 0;
			if (reverse) flags = ReadAlignment.FLAG_READ_REVERSE_STRAND;
			//System.out.println("Next sequence: "+readName+" first: "+first+" reverse: "+reverse+" flags: "+flags);
			ReadAlignment aln = new ReadAlignment(seqName.getName(), first, first+seq.getLength()-1, seq.getLength(), flags);
			aln.setReadNumber(i);
			aln.setReadCharacters(seq.getCharacters());
			alignments.add(aln);
		}
		return alignments;
	}
	private AssemblyGraph buildGoldStandardGraph(List<ReadAlignment> alignments, List<QualifiedSequence> sequences) {
		AssemblyGraph graph = new AssemblyGraph(sequences);
		//Sort by target genome location to calculate edges efficiently
		GenomicRegionComparator comparator = new GenomicRegionComparator(genome.getSequencesMetadata());
		Collections.sort(alignments, comparator);
		for(int i=0;i<alignments.size();i++) {
			ReadAlignment left = alignments.get(i);
			AssemblyVertex vertexLeft = graph.getVertex(left.getReadNumber(), left.isNegativeStrand());
			for(int j=i+1;j<alignments.size();j++) {
				ReadAlignment right = alignments.get(j);
				int cmp = comparator.compare(right, left);
				if(cmp>1) break;
				AssemblyVertex vertexRight = graph.getVertex(right.getReadNumber(), !right.isNegativeStrand());
				int overlap = left.getLast() - right.getFirst() + 1;
				AssemblyEdge edge = new AssemblyEdge(vertexLeft, vertexRight, overlap);
				graph.addEdge(edge);
				
				boolean relativeNegative = left.isNegativeStrand()!=right.isNegativeStrand();
				int relativeStart;
				
				if(left.getFirst()== right.getFirst() && left.getLast()<=right.getLast()) {
					//left is embedded in right
					if(right.isNegativeStrand()) {
						relativeStart = right.getLast() - left.getLast();
					} else {
						relativeStart = 0;
					}
					AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(left.getReadNumber(), left.getReadCharacters(), relativeNegative, right.getReadNumber(), relativeStart);
					graph.addEmbedded(embeddedEvent);
					if(right.getLast()<=left.getLast()) {
						//right is also embedded in left
						embeddedEvent = new AssemblyEmbedded(right.getReadNumber(), right.getReadCharacters(), relativeNegative, left.getReadNumber(), 0 );
						graph.addEmbedded(embeddedEvent);
					}
				} else if(right.getLast()<=left.getLast()) {
					//Right is embedded in left
					if(left.isNegativeStrand()) {
						relativeStart = left.getLast() - right.getLast();
					} else {
						relativeStart = right.getFirst()-left.getFirst();
					}
					AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(right.getReadNumber(), right.getReadCharacters(), relativeNegative, left.getReadNumber(), relativeStart );
					graph.addEmbedded(embeddedEvent);
				}
			}
		}
		log.info("Created gold standard assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges. Embedded: "+graph.getEmbeddedCount());
		//Build gold standard layouts it must be done after knowing which sequences are embedded
		String lastSeqName = null;
		List<AssemblyEdge> nextPath = new ArrayList<>();
		AssemblyVertex lastVertex = null;
		for(int i=0;i<alignments.size();i++) {
			ReadAlignment aln = alignments.get(i);
			if(graph.isEmbedded(aln.getReadNumber())) continue;
			if(!aln.getSequenceName().equals(lastSeqName)) {
				if(nextPath.size()>0) graph.addPath(nextPath);
				lastSeqName = aln.getSequenceName();
				nextPath = new ArrayList<>();
				lastVertex = null;
			}
			AssemblyVertex leftSequenceVertex = graph.getVertex(aln.getReadNumber(), !aln.isNegativeStrand());
			AssemblyVertex rightSequenceVertex = graph.getVertex(aln.getReadNumber(), aln.isNegativeStrand());
			AssemblyEdge edgeSequence = graph.getSameSequenceEdge(leftSequenceVertex);
			if(lastVertex!=null) {
				AssemblyEdge connectingEdge = graph.getEdge(lastVertex, leftSequenceVertex);
				if(connectingEdge==null) {
					log.info("Discontiguity in gold standard layout");
					if(nextPath.size()>0) graph.addPath(nextPath);
					nextPath = new ArrayList<>();
					lastVertex = null;
				} else {
					//System.out.println("Next layout edge between "+lastVertex.getSequenceIndex()+" start " +lastVertex.isStart()+" and "+leftSequenceVertex.getSequenceIndex()+" start "+leftSequenceVertex.isStart()+" read id: "+aln.getSequenceName()+" first: "+aln.getFirst());
					nextPath.add(connectingEdge);
					connectingEdge.setLayoutEdge(true);
				}
			}
			nextPath.add(edgeSequence);
			lastVertex = rightSequenceVertex;
		}
		if(nextPath.size()>0) graph.addPath(nextPath);
		return graph;
	}

	private void compareGraphs(AssemblyGraph goldStandardGraph, AssemblyGraph testGraph, PrintStream out) {
		
		int n = goldStandardGraph.getNumSequences();
		for(int i=0;i<n;i++) {
			QualifiedSequence sequence = goldStandardGraph.getSequence(i);
			//Check embedded status
			boolean gsE = goldStandardGraph.isEmbedded(i);
			boolean testE = testGraph.isEmbedded(i);
			if(gsE && testE) tpEmbSeqs++;
			else if (gsE) fnEmbSeqs++;
			else if (testE) {
				fpEmbSeqs++;
				List<AssemblyEmbedded> falseHosts = testGraph.getEmbeddedBySequenceId(i);
				if (logErrors) log.info("False embedded sequence "+logSequence(i, sequence)+" false hosts: "+falseHosts.size());
				for(AssemblyEmbedded embedded:falseHosts) {
					QualifiedSequence seqHost = testGraph.getSequence(embedded.getHostId()) ;
					if (logErrors) log.info("Next false host "+logSequence(embedded.getHostId(),seqHost));
				}
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
			log.warning("Test vertex not found for start of sequence "+logSequence(seqIndex, sequence));
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
	private String logVertex(AssemblyVertex vertex) {
		return ""+vertex.getUniqueNumber()+" "+logSequence(vertex.getSequenceIndex(), vertex.getRead());
	}
	private String logEdge(AssemblyEdge edge) {
		return "v1 "+logVertex(edge.getVertex1())+" v2: "+logVertex(edge.getVertex2())+" overlap: "+edge.getOverlap()+" mismatches: "+edge.getMismatches()+" coverage shared kmers: "+edge.getCoverageSharedKmers()+" cost: "+edge.getCost();
	}
	private void calculateComparisonStats(AssemblyGraph goldStandardGraph, AssemblyVertex gsVertex, AssemblyGraph testGraph,  AssemblyVertex testVertex, PrintStream out) {
		//Find path edge of this vertex
		List<AssemblyEdge> gsEdges = goldStandardGraph.getEdges(gsVertex);
		List<AssemblyEdge> testEdges = testGraph.getEdges(testVertex);
		boolean debug = gsVertex.getUniqueNumber()==-25800 || gsVertex.getUniqueNumber()==19204; 
		if(debug) {
			printEdgeList("Gold standard", gsVertex, gsEdges, out);
			printEdgeList("Test", testVertex, testEdges, out);
		}
		Map<Integer,Boolean> testEdgesMatched = new HashMap<Integer, Boolean>();
		for(AssemblyEdge edge:testEdges) {
			testEdgesMatched.put(edge.getConnectingVertex(testVertex).getUniqueNumber(), false);
		}
		boolean gsEmbedded = goldStandardGraph.isEmbedded(gsVertex.getSequenceIndex());
		for(AssemblyEdge gsEdge:gsEdges) {
			AssemblyVertex gsConnectingVertex = gsEdge.getConnectingVertex(gsVertex);
			boolean edgeEmbedded = gsEmbedded || goldStandardGraph.isEmbedded(gsConnectingVertex.getSequenceIndex());
			int number = gsConnectingVertex.getUniqueNumber();
			boolean match = testEdgesMatched.containsKey(number);
			if(match) {
				//True positive
				if(edgeEmbedded) tpEdgesEmbedded++;
				else tpEdgesNotEmbedded++;
				testEdgesMatched.put(number, true);
			}
			else {
				//False negative
				if(edgeEmbedded) fnEdgesEmbedded++;
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
					AssemblyEdge minCostTestEdge = null;
					AssemblyEdge maxOverlapTestEdge = null;
					AssemblyEdge layoutTestEdge = null;
					for(AssemblyEdge edge:testEdges) {
						if (edge.getConnectingVertex(testVertex).getUniqueNumber() == number) {
							layoutTestEdge = edge;
							distOverlapsTPPathEdges.processDatapoint(edge.getOverlap());
							distCostsTPPathEdges.processDatapoint(calculateCost(edge));
							distMismatchesTPPathEdges.processDatapoint(edge.getMismatches());
							distSharedKmersTPPathEdges.processDatapoint(edge.getEvidence().getNumDifferentKmers());
							distCoverageSharedKmersTPPathEdges.processDatapoint(edge.getCoverageSharedKmers());
							distMismatchesProportionTPPathEdges.processDatapoint((double)edge.getMismatches()/edge.getOverlap());
							distSharedKmersProportionTPPathEdges.processDatapoint((double)edge.getEvidence().getNumDifferentKmers()/edge.getOverlap());
							distCoverageSharedKmersProportionTPPathEdges.processDatapoint((double)edge.getCoverageSharedKmers()/edge.getOverlap());
						}
						if(!edge.isSameSequenceEdge() && (minCostTestEdge==null || minCostTestEdge.getCost()>edge.getCost())) minCostTestEdge=edge;
						if(!edge.isSameSequenceEdge() && (maxOverlapTestEdge==null || maxOverlapTestEdge.getOverlap()<edge.getOverlap())) maxOverlapTestEdge=edge;
					}
					//if(logErrors && minCostTestEdge!=layoutTestEdge) log.info("Min cost edge for vertex "+logVertex(testVertex)+" not in layout. Layout edge: "+logEdge(layoutTestEdge)+" min cost edge: "+logEdge(minCostTestEdge)+" layout embedded: "+goldStandardGraph.isEmbedded(gsConnectingVertex.getSequenceIndex())+" "+testGraph.isEmbedded(gsConnectingVertex.getSequenceIndex())+" min embedded "+goldStandardGraph.isEmbedded(minCostTestEdge.getConnectingVertex(testVertex).getSequenceIndex())+" "+testGraph.isEmbedded(minCostTestEdge.getConnectingVertex(testVertex).getSequenceIndex()));
					if(logErrors && maxOverlapTestEdge!=layoutTestEdge) log.info("Max overlap edge for vertex "+logVertex(testVertex)+" not in layout. Layout edge: "+logEdge(layoutTestEdge)+" ax overlap edge: "+logEdge(maxOverlapTestEdge)+" layout embedded: "+goldStandardGraph.isEmbedded(gsConnectingVertex.getSequenceIndex())+" "+testGraph.isEmbedded(gsConnectingVertex.getSequenceIndex())+" min embedded "+goldStandardGraph.isEmbedded(maxOverlapTestEdge.getConnectingVertex(testVertex).getSequenceIndex())+" "+testGraph.isEmbedded(maxOverlapTestEdge.getConnectingVertex(testVertex).getSequenceIndex()));
					
				} else {
					distOverlapsFNPathEdges.processDatapoint(gsEdge.getOverlap());
					distCostsFNPathEdges.processDatapoint(calculateCost(gsEdge));
					//log.info("Path edge not found between "+logVertex(gsVertex)+ " and "+logVertex(gsConnectingVertex)+" gsEdge: "+logEdge(gsEdge));
				}
				totalPathEdges++;
			}
		}
		for(AssemblyEdge edge:testEdges) {
			AssemblyVertex vertex = edge.getConnectingVertex(testVertex); 
			int number = vertex.getUniqueNumber();
			if(!testEdgesMatched.get(number)) {
				//False positive
				fpEdges++;
				distOverlapsFPEdges.processDatapoint(edge.getOverlap());
				distCostsFPEdges.processDatapoint(calculateCost(edge));
				distMismatchesFPEdges.processDatapoint(edge.getMismatches());
				distSharedKmersFPEdges.processDatapoint(edge.getEvidence().getNumDifferentKmers());
				distCoverageSharedKmersFPEdges.processDatapoint(edge.getCoverageSharedKmers());
				distMismatchesProportionFPEdges.processDatapoint((double)edge.getMismatches()/edge.getOverlap());
				distSharedKmersProportionFPEdges.processDatapoint((double)edge.getEvidence().getNumDifferentKmers()/edge.getOverlap());
				distCoverageSharedKmersProportionFPEdges.processDatapoint((double)edge.getCoverageSharedKmers()/edge.getOverlap());
				
				if (logErrors) log.info("False positive edge "+logEdge(edge));
			}
		}
	}
	private double calculateCost(AssemblyEdge edge) {
		if(edge.isSameSequenceEdge()) return edge.getCost();
		//int cost = edge.getCost();
		int cost = edge.getVertex1().getRead().getLength();
		cost+= edge.getVertex2().getRead().getLength();
		cost-= edge.getCoverageSharedKmers();
		return cost;
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

	public void printEdgeList(String text, AssemblyVertex v, List<AssemblyEdge> edges, PrintStream out) {
		out.println(text+" vertex "+logVertex(v));
		for(AssemblyEdge edge:edges) {
			out.println(logEdge(edge));
		}
		
	}
	
	private void compareLayouts(AssemblyGraph goldStandardGraph, AssemblyGraph testGraph, PrintStream out) {
		List<List<AssemblyEdge>> gsPaths = goldStandardGraph.getPaths();
		List<List<AssemblyEdge>> testPaths = testGraph.getPaths();
		errorsTPEdgeNoLayout = 0;
		errorsEdgeEmbeddedNoLayout = 0;
		errorsFPEdge = 0;
		errorsFNLayoutEdge = 0;
		totalGSLayoutEdges = 0;
		totalTestLayoutPaths = 0;
		totalTestLayoutEdges = 0;
		for(List<AssemblyEdge> gsPath:gsPaths) {
			totalGSLayoutEdges+=gsPath.size();
		}
		
		for(int i=0;i<testPaths.size();i++) {
			List<AssemblyEdge> nextPath = testPaths.get(i);
			if(nextPath.size()<=1) continue;
			log.info("Compare layouts. Limits next path. "+logEdge(nextPath.get(0))+" to "+logEdge(nextPath.get(nextPath.size()-1)));
			totalTestLayoutPaths++;
			totalTestLayoutEdges+=nextPath.size();
			List<AssemblyEdge> nextGSPath = null;
			int nextGSEdgeIdx = -1;
			int direction = 0;
			for(int j=0;j<nextPath.size();j++) {
				AssemblyEdge nextTestEdge = nextPath.get(j);
				boolean searchGSEdge = false;
				if(nextGSPath==null) {
					searchGSEdge = true;
				} else {
					if (direction == 0 && nextGSEdgeIdx<nextGSPath.size()-1 && sameEdges(nextGSPath.get(nextGSEdgeIdx+1), nextTestEdge)) direction = 1;
					else if (direction == 0) direction = -1;
					nextGSEdgeIdx+=direction;
					if(nextGSEdgeIdx>=0 && nextGSEdgeIdx<nextGSPath.size()) {
						AssemblyEdge nextGSEdge = nextGSPath.get(nextGSEdgeIdx);
						if(sameEdges(nextGSEdge, nextTestEdge)) {
							tpLayoutEdges++;
						} else {
							searchGSEdge = true;
						}
					} else {
						log.info("Compare layouts. Test path went over gs path. Edge leaving: "+logEdge(nextTestEdge));
						searchGSEdge = true;
					}
				}
				if(searchGSEdge) {
					direction = 0;
					int [] gsEdgeLocation = findGSEdgeLocation(gsPaths, nextTestEdge );
					if(gsEdgeLocation == null) {
						AssemblyEdge gsEdge = findEdge(goldStandardGraph, nextTestEdge);
						if(gsEdge==null) {
							errorsFPEdge++;
							log.info("Compare layouts. False positive edge: "+logEdge(nextTestEdge));
						}
						else if (goldStandardGraph.isEmbedded(gsEdge.getVertex1().getSequenceIndex()) || goldStandardGraph.isEmbedded(gsEdge.getVertex2().getSequenceIndex())) {
							errorsEdgeEmbeddedNoLayout++;
							log.info("Compare layouts. Embedded no GS layout "+logEdge(nextTestEdge));
						}
						else {
							errorsTPEdgeNoLayout++;
							log.info("Compare layouts. True positive no GS layout "+logEdge(nextTestEdge));
						}
						if(nextGSPath!=null && nextGSEdgeIdx>=0 && nextGSEdgeIdx<nextGSPath.size()) log.info("Last GS layout edge "+logEdge(nextGSPath.get(nextGSEdgeIdx)));
						nextGSPath = null;
						nextGSEdgeIdx = -1;
					} else  {
						if (nextGSPath!=null) {
							log.warning("Found path edge connecting to different gold standard paths. Edge "+logEdge(nextTestEdge));
						}
						tpLayoutEdges++;
						nextGSPath = gsPaths.get(gsEdgeLocation[0]);
						nextGSEdgeIdx = gsEdgeLocation[1];
					}
				}
			}
			if(nextGSEdgeIdx>0 && nextGSEdgeIdx<nextGSPath.size()-1) {
				nextGSEdgeIdx+=direction;
				errorsFNLayoutEdge++;
				AssemblyEdge nextGSEdge = nextGSPath.get(nextGSEdgeIdx);
				log.info("Compare layouts. Finished test path before GS path. last test edge: "+logEdge(nextPath.get(nextPath.size()-1))+" Next GS edge after end of test path: "+logEdge(nextGSEdge));	
			} else if (nextGSEdgeIdx==-1) {
				log.info("Compare layouts. Finished test path without concordance with GS path. last test edge: "+logEdge(nextPath.get(nextPath.size()-1)));
			}
		}
	}

	private AssemblyEdge findEdge(AssemblyGraph goldStandardGraph, AssemblyEdge nextTestEdge) {
		AssemblyVertex v1 = goldStandardGraph.getVertexByUniqueId(nextTestEdge.getVertex1().getUniqueNumber());
		AssemblyVertex v2 = goldStandardGraph.getVertexByUniqueId(nextTestEdge.getVertex2().getUniqueNumber());
		return goldStandardGraph.getEdge(v1, v2);
	}
	private int[] findGSEdgeLocation(List<List<AssemblyEdge>> gsPaths, AssemblyEdge nextTestEdge) {
		for(int i=0;i<gsPaths.size();i++) {
			List<AssemblyEdge> nextPath = gsPaths.get(i);
			for(int j=0;j<nextPath.size();j++) {
				AssemblyEdge nextGSEdge = nextPath.get(j);
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
		double [] d1 = distOverlapsTPPathEdges.getDistribution();
		double [] d2 = distCostsFPEdges.getDistribution();
		double [] d3 = distOverlapsFNPathEdges.getDistribution();
		
		
		out.println("Overlap distributions");
		out.println("Number\tTPpath\tFP\tFNPath");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distOverlapsTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]);
		}
		d1 = distCostsTPPathEdges.getDistribution();
		d2 = distCostsFPEdges.getDistribution();
		d3 = distCostsFNPathEdges.getDistribution();
		
		
		out.println("Cost distributions");
		out.println("Number\tTPpath\tFP\tFNPath");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distCostsTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]);
		}
		
		d1 = distMismatchesTPPathEdges.getDistribution();
		d2 = distMismatchesFPEdges.getDistribution();
		out.println("Mismatches distributions");
		out.println("Number\tTPpath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distMismatchesTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]);
		}
		d1 = distSharedKmersTPPathEdges.getDistribution();
		d2 = distSharedKmersFPEdges.getDistribution();
		out.println("Shared kmers distributions");
		out.println("Number\tTPpath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distSharedKmersTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]);
		}
		
		d1 = distCoverageSharedKmersTPPathEdges.getDistribution();
		d2 = distCoverageSharedKmersFPEdges.getDistribution();
		out.println("Coverage shared kmers distributions");
		out.println("Number\tTPpath\tFP");
		for(int i=0;i<d1.length;i++) {
			int min = i*(int)distCoverageSharedKmersTPPathEdges.getBinLength();
			out.println(min+"\t"+d1[i]+"\t"+d2[i]);
		}
		
		d1 = distMismatchesProportionTPPathEdges.getDistribution();
		d2 = distMismatchesProportionFPEdges.getDistribution();
		d3 = distSharedKmersProportionTPPathEdges.getDistribution();
		double [] d4 = distSharedKmersProportionFPEdges.getDistribution();
		double [] d5 = distCoverageSharedKmersProportionTPPathEdges.getDistribution();
		double [] d6 = distCoverageSharedKmersProportionFPEdges.getDistribution();
		out.println("Proportions distributions");
		out.println("Number\tMismatchesTPpath\tMismatchesFP\tsharedKmersTPpath\tsharedKmersFP\tCovSharedKmersTPpath\tCovSharedKmersFP");
		for(int i=0;i<d1.length;i++) {
			double min = distMismatchesProportionTPPathEdges.getBinLength()*i;
			out.println(ParseUtils.ENGLISHFMT_PROBABILITIES.format(min)+"\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]+"\t"+d4[i]+"\t"+d5[i]+"\t"+d6[i]);
		}
		//out.println("More\t"+d1[i]+"\t"+d2[i]+"\t"+d3[i]+"\t"+d4[i]+"\t"+d5[i]+"\t"+d6[i]);
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
		distOverlapsTPPathEdges.reset();
		distCostsTPPathEdges.reset();
		distMismatchesTPPathEdges.reset();
		distSharedKmersTPPathEdges.reset();
		distCoverageSharedKmersTPPathEdges.reset();
		distMismatchesProportionTPPathEdges.reset();
		distSharedKmersProportionTPPathEdges.reset();
		distCoverageSharedKmersProportionTPPathEdges.reset();
		
		
		distOverlapsFPEdges.reset();
		distCostsFPEdges.reset();
		distMismatchesFPEdges.reset();
		distSharedKmersFPEdges.reset();
		distCoverageSharedKmersFPEdges.reset();
		distMismatchesProportionFPEdges.reset();
		distSharedKmersProportionFPEdges.reset();
		distCoverageSharedKmersProportionFPEdges.reset();
		
		distOverlapsFNPathEdges.reset();
		distCostsFNPathEdges.reset();
		
		
	}
}
