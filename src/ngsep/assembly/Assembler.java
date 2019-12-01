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
import ngsep.alignments.ReadAlignment;
import ngsep.genome.GenomicRegionComparator;
import ngsep.genome.ReferenceGenome;
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
	
	
	public static final byte INPUT_FORMAT_FASTQ=0;
	public static final byte INPUT_FORMAT_FASTA=1;
	public static final byte INPUT_FORMAT_GRAPH=2;
	public static final int DEF_KMER_SIZE = 15;
	public static final int DEF_KMER_OFFSET = 15;
	public static final int DEF_MIN_KMER_PCT = 40;
	
	private byte inputFormat = INPUT_FORMAT_FASTQ;
	
	private String outFileGraph = null;
	private int kmerLength = DEF_KMER_SIZE;
	private int kmerOffset = DEF_KMER_OFFSET;
	private int minKmerPercentage = DEF_MIN_KMER_PCT;
	private ReferenceGenome targetGenome;
	
	
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
	public byte getInputFormat() {
		return inputFormat;
	}

	/**
	 * @param inputFormat the inputFormat to set
	 */
	public void setInputFormat(byte inputFormat) {
		this.inputFormat = inputFormat;
	}

	public void setInputFormat(String value) {
		this.setInputFormat((byte) OptionValuesDecoder.decode(value, Byte.class));
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
	
	/**
	 * @return the targetGenome
	 */
	public ReferenceGenome getTargetGenome() {
		return targetGenome;
	}

	/**
	 * @param targetGenome the targetGenome to set
	 */
	public void setTargetGenome(ReferenceGenome targetGenome) {
		this.targetGenome = targetGenome;
	}

	public void run(String inputFile, String outputFile) throws IOException {
		AssemblyGraph graph;
		AssemblyGraph goldStandardGraph=null;
		if(INPUT_FORMAT_GRAPH==inputFormat) {
			graph = AssemblyGraph.load(inputFile);
			log.info("Loaded assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges");
		} else {
			List<QualifiedSequence> sequencesQL = load(inputFile);
			log.info("Loaded "+sequencesQL.size()+" sequences");
			Collections.sort(sequencesQL, (l1, l2) -> l2.getLength() - l1.getLength());
			log.info("Sorted "+sequencesQL.size()+" sequences");
			
			List<CharSequence> finalSequences = Collections.unmodifiableList(extractReads(sequencesQL));
			if(targetGenome!=null) {
				goldStandardGraph = buildGoldStandardGraph(sequencesQL, finalSequences);
			}
			
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
		if(goldStandardGraph!=null) {
			compareGraphs(goldStandardGraph, graph);
		}

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
	public List<QualifiedSequence> load(String filename) throws IOException {
		if (INPUT_FORMAT_FASTQ == inputFormat) return loadFastq(filename);
		else if (INPUT_FORMAT_FASTA==inputFormat) return loadFasta(filename);
		else throw new IOException("the file not is a fasta or fastq file: " + filename);
	}

	/**
	 * Load the sequences of the Fasta file
	 * @param filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private List<QualifiedSequence> loadFasta(String filename) throws IOException {
		FastaSequencesHandler handler = new FastaSequencesHandler();
		QualifiedSequenceList seqsQL = handler.loadSequences(filename);
		List<QualifiedSequence> answer = new ArrayList<>();
		answer.addAll(seqsQL);
		return answer;
	}

	/**
	 * Load the sequences of the Fastq file
	 * 
	 * @param Filename the file path
	 * @return The sequences
	 * @throws IOException The file cannot opened
	 */
	private List<QualifiedSequence> loadFastq(String filename) throws IOException {
		List<QualifiedSequence> sequences = new ArrayList<>();
		try (FastqFileReader reader = new FastqFileReader(filename)) {
			if(targetGenome==null) reader.setLoadMode(FastqFileReader.LOAD_MODE_MINIMAL);
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
	
	private List<CharSequence> extractReads (List<QualifiedSequence> seqsQl) {
		List<CharSequence> sequences = new ArrayList<>();
		for (QualifiedSequence seq : seqsQl) {
			DNAMaskedSequence characters = (DNAMaskedSequence) seq.getCharacters();
			sequences.add(characters);
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
	
	private AssemblyGraph buildGoldStandardGraph(List<QualifiedSequence> sequencesQL, List<CharSequence> finalSequences) {
		//Create true alignments of simulated reads to te target genome
		List<ReadAlignment> alignments = new ArrayList<>();
		QualifiedSequenceList seqNames = targetGenome.getSequencesMetadata();
		for(int i=0;i<sequencesQL.size();i++) {
			QualifiedSequence seq = sequencesQL.get(i);
			String readName = seq.getName();
			String [] items = readName.split("_");
			QualifiedSequence seqName = seqNames.get(items[0]);
			int first = Integer.parseInt(items[1]);
			boolean reverse = items[2].charAt(0)=='1';
			int flags = 0;
			if (reverse) flags = ReadAlignment.FLAG_READ_REVERSE_STRAND;
			//System.out.println("Next sequence: "+readName+" first: "+first+" reverse: "+reverse+" flags: "+flags);
			ReadAlignment aln = new ReadAlignment(seqName.getName(), first, first+seq.getLength()-1, seq.getLength(), flags);
			aln.setReadNumber(i);
			aln.setReadCharacters(finalSequences.get(i));
			alignments.add(aln);
		}
		//Sort by target genome location to calculate edges efficiently 
		AssemblyGraph graph = new AssemblyGraph(finalSequences);
		GenomicRegionComparator comparator = new GenomicRegionComparator(seqNames);
		Collections.sort(alignments, comparator);
		for(int i=0;i<alignments.size();i++) {
			ReadAlignment left = alignments.get(i);
			AssemblyVertex vertexLeft = graph.getVertex(left.getReadNumber(), left.isNegativeStrand());
			
			for(int j=i+1;j<alignments.size();j++) {
				ReadAlignment right = alignments.get(j);
				int cmp = comparator.compare(right, left);
				if(cmp>1) break;
				boolean relativeNegative = left.isNegativeStrand()!=right.isNegativeStrand();
				int relativeStart;
				
				if(left.getFirst()== right.getFirst() && left.getLast()<=right.getLast()) {
					//left is embedded in right
					if(right.isNegativeStrand()) {
						relativeStart = right.getLast() - left.getLast();
					} else {
						relativeStart = 0;
					}
					AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(left.getReadNumber(), left.getReadCharacters(), relativeStart, relativeNegative );
					graph.addEmbedded(right.getReadNumber(), embeddedEvent);
					if(left.getLast()<=right.getLast()) {
						//right is also embedded in left
						embeddedEvent = new AssemblyEmbedded(right.getReadNumber(), right.getReadCharacters(), 0, relativeNegative );
						graph.addEmbedded(left.getReadNumber(), embeddedEvent);
					}
					break;
				}
				if(right.getLast()<=left.getLast()) {
					//Right is embedded in left
					if(left.isNegativeStrand()) {
						relativeStart = left.getLast() - right.getLast();
					} else {
						relativeStart = right.getFirst()-left.getFirst();
					}
					AssemblyEmbedded embeddedEvent = new AssemblyEmbedded(right.getReadNumber(), right.getReadCharacters(), relativeStart, relativeNegative );
					graph.addEmbedded(left.getReadNumber(), embeddedEvent);
				} else {
					AssemblyVertex vertexRight = graph.getVertex(right.getReadNumber(), !right.isNegativeStrand());
					int overlap = left.getLast() - right.getFirst() + 1;
					graph.addEdge(vertexLeft, vertexRight, left.getReadLength()+right.getReadLength() - overlap, overlap);
				}
				
			}
		}
		log.info("Created gold standard assembly graph with "+graph.getVertices().size()+" vertices and "+graph.getEdges().size()+" edges. Embedded: "+graph.getEmbeddedCount());
		graph.pruneEmbeddedSequences();
		log.info("Prunned graph. Edges: "+graph.getEdges().size());
		return graph;
	}

	private void compareGraphs(AssemblyGraph goldStandardGraph, AssemblyGraph testGraph) {
		//Compare embedded
		int tpEmbSeqs = 0;
		int fpEmbSeqs = 0;
		int fnEmbSeqs = 0;
		int tpEmbRel = 0;
		int fpEmbRel = 0;
		int fnEmbRel = 0;
		int tpEdges = 0;
		int fpEdges = 0;
		int fnEdges = 0;
		int n = goldStandardGraph.getNumSequences();
		for(int i=0;i<n;i++) {
			boolean gsE = goldStandardGraph.isEmbedded(i);
			boolean testE = testGraph.isEmbedded(i);
			if(gsE && testE) tpEmbSeqs++;
			else if (gsE) fnEmbSeqs++;
			else if (testE) fpEmbSeqs++;
			List<AssemblyEmbedded> embeddedGS = goldStandardGraph.getEmbedded(i);
			List<AssemblyEmbedded> embeddedTest = testGraph.getEmbedded(i);
			int tpM = calculateIntersection (embeddedGS,embeddedTest);
			tpEmbRel+=tpM;
			fpEmbRel+=(embeddedTest.size()-tpM);
			fnEmbRel+=(embeddedGS.size()-tpM);
			if(goldStandardGraph.isEmbedded(i)) continue;
			AssemblyVertex vGS = goldStandardGraph.getVertex(i, true);
			List<AssemblyEdge> edgesGS = goldStandardGraph.getEdges(vGS);
			int tpD1 = searchEdges(testGraph, vGS, edgesGS);
			tpEdges += tpD1;
			//The same sequence edge does not count
			int fnD1 = edgesGS.size()-1-tpD1;
			if(fnD1>=0) {
				fnEdges += fnD1;
			} else {
				log.warning("Less vertices than true positives for "+i+" with length: "+vGS.getRead().length()+" gsEdges: "+edgesGS.size()+" tp "+tpD1);
			}
			
			AssemblyVertex vTest = testGraph.getVertex(i, true);
			if(vTest!=null) {
				List<AssemblyEdge> testEdges = testGraph.getEdges(vTest);
				int fpD1 = testEdges.size()-1-tpD1; 
				if(fpD1>=0) {
					/*if(fpD1>0) {
						System.out.println("Printing edges comparison for sequence: "+i+" with length: "+vTest.getRead().length()+" start: "+vTest.isStart()+" testEdges: "+testEdges.size()+" gsEdges: "+edgesGS.size()+" tp "+tpD1);
						printEdgeList("GoldStandard",vGS, edgesGS);
						printEdgeList("Test",vTest, testEdges);
					}*/
					fpEdges += fpD1;
				} else {
					log.warning("Less vertices than true positives for "+i+" with length: "+vTest.getRead().length()+" testEdges: "+testEdges.size()+" tp "+tpD1);
				}
				
			}
			vGS = goldStandardGraph.getVertex(i, false);
			edgesGS = goldStandardGraph.getEdges(vGS);
			int tpD2 = searchEdges(testGraph, vGS, edgesGS);
			tpEdges += tpD2;
			int fnD2 = edgesGS.size()-1-tpD2;
			if(fnD2>=0) {
				fnEdges += fnD2;
			} else {
				log.warning("Less vertices than true positives for "+i+" with length: "+vGS.getRead().length()+" gsEdges: "+edgesGS.size()+" tp "+tpD2);
			}
			
			vTest = testGraph.getVertex(i, false);
			if(vTest!=null) {
				List<AssemblyEdge> testEdges = testGraph.getEdges(vTest);
				int fpD2 = testEdges.size()-1-tpD2; 
				if(fpD2>=0) {
					fpEdges += fpD2;
				} else {
					log.warning("Less vertices than true positives for "+i+" with length: "+vTest.getRead().length()+" testEdges: "+testEdges.size()+" tp "+tpD2);
				}
			}
		}
		double precision = (double)tpEmbSeqs/(tpEmbSeqs+fpEmbSeqs);
		double recall = (double)tpEmbSeqs/(tpEmbSeqs+fnEmbSeqs);
		log.info("EMBEDDED_SEQUENCES\t"+tpEmbSeqs+"\t"+fpEmbSeqs+"\t"+fnEmbSeqs+"\t"+precision+"\t"+recall);
		precision = (double)tpEmbRel/(tpEmbRel+fpEmbRel);
		recall = (double)tpEmbRel/(tpEmbRel+fnEmbRel);
		log.info("EMBEDDED_RELATIONS\t"+tpEmbRel+"\t"+fpEmbRel+"\t"+fnEmbRel+"\t"+precision+"\t"+recall);
		tpEdges/=2;
		fpEdges/=2;
		fnEdges/=2;
		double precisionEdges = (double)tpEdges/(tpEdges+fpEdges);
		double recallEdges = (double)tpEdges/(tpEdges+fnEdges);
		log.info("EDGES\t"+tpEdges+"\t"+fpEdges+"\t"+fnEdges+"\t"+precisionEdges+"\t"+recallEdges);	
	}

	public void printEdgeList(String text, AssemblyVertex v, List<AssemblyEdge> edges) {
		System.out.println(text);
		for(AssemblyEdge edge:edges) {
			AssemblyVertex vOut = edge.getConnectingVertex(v);
			System.out.println("Vertex "+vOut.getIndex()+"-"+vOut.isStart()+" length "+vOut.getRead().length()+" overlap "+edge.getOverlap());
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

	private int searchEdges(AssemblyGraph testGraph, AssemblyVertex vGS, List<AssemblyEdge> edgesGS) {
		int count = 0;
		AssemblyVertex testV1 = testGraph.getVertex(vGS.getIndex(), vGS.isStart());
		if(testV1==null) return count;
		List<AssemblyEdge> testEdgesV1 = testGraph.getEdges(testV1);
		for(AssemblyEdge edge:edgesGS) {
			AssemblyVertex vOutGS = edge.getConnectingVertex(vGS);
			if(vGS.getIndex()==vOutGS.getIndex()) continue;
			for(AssemblyEdge testEdge:testEdgesV1) {
				AssemblyVertex testOut = testEdge.getConnectingVertex(testV1);
				if(vOutGS.getIndex()==testOut.getIndex() && vOutGS.isStart()==testOut.isStart()) {
					count ++;
					break;
				}
			}
		}
		return count;
	}
}
