package ngsep.assembly;

import java.util.List;
import java.util.logging.Logger;

import ngsep.sequences.ByteArrayKmersMapImpl;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;
import ngsep.sequences.KmersExtractor;
import ngsep.sequences.KmersMap;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class GraphBuilderFMIndexKmers implements GraphBuilder {

	private Logger log = Logger.getLogger(GraphBuilderFMIndex.class.getName());
	private final static int TALLY_DISTANCE = 100;
	private final static int SUFFIX_FRACTION = 20;
	
	private static final int TIMEOUT_SECONDS = 30;

	private int kmerLength=15;
	private int numThreads;
	
private static int idxDebug = -1;
	
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
	
	@Override
	public AssemblyGraph buildAssemblyGraph(List<QualifiedSequence> sequences) {
		AssemblyGraph graph = new AssemblyGraph(sequences);
		log.info("Created graph vertices. Edges: "+graph.getEdges().size());
		KmersExtractor extractor = new KmersExtractor();
		extractor.setKmerLength(kmerLength);
		extractor.processQualifiedSequences(sequences);
		KmersMap map = extractor.getKmersMap();
		// Create FM-Index
		FMIndex fmIndexForward = new FMIndex();
		// TODO: Set tally distance and suffix fraction
		fmIndexForward.loadQualifiedSequences(sequences);
		log.info("Created FM-Index forward");
		QualifiedSequenceList reverseComplements = new QualifiedSequenceList();
		for(QualifiedSequence seq:sequences) {
			reverseComplements.add(new QualifiedSequence(seq.getName(),DNAMaskedSequence.getReverseComplement(seq.getCharacters())));
		}
		FMIndex fmIndexReverse = new FMIndex();
		fmIndexReverse.loadQualifiedSequences(reverseComplements);
		log.info("Created FM-Index reverse");
		List<AssemblyEdge> edges = calculateEdges(fmIndexForward, fmIndexReverse,map);
		for(AssemblyEdge edge:edges) graph.addEdge(edge);
		return graph;
	}

	private KmersMap generateKmersMap(List<QualifiedSequence> sequences) {
		KmersExtractor extractor = new KmersExtractor();
		extractor.setKmerLength(kmerLength);
		for (QualifiedSequence qseq: sequences) {
			extractor.countSequenceKmers(qseq.getCharacters());
		}
		return extractor.getKmersMap();
	}
	
	private List<AssemblyEdge> calculateEdges(FMIndex fmIndexForward, FMIndex fmIndexReverse, KmersMap map) {
		// TODO Auto-generated method stub
		return null;
	}

}
