package ngsep.assembly;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;

public class FmIndexEdgesFinder implements EdgesFinder {
	private final static int TALLY_DISTANCE = 100;
	private final static int SUFFIX_FRACTION = 25;

	private Map<Integer, EmbeddedSequence> embedded;
	private Map<Integer, List<Edge>> edges;

	public FmIndexEdgesFinder(List<DNAMaskedSequence> sequences) {
		Collections.sort(sequences,
				(DNAMaskedSequence l1, DNAMaskedSequence l2) -> l2.length()
						- l1.length());
		System.out.println("	building FMIndexes");
		long ini = System.currentTimeMillis();
		FMIndex index = new FMIndex();
		index.loadUnnamedSequences(sequences, TALLY_DISTANCE, SUFFIX_FRACTION);
		System.out.println("	build FMIndexes: "
				+ (System.currentTimeMillis() - ini) / (double) 1000 + " s");

		EmbeddedDetector embeddedDetector = EmbeddedDetector.NONE;
		embedded = embeddedDetector.getEmbedded(index, sequences);

		OverlappingDetector overlappingDetector = OverlappingDetector.NONE;
		edges = overlappingDetector.getEdges(index, sequences,
				embedded.keySet());
	}

	@Override
	public Map<Integer, List<Edge>> getEdges() {
		return edges;
	}

	@Override
	public Map<Integer, EmbeddedSequence> getEmbedded() {
		return embedded;
	}
}
