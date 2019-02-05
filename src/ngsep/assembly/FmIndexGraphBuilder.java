package ngsep.assembly;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import ngsep.sequences.FMIndex;

public class FmIndexGraphBuilder implements AssemblyGraphBuilder {
	private final static int TALLY_DISTANCE = 100;
	private final static int SUFFIX_FRACTION = 25;

	private AssemblyGraph graph;
	private Map<Integer, Integer> embeddedMap;

	@Override
	public void findOverlaps(List<CharSequence> sequences) {
		Collections.sort(sequences, (CharSequence l1, CharSequence l2) -> l2.length()
						- l1.length());
		System.out.println("	building FMIndexes");
		long ini = System.currentTimeMillis();
		FMIndex index = new FMIndex();
		index.loadUnnamedSequences(sequences, TALLY_DISTANCE, SUFFIX_FRACTION);
		System.out.println("	build FMIndexes: "
				+ (System.currentTimeMillis() - ini) / (double) 1000 + " s");

		
		/*EmbeddedDetector embeddedDetector = EmbeddedDetector.NONE;
		embedded = embeddedDetector.getEmbedded(index, sequences);

		OverlappingDetector overlappingDetector = new OverlapingDetector1();
		edges = overlappingDetector.getEdges(index, sequences, embedded);
		*/
	}

	@Override
	public AssemblyGraph getAssemblyGraph() {
		// TODO Auto-generated method stub
		return graph;
	}

	@Override
	public Map<Integer, Integer> getEmbeddedSequences() {
		// TODO Auto-generated method stub
		return embeddedMap;
	}
}
