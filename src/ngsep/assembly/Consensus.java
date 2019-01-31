package ngsep.assembly;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import ngsep.sequences.DNAMaskedSequence;

public interface Consensus {
	public List<CharSequence> makeConsensus(List<List<Integer>> paths, List<DNAMaskedSequence> sequences,
			Map<Integer, EmbeddedSequence> embedded, Map<Integer, List<Edge>> edges);

	public static Consensus NONE = (List<List<Integer>> paths, List<DNAMaskedSequence> sequences,
			Map<Integer, EmbeddedSequence> embedded, Map<Integer, List<Edge>> edges) -> {
		return new LinkedList<CharSequence>();
	};
}
