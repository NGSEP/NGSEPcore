package ngsep.assembly;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;

@FunctionalInterface
public interface OverlappingDetector {
	Map<Integer, List<Edge>> getEdges(FMIndex index, List<DNAMaskedSequence> sequences, Set<Integer> ignored);

	public static OverlappingDetector NONE = (FMIndex index, List<DNAMaskedSequence> sequences,
			Set<Integer> ignored) -> {
		return new HashMap<Integer, List<Edge>>();
	};
}
