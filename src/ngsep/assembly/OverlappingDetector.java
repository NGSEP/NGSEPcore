package ngsep.assembly;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;

@FunctionalInterface
public interface OverlappingDetector {
	Map<Integer, List<Edge>> getEdges(FMIndex index, List<DNAMaskedSequence> sequences, Map<Integer, Edge> emmbed);

	public static OverlappingDetector NONE = (FMIndex index, List<DNAMaskedSequence> sequences,
			Map<Integer, Edge> ignored) -> {
		return new HashMap<Integer, List<Edge>>();
	};
}
