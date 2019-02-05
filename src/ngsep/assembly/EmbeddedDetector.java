package ngsep.assembly;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ngsep.sequences.DNAMaskedSequence;
import ngsep.sequences.FMIndex;

@FunctionalInterface
public interface EmbeddedDetector {
	public Map<Integer, Edge> getEmbedded(FMIndex index, List<DNAMaskedSequence> sequences);

	public static EmbeddedDetector NONE = (FMIndex index, List<DNAMaskedSequence> sequences) -> {
		return new HashMap<Integer, Edge>();
	};
}
