package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;

@FunctionalInterface
public interface ConsensusBuilder {
	public List<CharSequence> makeConsensus(AssemblyGraph graph);

	public static ConsensusBuilder NONE = (AssemblyGraph graph) -> {
		return new ArrayList<>();
	};
}
