package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;


public interface ConsensusBuilder {
	public List<CharSequence> makeConsensus(List<List<Integer>> paths, List<CharSequence> sequences,
			AssemblyGraph graph);
	
	public static ConsensusBuilder NONE = (List<List<Integer>> paths, List<CharSequence> sequences,
			AssemblyGraph graph) -> {
		return new ArrayList<>();
	};
}
