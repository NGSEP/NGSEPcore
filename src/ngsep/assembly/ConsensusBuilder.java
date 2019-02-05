package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;


public interface ConsensusBuilder {
	public List<CharSequence> makeConsensus(List<List<Integer>> paths, List<CharSequence> sequences,
			Map<Integer, Integer> embedded, AssemblyGraph graph);
	
	public static ConsensusBuilder NONE = (List<List<Integer>> paths, List<CharSequence> sequences,
			Map<Integer, Integer> embedded, AssemblyGraph graph) -> {
		return new ArrayList<>();
	};
}
