package ngsep.assembly;

import java.util.List;

public interface GraphBuilder {
	/**
	 * Builds an assembly graph from the given sequences
	 * 
	 * @param sequences to process
	 */
	public AssemblyGraph buildAssemblyGraph(List<CharSequence> sequences);
}
