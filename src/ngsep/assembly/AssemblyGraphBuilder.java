package ngsep.assembly;

import java.util.List;
import java.util.Map;

public interface AssemblyGraphBuilder {
	/**
	 * Identifies overlaps among the given sequences
	 * @param sequences to process
	 */
	public void findOverlaps (List<CharSequence> sequences);
	
	/**
	 * Returns the assembly graph obtained from the identified overlaps
	 * @return AssemblyGraph 
	 */
	public AssemblyGraph getAssemblyGraph();
	
	/**
	 * Returns the map of embedded sequences
	 * @return Map<Integer, Integer> Map with index of embedded sequence as key and index of host sequence as value
	 */
	public Map<Integer, Integer> getEmbeddedSequences (); 
}
