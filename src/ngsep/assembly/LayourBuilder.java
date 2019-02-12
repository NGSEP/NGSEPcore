package ngsep.assembly;

@FunctionalInterface
public interface LayourBuilder {
	public void findPaths(AssemblyGraph graph);

	public static LayourBuilder NONE = (AssemblyGraph graph) -> {
		
	};
}
