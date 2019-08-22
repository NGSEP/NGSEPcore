package ngsep.assembly;

import static ngsep.assembly.TimeUtilities.timeGroup;
import static ngsep.assembly.TimeUtilities.timeIt;

import java.io.IOException;
import java.util.List;

public interface GraphBuilder {
	/**
	 * Builds an assembly graph from the given sequences
	 * 
	 * @param sequences to process
	 */
	public SimplifiedAssemblyGraph buildAssemblyGraph(List<CharSequence> sequences, AssemblyConfiguration config);

	public static void main(String[] args) throws Exception {
		List<CharSequence> sequences = timeIt("Load the sequences", () -> Assembler.load(args[0]));

		SimplifiedAssemblyGraph assemblyGraph = timeGroup("Build overlap Graph", () -> {
			AssemblyConfiguration config = (args.length > 3)
					? new AssemblyConfiguration(Double.valueOf(args[2]), Double.valueOf(args[3]))
					: new AssemblyConfiguration();
			return (new GraphBuilderFMIndex()).buildAssemblyGraph(sequences, config);

		});

		timeIt("Save the Graph", () -> {
			try {
				assemblyGraph.save(args[1]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		});

		System.out.println("---------Graph Properties--------------");
		assemblyGraph.printInfo();
		System.out.println("---------------------------------------");
		
	}
}
