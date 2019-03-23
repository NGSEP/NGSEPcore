package ngsep.assembly;

import java.util.List;

import ngsep.sequences.FMIndex;

public interface GraphBuilderOverlapFinder {

	public void calculate(List<CharSequence> seq, FMIndex index);

	public AssemblyGraph getGrap();

}
