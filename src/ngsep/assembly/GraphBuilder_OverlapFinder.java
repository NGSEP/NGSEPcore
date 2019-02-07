package ngsep.assembly;

import java.util.List;

import ngsep.sequences.FMIndex;

public interface GraphBuilder_OverlapFinder {

	public void calculate(List<CharSequence> seq, FMIndex index);

	public AssemblyGraph getGrap();

}
