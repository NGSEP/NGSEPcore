package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;

public class TestConsensusBuilder {

	public TestConsensusBuilder() {
		// TODO Auto-generated constructor stub
	}


	public static void main(String[] args) 
	{
		//ACGTACGTTGCGACCCAGACGTGCGACGTCAGCTCGCAGTCGACGCTCAGCGCATGCATCGCATCGCAGCTACGACTGCATGCG
		List<CharSequence> sequences = new ArrayList<CharSequence>();
		//sequences.add("ACGTACGTTGCGACCCAGACGTGCGA"); //0
		//sequences.add("GTACGTTGCGACCCAGACGTGCGACGTCAGCT"); //1
		sequences.add("TTGCGACCCAGACGTGCGACGTCAGCTCGCAGT"); //2
		   sequences.add("CGACCCAGTCGTGCGACGTCAGCTCGCAAGTCGACGCT"); //3
		////sequences.add("GTGCGACGTCAGCTCGCAGTCG"); //En anterior
		//sequences.add("CGACGTCAGCTCGCAGTCGACGCTCA"); //5
		//sequences.add("TCGCAGTCGACGCTCAGCGCATGCATCGCAT"); //6
		//sequences.add("CGACGCTCAGCGCATGCATCGCATCGC"); //7
		//sequences.add("GCGCATGCATCGCATCGCAGCTAC"); //8
		//sequences.add("ATGCATCGCATCGCAGCTACGACTGCATGCG"); //9
		AssemblyGraph graph = new AssemblyGraph(sequences);
		List<AssemblyEdge> path = new ArrayList<AssemblyEdge>();
		path.add(new AssemblyEdge(graph.getVertex(0, true), graph.getVertex(1, true), 30));
		graph.addPath(path);

		ConsensusBuilder consensusBuilder = new ConsensusBuilderBidirectionalFMIndex(2, 20, 1, 8, 7, 2);
		List<CharSequence> consensus = consensusBuilder.makeConsensus(graph);
	}

}
