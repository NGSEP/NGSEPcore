package ngsep.test.consensus;

import ngsep.assembly.*;
import java.util.ArrayList;
import java.util.List;

public class ConsensusTest 
{
	public static void main(String[] args) 
	{
		//ACGTACGTTGCGACCCAGACGTGCGACGTCAGCTCGCAGTCGACGCTCAGCGCATGCATCGCATCGCAGCTACGACTGCATGCG
		List<CharSequence> sequences = new ArrayList<CharSequence>();
		sequences.add("ACGTACGTTGCGACCCAGACGTGCGA"); //0
		sequences.add("GTACGTTGCGACCCAGACGTGCGACGTCAGCT"); //1
		sequences.add("TTGCGACCCAGACGTGCGACGTCAGCTCGCAGT"); //2
		sequences.add("CGACCCAGACGTGCGACGTCAGCTCGCAGTCGACGCT"); //3
		//sequences.add("GTGCGACGTCAGCTCGCAGTCG"); //En anterior
		sequences.add("CGACGTCAGCTCGCAGTCGACGCTCA"); //5
		sequences.add("TCGCAGTCGACGCTCAGCGCATGCATCGCAT"); //6
		sequences.add("CGACGCTCAGCGCATGCATCGCATCGC"); //7
		sequences.add("GCGCATGCATCGCATCGCAGCTAC"); //8
		sequences.add("ATGCATCGCATCGCAGCTACGACTGCATGCG"); //9
		AssemblyGraph graph = new AssemblyGraph(sequences);
		List<Integer> path = new ArrayList<Integer>();
		path.add(0);
		path.add(2);
		path.add(4);
		path.add(6);
		path.add(8);
		path.add(10);
		path.add(12);
		path.add(14);
		path.add(16);
		//graph.addPath(path);

		ConsensusBuilder consensusBuilder = new ConsensusBuilderBidirectionalAffineGap();
		//ConsensusBuilder consensusBuilder = new ConsensusBuilderBidirectionalConstantGap();
		List<CharSequence> consensus = consensusBuilder.makeConsensus(graph);
	}

}
