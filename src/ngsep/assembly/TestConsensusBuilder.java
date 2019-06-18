package ngsep.assembly;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

public class TestConsensusBuilder {

	public TestConsensusBuilder() {
		// TODO Auto-generated constructor stub
	}


	public static void main(String[] args) 
	{
		/*
		//ACGTACGTTGCGACCCAGACGTGCGACGTCAGCTCGCAGTCGACGCTCAGCGCATGCATCGCATCGCAGCTACGACTGCATGCG
		List<CharSequence> sequences = new ArrayList<CharSequence>();
		//sequences.add("ACGTACGTTGCGACCCAGACGTGCGA"); //0
		//sequences.add("GTACGTTGCGACCCAGACGTGCGACGTCAGCT"); //1
		sequences.add("TTGCGACCCAGACGTGCGACGTCAGCTCGCAGT"); //2
		sequences.add(   "CGACCCAGTCGTGCGACGTCAGCTCGCAAGTCGACGCT"); //3
		////sequences.add("GTGCGACGTCAGCTCGCAGTCG"); //En anterior
		sequences.add(                                  "CGACGTCAGCTCGCAGTCGACGCTCA"); //5
		             //TTGCGACCCAGTCGTGCGACGTCAGCTCGCAAGTC   GTC           GACGTCA

		//sequences.add("TCGCAGTCGACGCTCAGCGCATGCATCGCAT"); //6
		//sequences.add("CGACGCTCAGCGCATGCATCGCATCGC"); //7
		//sequences.add("GCGCATGCATCGCATCGCAGCTAC"); //8
		//sequences.add("ATGCATCGCATCGCAGCTACGACTGCATGCG"); //9
		AssemblyGraph graph = new AssemblyGraph(sequences);
		List<AssemblyEdge> path = new ArrayList<AssemblyEdge>();
		path.add(new AssemblyEdge(graph.getVertex(0, true), graph.getVertex(1, true), 31));
		path.add(new AssemblyEdge(graph.getVertex(1, true), graph.getVertex(2, true), 7));
		graph.addPath(path);
		
		List<AssemblyEmbedded> aeList = new ArrayList<AssemblyEmbedded>();
		AssemblyEmbedded ae = new AssemblyEmbedded("GTCGTGCCCGACGTCAGTTCGCAA", 0, false);
		aeList.add(ae);
		ae = new AssemblyEmbedded("ACGTCAGCTTCAAGCAAGTCGACGA", 0, false);
		aeList.add(ae);
		
		//SelfAlignment sa = new SelfAlignment("CGACCCAGTCGTGCGACGTCAGCTTCGCAAGTCGACGA", aeList, 2, 20, 1, 8);

		ConsensusBuilder consensusBuilder = new ConsensusBuilderBidirectionalFMIndex(2, 20, 1, 8, 7, 2);
		List<CharSequence> consensus = consensusBuilder.makeConsensus(graph);
		*/
		try
		{
			Assembler assembler = new Assembler(args[0], args[1]);
		}
		catch(Exception e)
		{
			try
			{
				StringWriter writer = new StringWriter();
				PrintWriter print = new PrintWriter(writer);
				e.printStackTrace(print);
				FileWriter out = new FileWriter("D:\\assembler.log");
				out.write(writer.toString());
				out.close();
			}
			catch(IOException i2)
			{
				i2.printStackTrace();
			}
		}
	}

}
