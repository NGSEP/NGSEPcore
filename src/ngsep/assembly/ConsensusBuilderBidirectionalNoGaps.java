package ngsep.assembly;

import java.util.ArrayList;
import java.util.List;

public class ConsensusBuilderBidirectionalNoGaps implements ConsensusBuilder {
	int match;
	int openGap;
	int extGap;
	int mismatch;
	int windowSize;
	int tolerance;
	AlignmentAffineGap aligner;
	boolean startConsensus = true;
	
	public ConsensusBuilderBidirectionalNoGaps(int match, int openGap, int extGap, int mismatch, int windowSize, double Rate_of_changes, double Rate_of_cuts, double Rate_of_cover) 
	{
		double rate_of_error = Rate_of_changes + Rate_of_cuts - Rate_of_changes * Rate_of_cuts;
		this.match = match;
		this.openGap = openGap;
		this.extGap = extGap;
		this.mismatch = mismatch;
		this.windowSize = windowSize;
		this.tolerance = (int) (Rate_of_cuts * (11.51292546 / rate_of_error));
		aligner = new AlignmentAffineGap(match, openGap, extGap, mismatch);
	}
	
	@Override
	public List<CharSequence> makeConsensus(AssemblyGraph graph) 
	{
		//List of final contigs
		List<CharSequence> consensusList = new ArrayList<CharSequence>();
		List<List<AssemblyEdge>> paths = graph.getPaths(); 
		for(int i = 0; i < paths.size(); i++)
		{
			List<AssemblyEdge> path = paths.get(i);
			CharSequence consensusPath = makeConsensus (graph, path);
			consensusList.add(consensusPath);
		}
		
		return consensusList;
	}
	
	private CharSequence makeConsensus(AssemblyGraph graph, List<AssemblyEdge> path) 
	{
		StringBuilder consensus = new StringBuilder();
		List<byte[]> consensusCounts = new ArrayList<byte[]>();
		AssemblyVertex lastVertex = null;
		int currentPos = 0;
		for(int j = 0; j < path.size(); j++)
		{
			if(j % 10 == 0)
				System.out.println("J = " + (String.format("%.2f", (j * 100.0 / path.size())) + "%"));
			//Needed to find which is the origin vertex
			AssemblyEdge edge = path.get(j);
			AssemblyVertex a = edge.getVertex1();
			AssemblyVertex b = edge.getVertex2();
			//If the first edge is being checked, compare to the second edge to find the origin vertex
			if(j == 0)
			{
				AssemblyEdge nextEdge = path.get(j + 1);
				//The common vertex is the second vertex of the path
				if(nextEdge.getVertex1().getIndex() == edge.getVertex1().getIndex() || nextEdge.getVertex2().getIndex() == edge.getVertex1().getIndex())
				{
					a = edge.getVertex2();
					b = edge.getVertex1();
				}
			}
			else if(lastVertex == b)
			{
				//The common vertex is the first vertex to be compared
				a = edge.getVertex2();
				b = edge.getVertex1();
			}
			if(j > 0 && lastVertex != a) 
			{
				throw new RuntimeException("Inconsistency found in path");
			}
			
			//If it's the first edge in the path, add the first read base count to the consensus
			if(j == 0) 
			{
				int seqId = (int) Math.floor(a.getIndex() / 2.0);
				List<AssemblyEmbedded> embeddedReads = graph.getEmbedded(seqId);
				if(embeddedReads == null)
					embeddedReads = new ArrayList<AssemblyEmbedded>();
				String seqA = a.isStart() ? a.getRead().toString() : reverseComplement(a.getRead().toString());
				boolean reverse = !a.isStart();
				consensusCounts.addAll(consensusReadEmbeddeds(seqA, embeddedReads, reverse));
				currentPos = seqA.length();
			} 
			else if(a.getRead() != b.getRead())
			{
				int seqId = (int) Math.floor(b.getIndex() / 2.0);
				List<AssemblyEmbedded> embeddedReads = graph.getEmbedded(seqId);
				if(embeddedReads == null)
					embeddedReads = new ArrayList<AssemblyEmbedded>();
				//If the second string isn't start, then the reverse complement is added to the consensus
				String nextSequence = b.isStart() ? b.getRead().toString() : reverseComplement(b.getRead().toString());
				boolean reverse = !b.isStart();
				
				List<byte[]> consensusEmbedded = consensusReadEmbeddeds(nextSequence, embeddedReads, reverse);
				
				currentPos = currentPos - edge.getOverlap();
				for(int i = 0; i < consensusEmbedded.size(); i++)
				{
					if(consensusCounts.size() == currentPos + i)
						consensusCounts.add(new byte[4]);
					byte[] consensusEmbeddedBase = consensusEmbedded.get(i);
					for(int k = 0; k < consensusEmbeddedBase.length; k++)
					{
						//Add the base count of the embeddeds count to the global consensus for each base position
						byte[] consensusBase = consensusCounts.get(currentPos + i);
						consensusBase[k] += consensusEmbeddedBase[k];
					}
				}
				currentPos = currentPos + nextSequence.length();
				
				if(nextSequence.length() - edge.getOverlap() <= tolerance) 
				{
					System.err.println("Non embedded edge has overlap: " + edge.getOverlap() + " and length: " + nextSequence.length());
				} 
			}
			lastVertex = b;
		}
		
		//Find the most common occurrence in each position and add the most common character to the consensus
		for(byte[] counts : consensusCounts)
		{
			int maxIndex = 0;
			for(int m = 1; m < counts.length; m++)
			{
				if(counts[m] > counts[maxIndex])
				{
					maxIndex = m;
				}
			}
			char maxChar = '-';
			switch (maxIndex)
			{
				case 0:
					maxChar = 'A';
					break;
				case 1:
					maxChar = 'T';
					break;
				case 2:
					maxChar = 'C';
					break;
				case 3:
					maxChar = 'G';
					break;
			}
			consensus.append(maxChar);
		}
		return consensus;
	}
	
	private List<byte[]> consensusReadEmbeddeds(String read, List<AssemblyEmbedded> embeddedReads, boolean reverse)
	{
		List<byte[]> consensusCounts = new ArrayList<byte[]>();
		//Align the read a with its embedded reads
		String[] alignedOriginRead = SelfAlignment.selfAlign(read, embeddedReads, reverse, match, openGap, extGap, mismatch);
		for(int k = 0; k < alignedOriginRead.length; k++)
		{
			//Get each aligned read
			String alignedRead = alignedOriginRead[k].toUpperCase();
			//Get each character of the aligned read
			for(int l = 0; l < read.length(); l++)
			{
				char alignedChar = alignedRead.charAt(l);
				//Add a new element to the consensus count if its size is less 
				//than the index of the current character relative to the substring that will be added
				if(consensusCounts.size() <= l)
				{
					//0 for A, 1 for T, 2 for C, 3 for G
					byte[] counts = new byte[4];
					consensusCounts.add(counts);
				}
				byte[] counts = consensusCounts.get(l);
				//Depending on the character, sum to the count in the current position
				switch (alignedChar)
				{
					case 'A':
						counts[0]++;
						break;
					case 'T':
						counts[1]++;
						break;
					case 'C':
						counts[2]++;
						break;
					case 'G':
						counts[3]++;
						break;
				}
			}
		}
		return consensusCounts;
	}
	
	private String reverseComplement(String s)
	{
		StringBuilder complementaryStrand = new StringBuilder();
		for(int i = 0; i < s.length(); i++)
		{
			complementaryStrand.append(complementaryBase(s.charAt(i)));
		}
		return complementaryStrand.reverse().toString();
	}
	
	private char complementaryBase(char b)
	{
		char complementaryBase;
		if(b == 'A')
			complementaryBase = 'T';
		else if(b == 'T')
			complementaryBase = 'A';
		else if(b == 'C')
			complementaryBase = 'G';
		else if(b == 'G')
			complementaryBase = 'C';
		else if (b == '-')
			complementaryBase = '-';
		else
			complementaryBase = 'N';
		return complementaryBase;
	}
}
