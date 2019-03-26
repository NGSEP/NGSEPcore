package ngsep.assembly;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import ngsep.sequences.FMIndexSingleSequence;

public class ConsensusBuilderBidirectionalFMIndex implements ConsensusBuilder {
	int match;
	int openGap;
	int extGap;
	int mismatch;
	int windowSize;
	int tolerance;
	AlignmentAffineGap aligner;
	boolean startConsensus = true;
	
	public ConsensusBuilderBidirectionalFMIndex(int match, int openGap, int extGap, int mismatch, int windowSize, int tolerance) 
	{
		this.match = match;
		this.openGap = openGap;
		this.extGap = extGap;
		this.mismatch = mismatch;
		this.windowSize = windowSize;
		this.tolerance = tolerance;
		aligner = new AlignmentAffineGap(match, openGap, extGap, mismatch);
	}
	
	@Override
	public List<CharSequence> makeConsensus(AssemblyGraph graph) 
	{
		List<CharSequence> consensusList = new ArrayList<CharSequence>();
		for(int i = 0; i < graph.getPaths().size(); i++)
		{
			List<AssemblyEdge> path = graph.getPaths().get(i);
			String consensus = "";
			startConsensus = true;
			for(int j = 0; j < path.size(); j++)
			{
				AssemblyEdge previousEdge = null;
				if(j > 0)
					previousEdge = path.get(j - 1);
				AssemblyEdge edge = path.get(j);
				AssemblyVertex a = edge.getVertex1();
				AssemblyVertex b = edge.getVertex2();
				if(previousEdge == null && path.size() > j + 1)
				{
					previousEdge = path.get(j + 1);
					if(previousEdge.getVertex1().getIndex() == edge.getVertex1().getIndex() || previousEdge.getVertex2().getIndex() == edge.getVertex1().getIndex())
					{
						a = edge.getVertex2();
						b = edge.getVertex1();
					}
				}
				else if (previousEdge != null)
				{
					if(previousEdge.getVertex1().getIndex() == edge.getVertex2().getIndex() || previousEdge.getVertex2().getIndex() == edge.getVertex2().getIndex())
					{
						a = edge.getVertex2();
						b = edge.getVertex1();
					}
				}
				String s1 = a.isStart() ? a.getRead().toString() : complementaryStrand(a.getRead().toString());
				s1 = s1.substring(s1.length() - edge.getOverlap() - tolerance);
				String s2 = b.isStart() ? b.getRead().toString() : complementaryStrand(b.getRead().toString());
				s2 = s2.substring(0, edge.getOverlap() + tolerance);
				System.out.println(s2);
				String[] alignments = getAlignment(s1, s2);
			}
			consensusList.add(consensus);
		}	
		return consensusList;
	}
	
	private String joinedString()
	{
		return null;
	}
	
	private String[] getAlignment(String s1, String s2)
	{
		FMIndexSingleSequence fm = new FMIndexSingleSequence(s1);
		String seq1 = "";
		String seq2 = "";
		int lastSearched = 0;
		int start1 = -1;
		int start2 = -1;
		for(int i = 0; i < Math.round(Math.floor(s2.length() / windowSize)); i++)
		{
			int windowStart = i * windowSize;
			int windowEnd = windowStart + windowSize;
			boolean found = false;
			String sub = s2.substring(windowStart, windowEnd);
			for(Integer x : fm.exactSearch(sub))
			{
				if(x >= lastSearched - tolerance && x <= lastSearched + tolerance)
				{
					if(x > lastSearched && lastSearched == 0)
					{
						seq1 += s1.substring(0, x);
						seq2 += IntStream.range(0, x).mapToObj(j -> "-").collect(Collectors.joining(""));
						seq1 += s1.substring(x, x + windowSize);
						seq2 += sub;
					}
					else if(start1 != -1 && start2 != -1)
					{
						String sub1 = s1.substring(start1, x);
						String sub2 = s2.substring(start2, windowStart);
						String[] alignments = aligner.getAlignment(sub1, sub2);
						seq1 += alignments[0];
						seq2 += alignments[1];
						seq1 += s1.substring(x, x + windowSize);
						seq2 += sub;
						start1 = -1;
						start2 = -1;
					}
					else if (x > lastSearched)
					{
						seq1 += s1.substring(lastSearched, x - 1);
						seq2 += IntStream.range(0, x - lastSearched).mapToObj(j -> "-").collect(Collectors.joining(""));
						seq1 += s1.substring(x, x + windowSize);
						seq2 += sub;
					}
					else if(x == lastSearched)
					{
						seq1 += s1.substring(x, x + windowSize);
						seq2 += sub;
					}
					lastSearched = x + windowSize;
					found = true;
					break;
				}
			}
			if(!found)
			{
				if(start1 == -1)
					start1 = lastSearched;
				if(start2 == -1)
					start2 = windowStart;
				lastSearched += windowSize;
			}
		}

		if(start1 != -1 && start2 != -1)
		{
			String sub1 = s1.substring(start1, s1.length() - 1);
			String sub2 = s2.substring(start2, s2.length() - 1);
			String[] alignments = aligner.getAlignment(sub1, sub2);
			seq1 += alignments[0];
			seq2 += alignments[1];
		}
		else
		{
			String sub1 = s1.substring(lastSearched, s1.length());
			String sub2 = s2.substring(s2.length() - (s2.length() % windowSize), s2.length());
			if(sub1.equals(sub2))
			{
				seq1 += sub1;
				seq2 += sub2;
			}
			else
			{
				String[] alignments = aligner.getAlignment(sub1, sub2);
				seq1 += alignments[0];
				seq2 += alignments[1];
			}
		}
		String[] alignments = new String[2];
		alignments[0] = seq1;
		alignments[1] = seq2;
		return alignments;
	}
	
	private String complementaryStrand(String s)
	{
		StringBuilder complementaryStrand = new StringBuilder();
		for(int i = 0; i < s.length(); i++)
		{
			complementaryStrand.append(complementaryBase(s.charAt(i)));
		}
		return complementaryStrand.toString();
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
