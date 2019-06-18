package ngsep.assembly;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
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
		//List of final contigs
		List<CharSequence> consensusList = new ArrayList<CharSequence>();
		for(int i = 0; i < graph.getPaths().size(); i++)
		{
			List<String> globalAlignment = new ArrayList<String>();
			List<AssemblyEdge> path = graph.getPaths().get(i);
			String consensus = "";
			int previousOverlap = 0;
			for(int j = 0; j < path.size(); j++)
			{
				//Needed to find which is the origin vertex
				AssemblyEdge previousEdge = null;
				if(j > 0)
					previousEdge = path.get(j - 1);
				AssemblyEdge edge = path.get(j);
				AssemblyVertex a = edge.getVertex1();
				AssemblyVertex b = edge.getVertex2();
				//If the first edge is being checked, compare to the second edge to find the origin vertex
				if(previousEdge == null && path.size() > 1)
				{
					previousEdge = path.get(j + 1);
					//The common vertex is the second vertex of the path
					if(previousEdge.getVertex1().getIndex() == edge.getVertex1().getIndex() || previousEdge.getVertex2().getIndex() == edge.getVertex1().getIndex())
					{
						a = edge.getVertex2();
						b = edge.getVertex1();
					}
				}
				else if (previousEdge != null)
				{
					//The common vertex is the first vertex to be compared
					if(previousEdge.getVertex1().getIndex() == edge.getVertex2().getIndex() || previousEdge.getVertex2().getIndex() == edge.getVertex2().getIndex())
					{
						a = edge.getVertex2();
						b = edge.getVertex1();
					}
				}
				//Ignore if both reads are equal
				if(!a.getRead().toString().equals(b.getRead().toString()))
				{
					//If the first string is start, then the reverse complement is added to the consensus
					String s1 = a.isStart() ? reverseComplement(a.getRead().toString()) : a.getRead().toString();
					String oriS1 = s1;
					//If the second string isn't start, then the reverse complement is added to the consensus
 					String s2 = !b.isStart() ? reverseComplement(b.getRead().toString()) : b.getRead().toString();
 					String oriS2 = s2;
 					//If it's starting the consensus, add the first characters that won't be compared to the consensus
 					if(j == 1 && s1.length() - edge.getOverlap() > tolerance)
 						consensus = s1.substring(0, s1.length() - edge.getOverlap() - tolerance);
 					//If the overlap + tolerance is greater than the length of the string, take the whole string
 					if(s1.length() - edge.getOverlap() > tolerance)
						s1 = s1.substring(s1.length() - edge.getOverlap() - tolerance);
					if(s2.length() - edge.getOverlap() > tolerance)
						s2 = s2.substring(0, edge.getOverlap() + tolerance);
					String[] alignments = getAlignment(s1, s2);
					
					if(j > 1 && previousOverlap + edge.getOverlap() + (2 * tolerance) < oriS1.length())
						consensus = consensus.concat(oriS1.substring(previousOverlap + tolerance, oriS1.length() - edge.getOverlap() - tolerance));
					System.out.println("J " + j);
					System.out.println("PREV " + previousOverlap);
					System.out.println("LEN " + oriS1.length());
					System.out.println("OL " + edge.getOverlap());
					System.out.println("ALI " + alignments[0].length());
					if(j == 483)	
					{
						int integ = 0;
					}
					String joined = joinedString(graph.getEmbedded(j), graph.getEmbedded(j+1), alignments);
					int startIndex = previousOverlap + tolerance - oriS1.length() + edge.getOverlap() + tolerance;
					System.out.println("IDX " + startIndex);
					System.out.println("JND " + joined.length());
					if(j > 1 && startIndex < joined.length())
						consensus = consensus.concat(joined.substring(previousOverlap + tolerance - oriS1.length() + edge.getOverlap() + tolerance));
					else if (j == 1)
						consensus = consensus.concat(joined);
					
					previousOverlap = edge.getOverlap();
					//If it's the last edge, add the remaining uncompared string from the last vertex
					if(j == path.size() - 2)
						consensus = consensus.concat(s2 = s2.substring(edge.getOverlap() + tolerance));
				}
			}
			
			consensusList.add(consensus);
		}
		
		return consensusList;
	}
	
	private String joinedString(List<AssemblyEmbedded> embedded1, List<AssemblyEmbedded> embedded2, String[] alignment)
	{
		StringBuilder finalString = new StringBuilder();
		int i = 0;
		if(startConsensus)
		{
			for(i = alignment[0].length() - 1; i >= 0; i--)
			{
				char b = alignment[0].charAt(i);
				if(b != '-')
					break;
			}
			startConsensus = false;
			i++;
		}
		
		i = 0;
		for (; i < alignment[0].length(); i++) 
		{
			char a = alignment[0].charAt(i);
			char b = alignment[1].charAt(i);
			if(a == '-' && b != '-')
			{
				finalString.append(b);
			}
			else if(a != '-' && b == '-')
			{
				finalString.append(a);
			}
			else if(a == b)
			{
				finalString.append(a);
			}
			else if(a != b)
			{
				if(embedded1 != null && embedded2 != null)
					finalString.append(embedded1.size() > embedded2.size() ? a : b);
				else if (embedded1 != null)
					finalString.append(a);
				else 
					finalString.append(b);
			}
		}
		return finalString.toString();
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
			Set<Integer> foundList = new TreeSet<Integer>();
			foundList.addAll(fm.exactSearch(sub));
			for(Integer x : foundList)
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
					else if(x < lastSearched)
						continue;
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
