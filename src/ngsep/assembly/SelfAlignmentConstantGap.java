package ngsep.assembly;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import ngsep.sequences.FMIndexSingleSequence;

/**
 * Alignment of a main read with its embedded reads
 */
public class SelfAlignmentConstantGap
{
	private int windowSize;
	AlignmentConstantGap aligner;
	
	public SelfAlignmentConstantGap(int match, int gap, int mismatch, int windowSize, int tolerance)
	{
		this.windowSize = windowSize;
		this.aligner = new AlignmentConstantGap(match, gap, mismatch);
	}
		
	/**
	 * Sets the alignment of the reads with the main read
	 * @param mainRead Main read that contains every embedded read
	 * @param embeddedReads List of embedded reads
	 * @param match Score for a match
	 * @param gap Score for a gap
	 * @param mismatch Score for a mismatch
	 * @return alignments Aligned strings
	 */
	public String[] selfAlign(String mainRead, List<AssemblyEmbedded> embeddedReads, boolean reverse)
	{
		StringBuilder[] alignmentsB = new StringBuilder[embeddedReads.size() + 1];
		for(int i = 0; i < embeddedReads.size(); i++)
		{
			AssemblyEmbedded embeddedRead = embeddedReads.get(i);
			String embedded = embeddedRead.getRead().toString();
			int start = embeddedRead.getStartPosition();
			//If the main and the embedded read have different senses, the embedded read must be reversed
			if(reverse != embeddedRead.isReverse())
			{
				embedded = reverseComplement(embedded);
			}
			//If the main read is reversed, the position 0 of the subalignment is the end of the main read,
			//so the start of the alignment is the end of the main read - the beginning of the embedding - the length 
			//of the embedded read
			if(reverse)
			{
				start = mainRead.length() - start - embedded.length();
			}
			FMIndexSingleSequence fm = new FMIndexSingleSequence(mainRead);
			StringBuilder[] alignedReads = getAlignment(fm, mainRead, embedded, start);
//			System.out.println("MAIN");
//			System.out.println(alignedReads[0]);
//			System.out.println(alignedReads[1]);
			if(i == 0)
			{
				alignmentsB[0] = new StringBuilder(alignedReads[0]);
			}
			else if (alignedReads.length > 0)
			{
				int offset = 0;
				for(int j = 0; j < alignmentsB[0].length(); j++)
				{
					char a = alignmentsB[0].charAt(j);
					char b = '-';
					if (j - offset >= alignedReads[0].length())
					{
						break;
					}
					else
					{
						b = alignedReads[0].charAt(j - offset);
					}
					
					if(a != b)
					{
						if(a == '-')
						{
							if((j > 0 && alignedReads[1].charAt(j - 1) == '_') || (alignedReads[1].charAt(j) == '_'))
								alignedReads[1].insert(j, '_');
							else
								alignedReads[1].insert(j, '-');
							offset++;
						}
						else if (b == '-')
						{
							for(int k = 0; k < alignmentsB.length; k++)
							{
								if(alignmentsB[k] != null)
								{
									if((j > 0 && alignmentsB[k].charAt(j - 1) == '_') || (alignmentsB[k].charAt(j) == '_'))
										alignmentsB[k].insert(j, '_');
									else
										alignmentsB[k].insert(j, '-');
								}
							}
						}
					}
				}
			}
			alignmentsB[i + 1] = new StringBuilder();
			alignmentsB[i + 1].append(alignedReads[1]);
		}
		String[] alignments = new String[alignmentsB.length];
		if(embeddedReads.size() == 0)
		{
			alignments = new String[]{mainRead};
		}
		else
		{
			for(int i = 0; i < alignmentsB.length; i++)
			{
				alignments[i] = alignmentsB[i].toString();
			}
		}
		System.out.println("ALIGNMENT");
		for(int i = 0; i < alignments.length; i++)
		{
			System.out.println("ALN " + i + " " + alignments[i]);
		}
		return alignments;
 	}
	
	private StringBuilder[] getAlignment(FMIndexSingleSequence fm, String s1, String s2, int start)
	{
		StringBuilder seq1 = new StringBuilder();
		StringBuilder seq2 = new StringBuilder();
		int lastSearched = start;
		int start1 = -1;
		int start2 = -1;
		int end = start + s2.length();
		seq1.append(s1.substring(0, start));
		for(int j = 0; j < start; j++)
		{
			seq2.insert(j, '_');
		}
		int windowStart = 0;
		for(int i = 0; i < Math.round(Math.floor(s2.length() / windowSize)); i++)
		{
			windowStart = i * windowSize;
			int windowEnd = windowStart + windowSize;
			boolean found = false;
			String sub = s2.substring(windowStart, windowEnd);
			Set<Integer> foundList = new TreeSet<Integer>();
			foundList.addAll(fm.exactSearch(sub));
			for(Integer x : foundList)
			{
				if(x >= lastSearched - windowSize && x <= lastSearched + windowSize)
				{
					if(x > lastSearched && lastSearched == 0)
					{
						seq1.append(s1.substring(0, x));
						seq2.append(IntStream.range(0, x).mapToObj(j -> "-").collect(Collectors.joining("")));
						seq1.append(s1.substring(x, x + windowSize));
						seq2.append(sub);
					}
					else if(start1 != -1 && start2 != -1)
					{
						String sub1 = s1.substring(start1, x);
						String sub2 = s2.substring(start2, windowStart);
						String[] alignments = aligner.getAlignment(sub1, sub2);
						seq1.append(alignments[0]);
						seq2.append(alignments[1]);
						seq1.append(s1.substring(x, x + windowSize));
						seq2.append(sub);
						start1 = -1;
						start2 = -1;
					}
					else if (x > lastSearched)
					{
						seq1.append(s1.substring(lastSearched, x));
						seq2.append(IntStream.range(0, x - lastSearched).mapToObj(j -> "-").collect(Collectors.joining("")));
						seq1.append(s1.substring(x, x + windowSize));
						seq2.append(sub);
					}
					else if(x == lastSearched)
					{
						seq1.append(s1.substring(x, x + windowSize));
						seq2.append(sub);
					}
					else if(x < lastSearched)
						continue;
					lastSearched = x + windowSize;
					found = true;
					break;
				}
			}
			if(!found && lastSearched <= end + windowSize)
			{
				if(start1 == -1)
					start1 = lastSearched;
				if(start2 == -1)
					start2 = windowStart;
				lastSearched += windowSize;
			}
		}
//		System.out.println("SEQ1 " + seq1);
//		System.out.println("SEQ2 " + seq2);
		if(start1 != -1 && start2 != -1)
		{
			String sub1 = "";
			if(s1.length() > lastSearched + windowSize)
			{
				sub1 = s1.substring(start1, lastSearched + windowSize);
			}
			else
			{
				sub1 = s1.substring(start1, s1.length());
			}
			String sub2 = "";
			if(s2.length() > windowStart + windowSize)
			{
				sub2 = s2.substring(start2, windowStart + windowSize);
			}
			else 
			{
				sub2 = s2.substring(start2, s2.length());
			}
			String[] alignments = aligner.getAlignment(sub1, sub2);
			seq1.append(alignments[0]);
			seq2.append(alignments[1]);
			lastSearched += windowSize;
		}
		else
		{
			String sub1 = s1.substring(lastSearched, s1.length());
			if(lastSearched + windowSize <= s1.length())
			{
				sub1 = s1.substring(lastSearched, lastSearched + windowSize);
			}
			String sub2 = s2.substring(s2.length() - (s2.length() % windowSize), s2.length());
			if(sub1.equals(sub2))
			{
				seq1.append(sub1);
				seq2.append(sub2);
			}
			else
			{
				String[] alignments = aligner.getAlignment(sub1, sub2);
				seq1.append(alignments[0]);
				seq2.append(alignments[1]);
			}
			lastSearched += windowSize;
		}
		if(lastSearched < s1.length())
		{
			seq1.append(s1.substring(lastSearched, s1.length()));
		}
		for(int j = lastSearched; j < s1.length(); j++)
		{
			seq2.append('_');
		}
		StringBuilder[] alignments = new StringBuilder[2];
		alignments[0] = seq1;
		alignments[1] = seq2;
		return alignments;
	}
	
	private Set<Integer> getMatches(String s1, String s2)
	{
		int lastIndex = 0;
		Set<Integer> result = new TreeSet<Integer>();
		while(lastIndex != -1) 
		{
		    lastIndex = s1.indexOf(s2,lastIndex);
		    if(lastIndex != -1)
		    {
		        result.add(lastIndex);
		        lastIndex += 1;
		    }
		}
		return result;
	}
	
	private static String reverseComplement(String s)
	{
		StringBuilder complementaryStrand = new StringBuilder();
		for(int i = 0; i < s.length(); i++)
		{
			complementaryStrand.append(complementaryBase(s.charAt(i)));
		}
		return complementaryStrand.reverse().toString();
	}
	
	private static char complementaryBase(char b)
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
