package ngsep.assembly;
import java.util.List;

/**
 * Alignment of a main read with its embedded reads
 */
public class SelfAlignment 
{
	String mainRead;
	String[] alignments;
	AlignmentAffineGap aligner;
	
	/**
	 * Self alignment constructor
	 * @param mainRead Main read that contains every embedded read
	 * @param embeddedReads List of embedded reads
	 * @param match Score for a match
	 * @param openGap Score for a gap opening
	 * @param extGap Score for a gap extension
	 * @param mismatch Score for a mismatch
	 */
	public SelfAlignment(String mainRead, List<AssemblyEmbedded> embeddedReads, int match, int openGap, int extGap, int mismatch) 
	{
		this.mainRead = mainRead;
		alignments = new String[embeddedReads.size() + 1];
		aligner = new AlignmentAffineGap(match, openGap, extGap, mismatch);
		selfAlign(embeddedReads);
	}
	
	/**
	 * Sets the alignment of the reads with the main read
	 * @param embeddedReads List of embedded reads
	 */
	private void selfAlign(List<AssemblyEmbedded> embeddedReads)
	{
		StringBuilder[] alignmentsB = new StringBuilder[embeddedReads.size() + 1];
		for(int i = 0; i < embeddedReads.size(); i++)
		{
			String s = embeddedReads.get(i).getRead().toString();
			String[] alignedReads = aligner.getAlignment(mainRead, s);
			if(i == 0)
			{
				alignmentsB[0] = new StringBuilder(alignedReads[0]);
				alignmentsB[1] = new StringBuilder(alignedReads[1]);
			}
			else
			{
				alignmentsB[i + 1] = new StringBuilder();
				StringBuilder s1 = new StringBuilder(alignedReads[0]);
				StringBuilder s2 = new StringBuilder(alignedReads[1]);
				int j = 0;
				int k = 0;
				while(j < s1.length())
				{
					char a = alignmentsB[0].charAt(k);
					char b = s1.charAt(j);
					if(a == b)
					{
						alignmentsB[i + 1].append(s2.charAt(j));
						j++;
					}
					else
					{
						if(a == '-')
						{
							alignmentsB[i + 1].append('-');
						}
						else if (b == '-')
						{
							for(int l = 0; l < i + 1; l++)
							{
								alignmentsB[l].insert(k, '-');
							}
							alignmentsB[i + 1].append(s2.charAt(j));
							j++;
						}
					}
					k++;
				}
			}
		}
		for(int i = 0; i < alignmentsB.length; i++)
		{
			System.out.println(alignmentsB[i]);
		}
 	}
}
