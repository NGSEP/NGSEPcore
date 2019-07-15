package ngsep.assembly;
import java.util.List;

/**
 * Alignment of a main read with its embedded reads
 */
public class SelfAlignment 
{
	/**
	 * Sets the alignment of the reads with the main read
	 * @param mainRead Main read that contains every embedded read
	 * @param embeddedReads List of embedded reads
	 * @param match Score for a match
	 * @param openGap Score for a gap opening
	 * @param extGap Score for a gap extension
	 * @param mismatch Score for a mismatch
	 * @return alignments Aligned strings
	 */
	public static String[] selfAlign(String mainRead, List<AssemblyEmbedded> embeddedReads, boolean reverse, int match, int openGap, int extGap, int mismatch)
	{
		StringBuilder[] alignmentsB = new StringBuilder[embeddedReads.size() + 1];
		AlignmentAffineGap aligner = new AlignmentAffineGap(match, openGap, extGap, mismatch);
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
			//so the start of the alignment is the end of the main read minus the beginning of the embedding - the length 
			//of the embedded read
			if(reverse)
			{
				start = mainRead.length() - start - embedded.length();
			}
			int end = start + embedded.length();
			String subMainRead = mainRead.substring(start, end);
			String[] alignedReads = aligner.getAlignment(subMainRead, embedded);
			if(i == 0)
			{
				alignmentsB[0] = new StringBuilder(mainRead);
			}
			alignmentsB[i + 1] = new StringBuilder();
			for(int j = 0; j < start; j++)
			{
				alignmentsB[i + 1].insert(j, '-');
			}
			alignmentsB[i + 1].append(alignedReads[1]);
			for(int j = end; j < mainRead.length(); j++)
			{
				alignmentsB[i + 1].insert(j, '-');
			}
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
		return alignments;
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
