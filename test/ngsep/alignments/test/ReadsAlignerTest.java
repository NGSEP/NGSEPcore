package ngsep.alignments.test;

import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMValidationError;
import junit.framework.TestCase;
import ngsep.alignments.PairEndsAlignments;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.ReadsAligner;

public class ReadsAlignerTest extends TestCase {
	private ReadsAligner readsAligner;

	public void setUpReadsAligner() {
		readsAligner= new ReadsAligner();
	}

	public void testCheckPairEnd() {
		setUpReadsAligner();

		//Useful flags tool https://broadinstitute.github.io/picard/explain-flags.html

		//Case 1: ProperPair
		int flags1=ReadAlignment.FLAG_PAIRED;
		ReadAlignment aln1=new ReadAlignment("sequenceName", 5000, 5249, 250, flags1);
		int flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		ReadAlignment aln2=new ReadAlignment("sequenceName", 5500, 5251, 250, flags2);


		PairEndsAlignments properAlignment = readsAligner.checkPairEnd(aln1, aln2);
		assertEquals("Incorrect flags",properAlignment.getAln1().getFlags(),	ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_PROPER+ReadAlignment.FLAG_FIRST_OF_PAIR);

		//Case 2: Long deletion
		flags1=ReadAlignment.FLAG_PAIRED;
		aln1=new ReadAlignment("sequenceName", 21000, 21249, 250, flags1);
		flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		aln2=new ReadAlignment("sequenceName", 22000, 21751, 250, flags2);


		properAlignment = readsAligner.checkPairEnd(aln1, aln2);

		assertEquals ("Incorrect flags",properAlignment.getAln1().getFlags(),ReadAlignment.FLAG_PAIRED);

		//Case 3a: Same strand
		flags1=ReadAlignment.FLAG_PAIRED;
		aln1=new ReadAlignment("sequenceName", 9000, 9249, 250, flags1);
		flags2=ReadAlignment.FLAG_PAIRED;
		aln2=new ReadAlignment("sequenceName", 9500, 9251, 250, flags2);

		
		properAlignment = readsAligner.checkPairEnd(aln1, aln2);
		assertEquals ("Incorrect flags",properAlignment.getAln1().getFlags(),ReadAlignment.FLAG_PAIRED);

		//Case 3b: Same strand
		flags1=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		aln1=new ReadAlignment("sequenceName", 9000, 9249, 250, flags1);
		flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		aln2=new ReadAlignment("sequenceName", 9500, 9251, 250, flags2);


		properAlignment = readsAligner.checkPairEnd(aln1, aln2);
		assertEquals ("Incorrect flags",properAlignment.getAln1().getFlags(),ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND);

	}
}
