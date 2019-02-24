package ngsep.alignments.test;

import junit.framework.TestCase;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.ReadsAligner;

public class ReadsAlignerTest extends TestCase {
	private ReadsAligner readsAligner;

	public void setUpReadsAligner() {
		readsAligner= new ReadsAligner();
	}

	public void checkPairEndTest() {
		setUpReadsAligner();

		//Useful flags tool https://broadinstitute.github.io/picard/explain-flags.html

		//Case 1: ProperPair
		int flags1=ReadAlignment.FLAG_PAIRED;
		ReadAlignment aln1=new ReadAlignment("sequenceName", 5000, 5249, 250, flags1);
		int flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		ReadAlignment aln2=new ReadAlignment("sequenceName", 5500, 5251, 250, flags2);

		ReadAlignment properAlignment = readsAligner.checkPairEnd(aln1, aln2);
		assertEquals("Incorrect flags",properAlignment.getFlags(),	ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_PROPER);

		//Case 2: Long erase
		flags1=ReadAlignment.FLAG_PAIRED;
		aln1=new ReadAlignment("sequenceName", 21000, 21249, 250, flags1);
		flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		aln2=new ReadAlignment("sequenceName", 22000, 21751, 250, flags2);

		properAlignment = readsAligner.checkPairEnd(aln1, aln2);
		assertEquals ("Incorrect flags",properAlignment.getFlags(),ReadAlignment.FLAG_PAIRED);

		//Case 3a: Same strand
		flags1=ReadAlignment.FLAG_PAIRED;
		aln1=new ReadAlignment("sequenceName", 9000, 9249, 250, flags1);
		flags2=ReadAlignment.FLAG_PAIRED;
		aln2=new ReadAlignment("sequenceName", 9500, 9251, 250, flags2);

		properAlignment = readsAligner.checkPairEnd(aln1, aln2);
		assertEquals ("Incorrect flags",properAlignment.getFlags(),ReadAlignment.FLAG_PAIRED);

		//Case 3b: Same strand
		flags1=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		aln1=new ReadAlignment("sequenceName", 9000, 9249, 250, flags1);
		flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		aln2=new ReadAlignment("sequenceName", 9500, 9251, 250, flags2);

		properAlignment = readsAligner.checkPairEnd(aln1, aln2);
		assertEquals ("Incorrect flags",properAlignment.getFlags(),ReadAlignment.FLAG_PAIRED);

	}

}
