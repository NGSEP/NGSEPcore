package ngsep.alignments.test;

import java.io.IOException;

import junit.framework.TestCase;
import ngsep.alignments.ReadAlignment;
import ngsep.alignments.ReadsAligner;

public class ReadsAlignerTest extends TestCase {
	public void testCheckPairEnd() throws IOException {
		ReadsAligner readsAligner = new ReadsAligner();
		//Useful flags tool https://broadinstitute.github.io/picard/explain-flags.html
		//Case 1: ProperPair
		int flags1=ReadAlignment.FLAG_PAIRED;
		ReadAlignment aln1=new ReadAlignment("sequenceName", 5000, 5249, 250, flags1);
		int flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		ReadAlignment aln2=new ReadAlignment("sequenceName", 5500, 5251, 250, flags2);
		aln2.setMateNegativeStrand(true);

		boolean proper =readsAligner.isValidPair(aln1, aln2,true);
		assertEquals(true, proper);
		proper =readsAligner.isValidPair(aln2, aln1,true);
		assertEquals(true, proper);


		//Case 2: Long deletion
		flags1=ReadAlignment.FLAG_PAIRED;
		aln1=new ReadAlignment("sequenceName", 21000, 21249, 250, flags1);
		flags2=ReadAlignment.FLAG_PAIRED+ReadAlignment.FLAG_READ_REVERSE_STRAND;
		aln2=new ReadAlignment("sequenceName", 22000, 21751, 250, flags2);
		aln2.setMateNegativeStrand(true);


		boolean unProper =readsAligner.isValidPair(aln1, aln2,true);
		assertEquals(false, unProper);
		unProper =readsAligner.isValidPair(aln2, aln1,true);
		assertEquals(false, unProper);

		unProper =readsAligner.isValidPair(aln1, aln2,false);
		assertEquals(true, unProper);
		unProper =readsAligner.isValidPair(aln2, aln1,false);
		assertEquals(true, unProper);
	}
}
