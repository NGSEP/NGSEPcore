package ngsep.alignments.test;

import junit.framework.TestCase;
import ngsep.alignments.PairwiseAlignerDynamicKmers;
import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.sequences.SimpleEditDistanceMeasure;

public class PairwiseAlignmentTest extends TestCase {
	public void testPairwiseAlignmentKmers() {
		String seq1 = "ATGGACATCGTACACAGTATGGACATTGTACACAGT";
		String seq2 = "ATGGACAGTACATAGTATGGACAGTACATAGT";
		SimpleEditDistanceMeasure aligner1 = new SimpleEditDistanceMeasure();
		String [] aln1 = aligner1.calculateAlignment(seq1, seq2);
		PairwiseAlignerDynamicKmers aligner2 = new PairwiseAlignerDynamicKmers();
		//PairwiseAlignerSimpleGap aligner2 = new PairwiseAlignerSimpleGap();
		String [] aln2 = aligner2.calculateAlignment(seq1, seq2);
		System.out.println(aln1[0]+"\n"+aln1[1]);
		System.out.println(aln2[0]+"\n"+aln2[1]);
		assertEquals(aln1[0], aln2[0]);
		assertEquals(aln1[1], aln2[1]);
		
		System.err.println("Testing alignments");
		seq1 = "GGGAAGAAAACGATGAATGTGTG";
		seq2 = "AGAATGGGAAGAAAACGATGAATGTGTG";
		String [] aln = aligner2.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		System.err.println("Testing weird cases");
		seq1 = "TTTTTTTTCACTTTCTTTTTTTTTT";
		seq2 = "GAAGGCAAATCCATTCGAACCGTAATTCGTTACTGA";
		aln = aligner2.calculateAlignment(seq1, seq2);
		//System.out.println(aln[0]+"\n"+aln[1]);
		assertNull(aln);
		seq1 =  "ACTCTTTTTCGCCAAACGGATGCGCCCTC";
		seq2 = "TACTCTTTTTCGCCAAAACGGATGCGCCCTCT";
		aln = aligner2.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		seq1 = "GGATTGGGGGGGTGAAGAGTGGGGGGTGGGGGAGGGGGGGGGAGCTGGGGGGGGGAAA";
		seq2 = "GCGCAACGCGTCAGTGGGCTGATCATTAACTATCCGCTGGATGACAGGATGCCATTGCTGTGG";
		aln = aligner2.calculateAlignment(seq1, seq2);
		//System.out.println(aln[0]+"\n"+aln[1]);
		assertNull(aln);
		seq1 = "TACAGTTGGAATGACCACAAATCACAAATGTCAGATACCCGAAGCGCAGCG";
		seq2 =  "ACAGTTGGAATGACCACAAATCACAATGTCAGATACCCGAAGCGCAGC";
		aln = aligner2.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
	}
}
