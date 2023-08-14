package ngsep.alignments.test;

import junit.framework.TestCase;
import ngsep.alignments.PairwiseAlignerDynamicKmers;
import ngsep.alignments.PairwiseAlignerStaticBanded;
import ngsep.sequences.SimpleEditDistanceMeasure;

public class PairwiseAlignmentTestBanded extends TestCase {
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
		//assertEquals(aln1[0], aln2[0]);
		//assertEquals(aln1[1], aln2[1]);
		PairwiseAlignerStaticBanded aligner3 = new PairwiseAlignerStaticBanded(); 
		String [] aln3 = aligner3.calculateAlignment(seq1, seq2); 
		System.out.println(aln3[0]+"\n"+aln3[1]);


		
		System.err.println("Testing alignments");
		seq1 = "GGGAAGAAAACGATGAATGTGTG";
		seq2 = "AGAATGGGAAGAAAACGATGAATGTGTG";
		String [] aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		
		seq1 = "CACGGTTTATCCTGCATGCTCCGTCTTCT";
		seq2 = "GATTATCTGCATGCTCCGTCTTCT";
		aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		
		System.err.println("Testing weird cases");
		seq1 = "TTTTTTTTCACTTTCTTTTTTTTTT";
		seq2 = "GAAGGCAAATCCATTCGAACCGTAATTCGTTACTGA";
		aln = aligner3.calculateAlignment(seq1, seq2);
		//System.out.println(aln[0]+"\n"+aln[1]);
		assertNull(aln);
		seq1 =  "ACTCTTTTTCGCCAAACGGATGCGCCCTC";
		seq2 = "TACTCTTTTTCGCCAAAACGGATGCGCCCTCT";
		aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		seq1 = "GGATTGGGGGGGTGAAGAGTGGGGGGTGGGGGAGGGGGGGGGAGCTGGGGGGGGGAAA";
		seq2 = "GCGCAACGCGTCAGTGGGCTGATCATTAACTATCCGCTGGATGACAGGATGCCATTGCTGTGG";
		aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]+" "+seq1.length());
		//assertNull(aln);
		seq1 = "TACAGTTGGAATGACCACAAATCACAAATGTCAGATACCCGAAGCGCAGCG";
		seq2 =  "ACAGTTGGAATGACCACAAATCACAATGTCAGATACCCGAAGCGCAGC";
		aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
	}
}
