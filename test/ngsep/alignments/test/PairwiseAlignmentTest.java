package ngsep.alignments.test;

import java.io.IOException;
import java.util.List;
import java.util.Random;

import junit.framework.TestCase;
import ngsep.alignments.PairwiseAligner;
import ngsep.alignments.PairwiseAlignerAffineGap;
import ngsep.alignments.PairwiseAlignerDynamicKmers;
import ngsep.alignments.PairwiseAlignerSimpleGap;
import ngsep.alignments.PairwiseAlignerStaticBanded;
import ngsep.alignments.PairwiseAlignment;
import ngsep.sequences.DNASequence;
import ngsep.sequences.LimitedSequence;
import ngsep.sequences.QualifiedSequence;

import ngsep.sequences.io.*;

public class PairwiseAlignmentTest extends TestCase {
	private CharSequence seq1 = ""; 
	private CharSequence seq2 = ""; 
	private int match =1; 

    private int mismatch = 1; 

    private int  gapStart = 4;
    private int  gapExtension = 2;
	
	
	private int FindScore(CharSequence algseq1, CharSequence algseq2) {
		int cost = 0; 
		int status = 0;
		for (int i = 0; i<algseq1.length(); i++) {
			
			if (algseq1.charAt(i) == LimitedSequence.GAP_CHARACTER) {
				cost += insertionCost(status);
				status = 1;
				
			}else if (algseq2.charAt(i) == LimitedSequence.GAP_CHARACTER) {
				cost += deletionCost(status);
				status = 2;
			}else {
				cost+=matchMismatchCost(algseq1.charAt(i),algseq2.charAt(i));
				status = 0;
			}
		}
		return cost; 
		
		
	}
	
    private int insertionCost(int status){
    	if (status == 1) return -gapExtension;
        return -gapStart; 

    }
    private int deletionCost(int status){
    	if (status == 2) return -gapExtension;
        return -gapStart; 

    }
    private int matchMismatchCost(char l1, char l2){
        if (l1 == l2){
            return match;

        }else{
            return -mismatch; 
        }

    }
	
	
	public List<QualifiedSequence> loadSequences() throws IOException {
		String file = "./dataTest/exampleFLOGenes.fa"; 
		
		FastaSequencesHandler reader = new FastaSequencesHandler(); 
		reader.setSequenceType(DNASequence.class);
		reader.setKeepLowerCase(false);
		return reader.loadSequences(file);
	}
	
	//Create test with JUnit cases 
	public void testPairwiseAlignmentStaticBanded() throws IOException {
		
		//Create Aligners
		PairwiseAlignerStaticBanded aligner1 = new PairwiseAlignerStaticBanded(); 
		PairwiseAlignerSimpleGap aligner2 = new PairwiseAlignerSimpleGap();
				
		//Get random sequence
		Random rand = new Random();
		List<QualifiedSequence> sequences = loadSequences();
		double deltaMax = 10*0.01;
		double deltaMin = (0)*0.01;
		int cases = 50; 
		int smallError = 0;
		for (int c = 0; c<cases; c++) {
				
			
			int i = rand.nextInt(sequences.size()); 
			int j = rand.nextInt(sequences.size());
			//i=36;
			//j=44;
			seq1 = sequences.get(i).getCharacters(); 
			seq2 = sequences.get(j).getCharacters(); 
			//System.err.println("Testing alignments");
			//Banded Alignment
			//Set K-band value-Positive number
			
			System.out.println("Starting alignments. P1: "+i+" Name1: "+sequences.get(i).getName()+ " L1: "+seq1.length()+" P2: "+j+" Name 2: "+sequences.get(j).getName()+" L2: "+seq2.length());
			PairwiseAlignment aln1 = aligner1.calculateAlignment(seq1, seq2);
			System.out.println(aln1.getAlignedSequence1()+"\n"+aln1.getAlignedSequence2());
			
			//Score static banded

			int maxScoreAligner1 = FindScore(aln1.getAlignedSequence1(),aln1.getAlignedSequence2()); 

		
			//Simple Gap Alignment
			PairwiseAlignment aln2 = aligner2.calculateAlignment(seq1,seq2);
			System.out.println(aln2.getAlignedSequence1()+"\n"+aln2.getAlignedSequence2());
			
			//Score simple gap
			int maxScoreAligner2 = FindScore(aln2.getAlignedSequence1(),aln2.getAlignedSequence2()); 
			
			//System.out.println(maxScoreAligner1+" "+maxScoreAligner2);
			
			double delta = maxScoreAligner2- maxScoreAligner1 ; 
			System.out.println("L1: "+seq1.length()+" L2: "+seq2.length()+" score 1: "+maxScoreAligner1+" score 2: "+maxScoreAligner2+" Delta: "+ delta/Math.abs(maxScoreAligner2));
			
			assertTrue((delta/Math.abs(maxScoreAligner2))<0.1);

		}
		
		/*
		
		aln1 = aligner1.calculateAlignment(seq1, seq2);
		System.out.println(aln1[0]+"\n"+aln1[1]);
		
		aligner3.setK(seq2.length()-seq1.length());
		String [] aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		
		
		aln1 = aligner1.calculateAlignment(seq1, seq2);
		System.out.println(aln1[0]+"\n"+aln1[1]);
		aligner3.setK(seq1.length()-seq2.length());
		aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		
		System.err.println("Testing weird cases");
	
		
		aln1 = aligner1.calculateAlignment(seq1, seq2);
		System.out.println(aln1[0]+"\n"+aln1[1]);
		
		aligner3.setK(seq2.length()-seq1.length());
		aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		//assertNull(aln);
		
		
		aln1 = aligner1.calculateAlignment(seq1, seq2);
		System.out.println(aln1[0]+"\n"+aln1[1]);
		
		aligner3.setK(seq2.length()-seq1.length());
		aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		
		
		
		aln1 = aligner1.calculateAlignment(seq1, seq2);
		System.out.println(aln1[0]+"\n"+aln1[1]);
		
		aligner3.setK(seq2.length()-seq1.length());
		aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		
		//assertNull(aln);
		
		
		aln1 = aligner1.calculateAlignment(seq1, seq2);
		System.out.println(aln1[0]+"\n"+aln1[1]);
		
		aligner3.setK(seq1.length()-seq2.length());
		aln = aligner3.calculateAlignment(seq1, seq2);
		System.out.println(aln[0]+"\n"+aln[1]);
		*/
		
		//Old tests
		
		/*
		String seq1 = "ATGGACATCGTACACAGTATGGACATTGTACACAGT";
		String seq2 = "ATGGACAGTACATAGTATGGACAGTACATAGT";
		*/
		
		/*
		seq1 = "TACAGTTGGAATGACCACAAATCACAAATGTCAGATACCCGAAGCGCAGCG";
		seq2 =  "ACAGTTGGAATGACCACAAATCACAATGTCAGATACCCGAAGCGCAGC";
		*/
		
		/*
		seq1 = "GGATTGGGGGGGTGAAGAGTGGGGGGTGGGGGAGGGGGGGGGAGCTGGGGGGGGGAAA";
		seq2 = "GCGCAACGCGTCAGTGGGCTGATCATTAACTATCCGCTGGATGACAGGATGCCATTGCTGTGG";
		*/
		
		/*
		seq1 =  "ACTCTTTTTCGCCAAACGGATGCGCCCTC";
		seq2 = "TACTCTTTTTCGCCAAAACGGATGCGCCCTCT";
		*/
		
		/*
		seq1 = "TTTTTTTTCACTTTCTTTTTTTTTT";
		seq2 = "GAAGGCAAATCCATTCGAACCGTAATTCGTTACTGA";
		 */
		
		/*
		seq1 = "GGGAAGAAAACGATGAATGTGTG";
		seq2 = "AGAATGGGAAGAAAACGATGAATGTGTG";
		 */
		
		/*
		seq1 = "CACGGTTTATCCTGCATGCTCCGTCTTCT";
		seq2 = "GATTATCTGCATGCTCCGTCTTCT";
		 */
	}
	public void testPairwiseAlignmentDynamicKmers() throws IOException {
		
		//Create Aligners
		PairwiseAlignerDynamicKmers aligner1 = new PairwiseAlignerDynamicKmers(); 
		PairwiseAlignerSimpleGap aligner2 = new PairwiseAlignerSimpleGap();
				
		//Get random sequence
		Random rand = new Random();
		List<QualifiedSequence> sequences = loadSequences();
		double deltaMax = 10*0.01;
		double deltaMin = (0)*0.01;
		int cases = 1; 
		int smallError = 0;
		for (int c = 0; c<cases; c++) {
				
			
			int i = rand.nextInt(sequences.size()); 
			int j = rand.nextInt(sequences.size());
			//i=19;
			//j=35;
			seq1 = sequences.get(i).getCharacters(); 
			seq2 = sequences.get(j).getCharacters(); 
			//System.err.println("Testing alignments");
			//Banded Alignment
			//Set K-band value-Positive number
			
			System.out.println("Starting alignments. P1: "+i+" Name1: "+sequences.get(i).getName()+ " L1: "+seq1.length()+" P2: "+j+" Name 2: "+sequences.get(j).getName()+" L2: "+seq2.length());
			PairwiseAlignment aln1 = aligner1.calculateAlignment(seq1, seq2);
			System.out.println(aln1.getAlignedSequence1()+"\n"+aln1.getAlignedSequence2());
			
			//Score static banded

			int maxScoreAligner1 = FindScore(aln1.getAlignedSequence1(),aln1.getAlignedSequence2()); 

		
			//Simple Gap Alignment
			PairwiseAlignment aln2 = aligner2.calculateAlignment(seq1,seq2);
			System.out.println(aln2.getAlignedSequence1()+"\n"+aln2.getAlignedSequence2());
			
			//Score simple gap
			int maxScoreAligner2 = FindScore(aln2.getAlignedSequence1(),aln2.getAlignedSequence2());  
			
			//System.out.println(maxScoreAligner1+" "+maxScoreAligner2);
			
			double delta = maxScoreAligner2- maxScoreAligner1; 
			System.out.println("L1: "+seq1.length()+" L2: "+seq2.length()+" score 1: "+maxScoreAligner1+" score 2: "+maxScoreAligner2+" Delta: "+ delta/Math.abs(maxScoreAligner2));
			
			assertTrue((delta/Math.abs(maxScoreAligner2))<0.1);

		}
	}
	public void testRuntime() throws IOException {
		
		List<QualifiedSequence> sequences = loadSequences();
		System.out.println("Static banded");
		calculateRuntime(sequences, new PairwiseAlignerStaticBanded());
		System.out.println("Simple gap");
		calculateRuntime(sequences, new PairwiseAlignerSimpleGap());
		System.out.println("Affine gap");
		calculateRuntime(sequences, new PairwiseAlignerAffineGap(5000));
		System.out.println("Dynamic kmers");
		calculateRuntime(sequences, new PairwiseAlignerDynamicKmers());
		
		
	}

	public void calculateRuntime(List<QualifiedSequence> sequences, PairwiseAligner aligner) {
		int n = sequences.size();
		long s1 = System.currentTimeMillis();
		for(int i=0;i<20;i++) {
			String seq1 = sequences.get(i).getCharacters().toString();
			for(int j=i;j<20;j++) {
				String seq2 = sequences.get(j).getCharacters().toString();
				aligner.calculateAlignment(seq1, seq2);
			}
		}
		long s2 = System.currentTimeMillis();
		System.out.println("Time: "+(s2-s1));
	}
}
