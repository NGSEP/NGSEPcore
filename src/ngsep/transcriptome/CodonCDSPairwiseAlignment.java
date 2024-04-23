package ngsep.transcriptome;

import java.util.List;

import ngsep.alignments.PairwiseAligner;
import ngsep.alignments.PairwiseAlignment;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.io.FastaSequencesHandler;

public class CodonCDSPairwiseAlignment implements PairwiseAligner {
	private double pctIdentity;
	
	@Override
	public PairwiseAlignment calculateAlignment(CharSequence sequence1, CharSequence sequence2) {
		String cds1 = sequence1.toString();
		String cds2 = sequence2.toString();
		PairwiseAlignment alignment = new PairwiseAlignment(cds1, cds2);
		int proteinLengthCDS1=cds1.length()/3;
		int proteinLengthCDS2=cds2.length()/3;
		int[][] scores = new int[proteinLengthCDS1+1][proteinLengthCDS2+1];
		byte [][] direction = new byte [proteinLengthCDS1+1][proteinLengthCDS2+1];
		int indelPenalty = -4;
		int matchScore = 2;
		for (int i = 0; i < scores.length; i++) {	
			for (int j = 0; j < scores[i].length; j++) {
				if(i==0 && j==0) {
					scores[i][j]=0;
					direction[i][j] = 0;
				} else if (i==0) {
					//scores[i][j]=scores[i][j-1]+indelPenalty;
					scores[i][j]=0;
					direction[i][j] = 1;
				} else if (j==0) {
					//scores[i][j]=scores[i-1][j]+indelPenalty;
					scores[i][j]=0;
					direction[i][j] = 2;
				} else {
					scores[i][j]=scores[i-1][j-1];
					direction[i][j] = 0;
					String codon1 = cds1.substring(3*(i-1), 3*i);
					String codon2 = cds2.substring(3*(j-1), 3*j);
					boolean equal = codon1.equals(codon2);
					if(equal) scores[i][j]+=matchScore;
					else scores[i][j]+=calculateMismatchPenalty(codon1, codon2);
					if(scores[i][j-1]+indelPenalty > scores[i][j]) {
						scores[i][j] = scores[i][j-1]+indelPenalty;
						direction[i][j] = 1;
					}
					if(scores[i-1][j]+indelPenalty > scores[i][j]) {
						scores[i][j] = scores[i-1][j]+indelPenalty;
						direction[i][j] = 2;
					}
				}
			}
		}
		//printMatrix(scores);
		StringBuilder aln1 = new StringBuilder();
		StringBuilder aln2 = new StringBuilder();
		int maxI = scores.length-1;
		int maxJ = scores[maxI].length-1;
		int score = scores[maxI][maxJ];
		alignment.setScore(score);
		/*for(int i=maxI-1;i>=0.5*scores.length;i--) {
			if(scores[i][maxJ]>score) {
				maxI = i;
				score = scores[i][maxJ];
			}
		}
		for(int j=maxJ-1;j>=0.5*scores[maxI].length;j--) {
			if(scores[scores.length-1][j]>score) {
				maxI = scores.length-1;
				maxJ = j;
				score = scores[maxI][maxJ];
			}
		}*/
		String gapCodon = "---";
		int countIdentical = 0;
		int i= maxI;
		int j= maxJ;
		while(i>0 || j>0) {
			byte d = direction[i][j];
			if(d==0) {
				String codon1 = cds1.substring(3*(i-1), 3*i);
				String codon2 = cds2.substring(3*(j-1), 3*j);
				if(codon1.equals(codon2)) countIdentical+=3;
				addCodon(codon1,aln1);
				addCodon(codon2,aln2);
				i--;
				j--;
			} else if (d==1) {
				addCodon(gapCodon,aln1);
				addCodon(cds2.substring(3*(j-1), 3*j),aln2);
				j--;
			} else {
				addCodon(cds1.substring(3*(i-1), 3*i),aln1);
				addCodon(gapCodon,aln2);
				i--;
			}
		}
		
		String alignment1 = aln1.reverse().toString();
		String alignment2 = aln2.reverse().toString();
		if(alignment1.length()>0) {
			pctIdentity = 100.0*countIdentical;
			pctIdentity /= alignment1.length();
		} else pctIdentity = 0;
		alignment.setAlignedSequences(alignment1, alignment2);
		return alignment;
	}

	private void printMatrix(int[][] scores) {
		for(int i=0;i<scores.length;i++) {
			for(int j=0;j<scores[i].length;j++) {
				System.out.print("\t"+scores[i][j]);
			}
			System.out.println();
		}
		
	}

	private int calculateMismatchPenalty(String codon1, String codon2) {
		int penalty = 0;
		for(int i=0;i<codon1.length();i++) {
			if(codon1.charAt(i)!= codon2.charAt(i)) penalty--;
		}
		return penalty;
	}

	private void addCodon(String codon, StringBuilder aln) {
		aln.append(codon.charAt(2));
		aln.append(codon.charAt(1));
		aln.append(codon.charAt(0));
	}

	public double getPctIdentity() {
		return pctIdentity;
	}

	
	public static void main(String[] args) throws Exception {
		CodonCDSPairwiseAlignment aligner = new CodonCDSPairwiseAlignment();
		FastaSequencesHandler handler = new FastaSequencesHandler();
		List<QualifiedSequence> seqs = handler.loadSequences(args[0]);
		PairwiseAlignment aln = aligner.calculateAlignment(seqs.get(0).getCharacters(), seqs.get(1).getCharacters());
		System.out.println(aln.getAlignedSequence1());
		System.out.println(aln.getAlignedSequence2());
		
	}
	
	
}
