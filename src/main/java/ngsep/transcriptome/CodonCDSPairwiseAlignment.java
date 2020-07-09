package ngsep.transcriptome;

public class CodonCDSPairwiseAlignment {
	public String[] pairwiseAlignment(String cds1, String cds2) {
		int proteinLengthCDS1=cds1.length()/3;
		int proteinLengthCDS2=cds2.length()/3;
		int[][] scores = new int[proteinLengthCDS1+1][proteinLengthCDS2+1];
		byte [][] direction = new byte [proteinLengthCDS1+1][proteinLengthCDS2+1];
		int indelPenalty = -2;
		int mismatchPenalty = -1;
		int matchScore = 1;
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
					else scores[i][j]+=mismatchPenalty;
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
		StringBuilder aln1 = new StringBuilder();
		StringBuilder aln2 = new StringBuilder();
		int maxI = scores.length-1;
		int maxJ = scores[maxI].length-1;
		int score = scores[maxI][maxJ];
		for(int i=maxI-1;i>=0.5*scores.length;i--) {
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
		}
		
		String gapCodon = "---";
		
		int i= maxI;
		int j= maxJ;
		while(i>0 || j>0) {
			byte d = direction[i][j];
			if(d==0) {
				addCodon(cds1.substring(3*(i-1), 3*i),aln1);
				addCodon(cds2.substring(3*(j-1), 3*j),aln2);
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
		
		String [] answer = new String[2];
		answer[0] = aln1.reverse().toString();
		answer[1] = aln2.reverse().toString();
		return answer;
	}

	private void addCodon(String codon, StringBuilder aln) {
		aln.append(codon.charAt(2));
		aln.append(codon.charAt(1));
		aln.append(codon.charAt(0));
	}
	
	
}
