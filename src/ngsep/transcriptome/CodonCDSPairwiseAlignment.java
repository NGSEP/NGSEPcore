package ngsep.transcriptome;

public class CodonCDSPairwiseAlignment {
	public String[] pairwiseAlignment(String cds1, String cds2) {
		int proteinLengthCDS1=cds1.length()/3;
		int proteinLengthCDS2=cds2.length()/3;
		int[][] scores = new int[proteinLengthCDS1+1][proteinLengthCDS2+1];
		byte [][] direction = new byte [proteinLengthCDS1+1][proteinLengthCDS2+1];
		int indelDistance = 2;
		int mismatchDistance = 1;
		for (int i = 0; i < scores.length; i++) {	
			for (int j = 0; j < scores[i].length; j++) {
				if(i==0 && j==0) {
					scores[i][j]=0;
					direction[i][j] = 0;
				} else if (i==0) {
					scores[i][j]=scores[i][j-1]+indelDistance;
					direction[i][j] = 1;
				} else if (j==0) {
					scores[i][j]=scores[i-1][j]+indelDistance;
					direction[i][j] = 2;
				} else {
					scores[i][j]=scores[i-1][j-1];
					direction[i][j] = 0;
					String codon1 = cds1.substring(3*(i-1), 3*i);
					String codon2 = cds2.substring(3*(j-1), 3*j);
					boolean equal = codon1.equals(codon2);
					if(!equal) scores[i][j]+=mismatchDistance;
					if(scores[i][j-1]+indelDistance < scores[i][j]) {
						scores[i][j] = scores[i][j-1]+indelDistance;
						direction[i][j] = 1;
					}
					if(scores[i-1][j]+indelDistance < scores[i][j]) {
						scores[i][j] = scores[i-1][j]+indelDistance;
						direction[i][j] = 2;
					}
				}
			}
		}
		StringBuilder aln1 = new StringBuilder();
		StringBuilder aln2 = new StringBuilder();
		int i = scores.length-1;
		int j = scores[i].length-1;
		String gapCodon = "---";
		
		
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
